// stc_engine.h -- the permutation engine of STcompare's Viladomat null (dev/engine-spec.md, sections 2,
// 3 and 6): legacy-exact surrogates, null correlations and delta* for many genes at once, on worker
// threads, with optional Besag-Clifford (adaptive) stopping.
//
// Plain C++: no R API anywhere in this file or in stc_engine.cpp, so every function can run on a worker
// thread. The host (src/stc_engine_rcpp.cpp) copies the R inputs into the engine on R's main thread and
// passes a callback that polls for user interrupts; the engine calls that callback only from the
// thread that called it.
//
// Data flow (one Engine per call of the R glue, R/engine.R):
//   session   coordinates (x1 = long = pos[, 2], x2 = lat = pos[, 1], locfit's order), the variogram
//             subsample `ids`, its geoR pair table (built in R by .stc_variog_plan()), and the smoother
//             operators of every distinct delta of every task (built once, in parallel over deltas).
//   pool      the data vectors (genes) as an N x V matrix; tasks and their targets refer to its columns.
//   task      a source vector permuted and smoothed with a delta grid (indices into the operators),
//             correlated with one or more target vectors. One task gives, for every permutation b,
//             delta* (a grid position), the surrogate (optional) and one null correlation per target.
//   unit      the tasks that stop together: in gene-wise comparisons the two directions of a gene
//             (permute X, correlate with Y; permute Y, correlate with X), which run in lockstep on the
//             same permutations. A failed task fails its unit (legacy NA row).
//   batch     permutations b_from..b_to (1-based) for a set of active units: the permutation indices
//             (generated in R from the legacy stream) and the legacy noise (generated here, once per
//             batch and shared by all tasks, or supplied by the host) are read; the work items
//             (task x sub-chunk of permutations) are scheduled dynamically over std::thread workers; then
//             the calling thread scans the new nulls of each unit in permutation order for the stopping
//             rule and retires stopped units.
//
// Every (task, permutation) result is computed by a fixed sequence of operations that does not depend
// on which other permutations share its work item, so results do not depend on the number of threads,
// the sub-chunk size or the batch schedule.
#ifndef STC_ENGINE_H
#define STC_ENGINE_H

#include <atomic>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "stc_smoother.h"
#include "stc_stats.h"
#include "stc_variog.h"
#include "stc_fp.h"

namespace stc {

// Why a task, a direction or a unit has no result (the legacy code returns an NA row for all of them,
// except where noted).
enum EngineStatus {
  ENG_OK = 0,
  ENG_SOURCE_NONFINITE = 1,     // NA/NaN/Inf in the permuted vector (geoR::variog fails on it)
  ENG_TARGET_NONFINITE = 2,     // NA/NaN/Inf in a vector the surrogates are correlated with
  ENG_SOURCE_CONSTANT = 3,      // the permuted vector has zero variance (every lm() slope is NA)
  ENG_TARGET_CONSTANT = 4,      // a correlated vector has zero variance (cor() is NA)
  ENG_PLAN_FAILED = 5,          // no variogram pair table (geoR::variog fails), or fewer than 2 bins
                                // (lm() gives an NA slope)
  ENG_SMOOTHER_FAILED = 6,      // a delta of the grid has no smoother: N * delta < 2 (locfit fails or
                                // kills R), out of vertex space, or delta <= 0
  ENG_OLS_DEGENERATE = 7,       // lm() would report an NA slope at some (permutation, delta); the
                                // legacy code then fails in geoR::variog
  ENG_SURROGATE_NONFINITE = 8,  // a rescaled field has a non-finite value on the subsample
                                // (geoR::variog fails)
  ENG_RSS_NAN = 9,              // every RSS of a permutation is NaN (which.min() gives integer(0))
  ENG_COR_NA = 10               // a null correlation is NA or NaN (the legacy p-value is NA)
};

const char* engine_status_text(int status);

enum UnitState {
  UNIT_ACTIVE = 0,      // more permutations can be added
  UNIT_STOPPED_H = 1,   // a direction reached h exceedances at permutation len (Besag-Clifford)
  UNIT_STOPPED_CAP = 2, // len reached n_max (fixed-B runs always end here)
  UNIT_FAILED = 3       // a task or direction failed at permutation fail_b (0: before any permutation)
};

// Run fn(item, worker) for item = 0..n_items-1 on up to n_threads std::thread workers (worker =
// 0..n_threads-1; items are taken from an atomic counter, so the schedule is dynamic). The calling
// thread waits and calls interrupted() about every 100 ms; when it returns true the workers finish
// their current item and take no more. The first exception thrown by fn is captured, the remaining
// items are abandoned, and it is rethrown on the calling thread after every worker has joined.
// Returns false if interrupted. fn must not call the R API.
bool parallel_for(int n_items, int n_threads, const std::function<void(int, int)>& fn,
                  const std::function<bool()>& interrupted);

struct EngineOptions {
  int n_threads = 1;
  int chunk = 16;               // permutations per work item
  bool keep_surrogates = false; // store the surrogate fields (returnPermutations)
  int cor_mode = COR_PLAIN;     // see stc_stats.h; .stc_cor_mode() picks it on R's main thread
  double seed = 0.0;            // legacy noise of permutation b: set.seed((int)(seed + b)) under
                                // L'Ecuyer-CMRG, then one rnorm(N) per grid position
  bool noise_supplied = false;  // the host supplies the noise of every batch instead
};

// One distinct delta value and its smoother.
struct EngineDelta {
  double delta = 0.0;
  int status = SMOOTH_BAD_INPUT;  // SmoothStatus
  int k = 0, nv = 0, nvm = 0, depth = 0;
  SmoothOperator op;              // m = op.m real vertices
};

struct EngineTask {
  int source = 0;                 // pool column that is permuted
  std::vector<int> grid;          // delta indices (into Engine::deltas); position k = noise block k
  std::vector<int> targets;       // pool columns the surrogates are correlated with
  std::vector<double> rabs;       // |observed r| per target (exceedances: |null| >= rabs)
  int unit = 0;
  std::vector<double> tvar;       // target variogram: the source on ids
  double cost = 0.0;              // relative cost of one permutation (scheduling only)
  int status = ENG_OK;            // pre-check (before any permutation)
  std::string message;
  std::vector<int> target_status; // pre-check of each target
  // Run-time failure with the smallest permutation index (1-based) and, within it, grid position; the
  // permutation index is also kept in Engine::fail_b_ (atomic) for the workers' skip test.
  int fail_k = -1;
  int fail_status = ENG_OK;
  std::vector<int> dir_fail_b;    // per target: the smallest permutation whose null is NA (INT_MAX)
  // Outputs for permutations 1..cap (slot b - 1).
  int cap = 0;
  std::vector<int> dstar;         // grid position of delta*, -1 where the permutation failed
  std::vector<double> nulls;      // cap x T, permutation-major: nulls[(b - 1) * T + t]
  std::vector<double> surr;       // cap x N, surrogate b at (b - 1) * N (keep_surrogates only)
};

struct EngineUnit {
  std::vector<int> tasks;
  std::vector<int> dir_task, dir_target;  // directions: (task, position in its targets)
  bool fail_on_direction = true;  // an NA null correlation fails the unit (gene-wise mode); in the
                                  // within-sample mode it only fails that pair, in R
  int len = 0;                    // permutations kept (L once stopped)
  std::vector<int> counts;        // exceedances of each direction in permutations 1..len
  int state = UNIT_ACTIVE;
  int status = ENG_OK;
  std::string message;
  int fail_b = 0;                 // permutation of the failure (0: pre-check)
  int fail_task = -1;             // failing task, or -1
  int fail_dir = -1;              // failing direction, or -1
};

class Engine {
 public:
  static const std::uint32_t kMagic = 0x53544345u;
  std::uint32_t magic = kMagic;

  // ---- session ----
  int n = 0;                       // pixels
  std::vector<double> x1, x2;      // long, lat
  std::vector<int> ids;            // variogram subsample, 0-based
  VariogTable vt;                  // geoR pair table on ids (ptr, pi, pj, n, nbins)
  bool plan_ok = false;
  std::string plan_reason;
  std::vector<EngineDelta> deltas;
  EngineOptions opt;

  // ---- tasks ----
  bool tasks_defined = false;
  int pool_cols = 0;
  std::vector<double> pool;        // n x pool_cols, column-major
  std::vector<CorTarget> pool_cor; // correlation target data of every pool column
  std::vector<int> pool_status;    // ENG_OK, ENG_SOURCE_NONFINITE or ENG_SOURCE_CONSTANT
  std::vector<EngineTask> tasks;
  std::vector<EngineUnit> units;
  int kmax = 0;                    // longest grid

  // Session: copies the coordinates, the subsample and the pair table, then builds the smoother of
  // every delta (in parallel over deltas). Returns false if interrupted.
  bool build(const double* x1_, const double* x2_, int n_, const int* ids_, int n_ids,
             const VariogTable& table, bool plan_ok_, const std::string& plan_reason_,
             const double* delta_values, int n_deltas, const EngineOptions& options,
             const std::function<bool()>& interrupted);

  // Tasks and units (once per session). pool: n x V column-major. For task t: task_source[t],
  // task_grid[t] (delta indices), task_targets[t] (pool columns), task_rabs[t] (one per target),
  // task_unit[t]. unit_fail_on_dir: one per unit (units are numbered 0..n_units-1 and each needs a
  // task). Runs the pre-checks and the target variograms; units with a failed pre-check are FAILED.
  void define(const double* pool_, int V, const std::vector<int>& task_source,
              const std::vector<std::vector<int>>& task_grid,
              const std::vector<std::vector<int>>& task_targets,
              const std::vector<std::vector<double>>& task_rabs,
              const std::vector<int>& task_unit, const std::vector<int>& unit_fail_on_dir);

  // One batch: permutations b_from..b_to (1-based, inclusive) for the listed units, each of which must
  // have len == b_from - 1. perm: n x nb permutation indices, 0-based and validated by the caller
  // (nb = b_to - b_from + 1). noise: n x kn x nb (column-major) when opt.noise_supplied, with
  // kn >= the longest grid of the listed units; otherwise nullptr (generated here). A unit is run
  // only if it is not failed, every count is below h, and len < n_max; after the batch each unit's
  // new nulls are scanned in permutation order: a failure at b fails the unit, the first b at which a
  // direction reaches h stops it (len = b, outputs past b are discarded), and len = n_max caps it.
  // Requires b_to <= n_max. Returns false if interrupted (the units are then unchanged).
  bool run_batch(const std::vector<int>& unit_list, int b_from, int b_to, const int* perm,
                 const double* noise, int kn, double h, int n_max,
                 const std::function<bool()>& interrupted);

  void set_threads(int n_threads, int chunk);

 private:
  struct Workspace {
    std::vector<double> Y;        // n x C, gathered and centred columns
    std::vector<double> centre;   // C
    std::vector<double> Z;        // sum over the grid of m x C: Wn (y - c) of every grid position
    std::vector<std::size_t> zoff;// offset of grid position k in Z
    std::vector<double> Xd;       // n_ids x C: smoothed values on the subsample
    std::vector<double> H;        // n_ids x C: rescaled values on the subsample
    std::vector<double> G, G2;    // nbins x C: variograms of Xd and H
    std::vector<double> s0, s1;   // kmax x C: sqrt(|intercept|), sqrt(|slope|)
    std::vector<double> rss;      // kmax x C
    std::vector<int> fail_k, fail_st;  // C: first failing grid position and status, or -1
    std::vector<double> surr;     // n: one surrogate (when they are not kept)
    std::vector<double> r;        // T: null correlations of one surrogate
    std::vector<int> rst;         // T: their statuses
    std::vector<const double*> yptr;     // T: target vectors
    std::vector<const CorTarget*> tptr;  // T: target data
  };

  void ensure_capacity(EngineTask& t, int cap);
  void size_workspaces(const std::vector<int>& task_list, int n_workers);
  void process_item(int task_index, int p0, int p1, int b_from, Workspace& w);
  void record_task_failure(int task_index, int b, int k, int status);
  void scan_unit(EngineUnit& u, int b_from, int b_to, double h, int n_max);
  void fail_unit(EngineUnit& u, int b, int task, int dir, int status, const std::string& msg);

  std::vector<Workspace> ws_;
  std::unique_ptr<std::atomic<int>[]> fail_b_;  // per task: smallest failing permutation (INT_MAX)
  std::mutex fail_mutex_;
  // batch inputs (valid during run_batch)
  const int* perm_ = nullptr;
  const double* noise_ = nullptr;
  int noise_k_ = 0;               // noise blocks per permutation
  std::vector<double> noise_buf_;
};

}  // namespace stc

#endif  // STC_ENGINE_H
