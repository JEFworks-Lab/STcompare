// stc_fp.h -- no floating-point contraction (fused multiply-add) from here to the end of the
// translation unit, and a compile error under -ffast-math and similar flags (see the end).
//
// Include this header AFTER all system and Rcpp headers, before any floating-point code of
// STcompare. Every STcompare header does so, so a source file that includes one of them is covered
// from that point on.
//
// Why: the engine reproduces locfit's adaptive tree (distances and vertex bandwidths), geoR's
// variogram sums and R's own arithmetic exactly. A fused multiply-add rounds once where the
// reference code rounds twice. That changes last bits, and through a vertex bandwidth or a pair
// distance it can change a tree or a bin count. Compilers contract by default:
//   - clang (macOS, and LLVM toolchains elsewhere): -ffp-contract=on, within an expression;
//   - GCC: -ffp-contract=fast, even across statements, on targets with FMA instructions
//     (arm64; x86-64 only with -mfma or -march=native).
// The installed locfit and geoR binaries on macOS arm64 contain no FMA instructions, so the
// replicas below must not contain any either. Where R itself was compiled with contraction (R's
// cor() on macOS arm64), the code calls std::fma() explicitly; explicit calls are not affected by
// these pragmas.
//
// Why pragmas and not "PKG_CXXFLAGS = -ffp-contract=off" in src/Makevars: R CMD check reports every
// -f... flag in PKG_CXXFLAGS as a non-portable flag, which is a WARNING (tools:::.check_make_vars).
// R CMD check does not inspect these pragmas (tools:::.check_pragmas only looks for pragmas that
// suppress diagnostics).
//
// GCC: '#pragma GCC optimize' gives every function defined after it the option. Standard-library
// inline functions defined in headers included before it are still inlined into our functions
// (verified with GCC 13), which is why the system headers must come first. clang: the pragma applies
// to every expression that follows it, including inline functions defined in later headers.
//
// Value-changing optimisations (-ffast-math, -Ofast, -ffinite-math-only, -fassociative-math,
// -freciprocal-math, -funsafe-math-optimizations, clang's -ffp-model=fast) are refused with a
// compile error. They are sometimes put in ~/.R/Makevars for speed, and they would silently break
// this code: reassociation changes summation orders, so the results are no longer exact, and
// finite-math-only lets the compiler delete the NaN checks, so missing values give numbers (cor = -1)
// and wrong statuses instead of NA (measured with -O2 -ffast-math: 168 of 1083 component-test
// expectations failed). A pragma cannot undo them reliably: clang has none for finite-math-only, and
// GCC applies "#pragma GCC optimize" only to functions defined after it, not to the standard
// library's. The compilers announce these modes with the macros tested here (GCC defines
// __ASSOCIATIVE_MATH__ and __RECIPROCAL_MATH__ for the partial flags; clang defines no macro for
// them, so only its -ffast-math, -Ofast, -ffp-model=fast and -ffinite-math-only are caught).
#ifndef STC_FP_H
#define STC_FP_H

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__) || \
    defined(__ASSOCIATIVE_MATH__) || defined(__RECIPROCAL_MATH__)
#error "STcompare must not be compiled with -ffast-math, -Ofast, -ffinite-math-only or other unsafe floating-point optimisations: its compiled code reproduces R, locfit and geoR exactly and needs IEEE NaN semantics. Remove these flags (for example from CXXFLAGS or CXX17FLAGS in ~/.R/Makevars) and reinstall."
#endif

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#endif  // STC_FP_H
