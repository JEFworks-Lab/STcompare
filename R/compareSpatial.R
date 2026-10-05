#' Compare the spatial expression patterns of genes between two samples
#'
#' @description For every gene measured in two samples on one shared grid of
#'   pixels, \code{compareSpatial()} tests whether the gene's expression
#'   patterns are correlated beyond what the spatial autocorrelation of each
#'   pattern alone would produce (spatial correlation test), and how often the
#'   two samples' values at the same pixel are within a fold change of each
#'   other (spatial similarity). It returns one row per gene.
#'
#' @details \strong{Input.} The two samples must share their pixels: rasterize
#'   them together, in one call such as
#'   \code{SEraster::rasterizeGeneExpression(list(x, y), ...)}, so that a pixel
#'   name means the same location in both. The pixels used are those whose
#'   names (column names) are in both objects, in the column order of
#'   \code{x}, and their coordinates (\code{SpatialExperiment::spatialCoords()})
#'   must agree.
#'
#'   \strong{Spatial correlation test.} The observed statistic is Pearson's
#'   correlation \code{r} of the gene's values in \code{x} and \code{y} over
#'   the shared pixels. Its null distribution comes from spatially
#'   autocorrelated permutations (Viladomat et al. 2014): the values of
#'   \code{x} are shuffled across the pixels, smoothed with a Gaussian kernel
#'   whose window holds a share \eqn{\delta} of the pixels, and rescaled with
#'   added noise so that their variogram matches the variogram of \code{x}; the
#'   \eqn{\delta} of \code{delta} whose variogram matches best (delta star) is
#'   kept, the values of this surrogate are replaced by the gene's own values
#'   in the surrogate's rank order (\code{surrogate = "remap"}, the default;
#'   see \emph{Rank-remapped surrogates}), and its correlation with \code{y} is
#'   one null correlation. The same is done with \code{y} permuted and
#'   correlated with \code{x}. If there are more than 1000 shared pixels, the variograms use a
#'   random subsample of 1000 of them, drawn from \code{seed} with R's default
#'   random number generator whatever \code{RNGkind()} is set to. The smoother
#'   and the variogram are those of the other correlation functions of the
#'   package (\code{\link{spatialCorrelationGeneExp}()}).
#'
#'   \strong{Adaptive p-values.} The permutations of a gene stop as soon as
#'   its p-value is clearly not small (Besag and Clifford 1991). Both
#'   directions (\code{x} permuted, \code{y} permuted) run in lockstep, and a
#'   gene stops at the first permutation \eqn{L} at which either direction
#'   has \code{exceedances} null correlations whose absolute value is at least
#'   \eqn{|r|}, or at \code{nPermutations}. Its p-value is
#'   \code{exceedances / L} if it stopped early, and otherwise
#'   \eqn{(\max(b_x, b_y) + 1) / (L + 1)}{(max(bX, bY) + 1) / (L + 1)}, where
#'   \eqn{b_x}{bX} and \eqn{b_y}{bY} count the exceedances of each direction.
#'   This equals the larger of the two directions' own sequential p-values, so
#'   it needs no assumption about how the two directions are related: if the
#'   observed correlation is distributed like the null correlations of either
#'   direction, \eqn{P(p \le \alpha) \le \alpha}{P(p <= alpha) <= alpha} at
#'   every \eqn{\alpha}. The surrogates approximate that null distribution;
#'   with \code{surrogate = "gaussian"} the approximation fails for genes
#'   detected in few pixels (see \emph{Rarely detected genes}). A
#'   gene that is clearly not significant stops after a few dozen
#'   permutations, and a significant gene runs to \code{nPermutations}, which
#'   sets the smallest possible p-value, \code{1 / (nPermutations + 1)}. With
#'   \code{exceedances = Inf} every gene gets exactly \code{nPermutations}
#'   permutations and \code{p = max(pX, pY)}. The p-values are adjusted
#'   across the tested genes with \code{adjustMethod}.
#'
#'   \strong{Rank-remapped surrogates.} With \code{surrogate = "remap"} (the
#'   default), the values of each surrogate are replaced by the gene's own
#'   values, placed in the rank order of the surrogate (the amplitude
#'   adjustment of the AAFT surrogates of Theiler et al. 1992, in a single
#'   step): the surrogate keeps the spatial arrangement of the smoothed,
#'   rescaled permutation and has exactly the distribution of values of the
#'   gene, zeros included. Delta star is chosen before the remapping, and for a
#'   gene with a spatial pattern the remapping is a monotone distortion of the
#'   surrogate's amplitudes, which changes its variogram a little. For a gene
#'   without spatial structure the ranks of its surrogates are close to
#'   uniformly random permutations, so its remapped surrogates are plain
#'   random permutations of its values and the test is the exact permutation
#'   test, whatever the gene's distribution (exactly so when the other gene has
#'   no spatial structure either). Many rearrangements of a gene detected in
#'   few pixels give exactly the observed correlation; as in the exact test,
#'   such ties count as exceedances (a null correlation within a relative 1e-9
#'   of \eqn{|r|} counts, because \code{cor()} rounds each rearrangement
#'   differently). In simulations on the kidney and brain grids
#'   (\code{bench/calibration-results.md} in the source repository), remapped
#'   surrogates kept \eqn{P(p \le \alpha) \le \alpha}{P(p <= alpha) <= alpha}
#'   for independent genes detected in 1 to 100 percent of the pixels, had
#'   the same power as gaussian surrogates on correlated Gaussian fields,
#'   recovered the published kidney and brain genes at least as well, and
#'   cost about 7 percent more per permutation.
#'   \code{surrogate = "gaussian"} uses the surrogates as the smoothing and
#'   the added noise leave them. They are the surrogates of the other
#'   correlation functions of the package and of the published analyses, and
#'   are kept for comparability with them: their values are close to normally
#'   distributed, which gives p-values that are far too small for genes
#'   detected in few pixels (next paragraph) and conservative ones for skewed
#'   genes with a spatial pattern.
#'
#'   \strong{Rarely detected genes.} Gaussian surrogates are smoothed and mixed
#'   with Gaussian noise, so their values are close to normally distributed.
#'   A gene detected in few pixels has a few high values among many zeros, and
#'   the correlation of two such genes is dominated by the few pixels where
#'   both are detected: chance coincidences give extreme correlations more
#'   often than gaussian surrogates do, so small p-values are too small, ten
#'   times or more for genes detected in one or two pixels. With
#'   \code{surrogate = "gaussian"}, genes detected in fewer than a share
#'   \code{minDetected} of the shared pixels of either sample are therefore
#'   not tested for correlation (\code{status = "skipped"}). A gene is
#'   detected at a pixel where its value is above its lowest value over the
#'   shared pixels (for counts and normalized counts: where it is not zero).
#'   With gaussian surrogates the default asks for \eqn{\sqrt{N}}{sqrt(N)} of
#'   the \eqn{N} shared pixels (18 of 311, 47 of 2170): the departure from the
#'   surrogates grows roughly as \eqn{N / (k_x k_y)}{N / (kx * ky)} for genes
#'   detected in \eqn{k_x}{kx} and \eqn{k_y}{ky} pixels, and this removes the
#'   genes for which it is largest; sparse or zero-inflated genes above it can
#'   still get p-values that are somewhat too small, so with gaussian
#'   surrogates test the genes with a spatial pattern in both samples, or
#'   raise \code{minDetected}. Remapped surrogates have the gene's own values,
#'   and the problem disappears: in the simulations their p-values were
#'   calibrated down to genes detected in a single pixel, so with
#'   \code{surrogate = "remap"} the default \code{minDetected} is 0 and every
#'   gene that is not constant in a sample is tested. Such genes cost little
#'   (a gene without spatial structure stops after about ten permutations),
#'   and a gene detected in a single pixel cannot get a p-value below
#'   \eqn{1 / N}: the exact test is discrete.
#'
#'   \strong{Spatial similarity.} The same computation as
#'   \code{\link{spatialSimilarity}()}: for each gene, pixels are kept when
#'   their value is above the gene's \code{minQuantile} quantile in \code{x} or
#'   in \code{y}; zeros among them are replaced by 1e-4; \code{similarity} is
#'   the share of kept pixels with
#'   \eqn{|\log_2(y / x)| \le}{|log2(y / x)| <=} \code{foldChange}. It is
#'   \code{NA} when fewer than a share \code{minPixels} of the pixels are kept.
#'   A fold change needs values that are not negative (counts or normalized
#'   expression on a linear scale): for a gene with negative values the
#'   similarity is \code{NA}, and a warning lists such genes.
#'
#'   \strong{Reproducibility.} The permutations and the noise of each gene and
#'   direction come from their own random stream, which depends only on
#'   \code{seed}, the gene's name, the direction and the permutation number.
#'   The results therefore do not depend on \code{nThreads}, on the order of
#'   the genes, or on which other genes are in the call (except \code{padj}),
#'   and the global random number generator (\code{.Random.seed}) is left
#'   unchanged. The permutations are drawn over the shared pixels in the
#'   column order of \code{x} (the order of \code{y} does not matter), and
#'   each direction has its own streams: reordering the pixels of \code{x}, or
#'   swapping \code{x} and \code{y}, gives other permutations and p-values
#'   that differ within Monte Carlo error. Results are deterministic on a
#'   given platform; the normal deviates use the C library's \code{log()} and
#'   can differ in the last bits between platforms.
#'
#'   \strong{Run time.} It is proportional to the total number of permutations
#'   (the sum of the \code{nPermutations} column), and the work of all genes is
#'   shared by \code{nThreads} threads. With the defaults, genes that are
#'   clearly not significant use a few dozen permutations each and significant
#'   genes 10000, so the run time is dominated by the significant genes. A
#'   smaller \code{nPermutations}, such as 1000, was 6 to 7 times faster on the
#'   published data and still resolves p-values down to 1 / 1001.
#'
#'   \strong{Skipped and failed genes.} A gene that is constant in a sample
#'   (for example never detected) or detected in fewer pixels than
#'   \code{minDetected} asks for (by default only with
#'   \code{surrogate = "gaussian"}) is not tested for correlation: it gets
#'   \code{status = "skipped"}, \code{NA} p-values and the reason in
#'   \code{message} (its similarity is still computed), and one message counts
#'   such genes. A gene that cannot be tested gets \code{status = "failed"},
#'   \code{NA} in the statistics that could not be computed, and the reason in
#'   \code{message}; one warning lists the failed genes, and the other genes
#'   are not affected. This happens for genes with missing or infinite values
#'   on the shared pixels (every statistic is \code{NA}) and when the
#'   permutations of a gene fail (for example a degenerate variogram fit).
#'
#' @param x,y Two \code{SpatialExperiment} objects with the same pixels, for
#'   example the elements of the list returned by one
#'   \code{SEraster::rasterizeGeneExpression()} call. Alternatively \code{x}
#'   is a list of two such objects and \code{y} is \code{NULL}, such as the
#'   list returned by \code{SEraster::rasterizeGeneExpression(list(a = ...,
#'   b = ...))}; the names of the list, if any, are used as the sample labels
#'   (otherwise \code{"x"} and \code{"y"}).
#' @param assay The assay to use: a name or an index, used for both objects,
#'   or two names or two indices, one per object. Default: the first assay.
#' @param genes The genes (row names) to compare. Default \code{NULL}: every
#'   gene present in both objects (a message says how many were dropped).
#'   Genes given here must be present in both objects.
#' @param tests Which tests to run: \code{"correlation"} (the spatial
#'   correlation test with permutations) and/or \code{"similarity"} (no
#'   permutations).
#' @param nPermutations The largest number of permutations per gene; the
#'   smallest possible p-value is \code{1 / (nPermutations + 1)}. Default
#'   10000.
#' @param exceedances The number of null correlations at least as extreme as
#'   the observed one after which a gene stops early (\eqn{h} of Besag and
#'   Clifford). The relative standard error of a p-value is about
#'   \code{1 / sqrt(exceedances)}. Default 10; \code{Inf} gives every gene
#'   exactly \code{nPermutations} permutations.
#' @param delta The candidate smoothing bandwidths, as shares of the pixels in
#'   the kernel's window (values in (0, 1]). Values that leave fewer than 2
#'   pixels in the window (\eqn{\lfloor N\delta \rfloor < 2}{floor(N * delta)
#'   < 2} for \eqn{N} shared pixels) are dropped with a message. The default
#'   adds 0.01 and 0.05 to \code{seq(0.1, 0.9, 0.1)}, the grid the authors
#'   used for the kidney and MERFISH analyses: with the shorter grid, delta
#'   star is often its smallest value.
#' @param maxDistPrctile The quantile of the pairwise pixel distances up to
#'   which the variograms are computed. Default 0.25.
#' @param minDetected The smallest share of the shared pixels at which a gene
#'   must be detected, in each sample, to be tested for correlation (see
#'   \emph{Rarely detected genes} in Details); \code{0} tests every gene that
#'   is not constant. Default \code{NULL}, which depends on \code{surrogate}:
#'   \code{0} (no filter) with \code{"remap"}, whose p-values are calibrated
#'   for rarely detected genes, and \eqn{\sqrt{N}}{sqrt(N)} of the \eqn{N}
#'   shared pixels (18 of 311, 47 of 2170) with \code{"gaussian"}, whose
#'   p-values are far too small for such genes. The number of pixels asked
#'   for is recorded in \code{attr(result, "params")$minDetectedPixels}.
#' @param surrogate How the surrogates of the correlation test take their
#'   values (see \emph{Rank-remapped surrogates} in Details): \code{"remap"}
#'   (the default) gives each surrogate the gene's own values, in the
#'   surrogate's rank order, so that it has exactly the gene's distribution
#'   of values, which calibrates the p-values of sparse and skewed genes;
#'   \code{"gaussian"} uses them as the smoothing and the added noise leave
#'   them, the surrogates of the legacy functions and of the published
#'   analyses.
#' @param foldChange The similarity band: pixels with
#'   \eqn{|\log_2(y / x)| \le}{|log2(y / x)| <=} \code{foldChange} are
#'   similar. Default 1 (within two-fold).
#' @param minQuantile The quantile of each gene's values in each sample used
#'   as its expression threshold for the similarity. Default 0.05.
#' @param minPixels The smallest share of the pixels that must pass the
#'   threshold for the similarity to be computed. Default 0.1.
#' @param adjustMethod The multiple-testing adjustment of \code{p} across
#'   genes, one of \code{stats::p.adjust.methods}. Default \code{"BH"}.
#' @param seed An integer seed for the permutations and the noise. Default 0.
#' @param nThreads The number of threads. Default
#'   \code{getOption("STcompare.nThreads", 1L)}, so it can be set once per
#'   session with \code{options(STcompare.nThreads = 8)}. The results do not
#'   depend on it.
#' @param progress Show a progress line (with an estimate of the remaining
#'   time) while the permutations run. Default \code{interactive()}. It is
#'   written as a message, so \code{suppressMessages()} also hides it.
#' @param verbose Show informative messages: genes or deltas that were
#'   dropped, genes whose delta star is at the edge of the grid, p-values
#'   limited by \code{nPermutations}, and the total time. Default \code{TRUE}.
#' @param keepNulls Keep the null correlations and delta stars of every
#'   permutation in \code{attr(result, "details")}. Default \code{FALSE}.
#'
#' @return A data frame of class \code{"STcompareResult"} with one row per gene
#'   (row names are the genes). Columns of a test that was not run are
#'   absent.
#' \describe{
#'   \item{\code{gene}}{The gene.}
#'   \item{\code{nPixels}}{The number of shared pixels.}
#'   \item{\code{r}}{Pearson's correlation of the gene in \code{x} and
#'   \code{y}.}
#'   \item{\code{pNaive}}{The p-value of \code{cor.test()}, which assumes
#'   independent pixels; for reference only: under spatial autocorrelation it
#'   is far too small.}
#'   \item{\code{p}}{The permutation p-value (see Details).}
#'   \item{\code{padj}}{\code{p} adjusted across the genes with
#'   \code{adjustMethod} (skipped and failed genes are not counted).}
#'   \item{\code{pX}, \code{pY}}{\eqn{(b + 1) / (L + 1)} for the direction
#'   that permutes \code{x} and for the one that permutes \code{y}, at the
#'   common number of permutations \eqn{L}. They are descriptive: when a gene
#'   stops early, the direction that did not reach \code{exceedances} has an
#'   inflated value. Use \code{p}.}
#'   \item{\code{nPermutations}}{The number of permutations used (\eqn{L}).}
#'   \item{\code{stop}}{Why the permutations stopped: \code{"exceedances"}
#'   (early), \code{"limit"} (\code{nPermutations} reached),
#'   \code{"skipped"} (not tested) or \code{"failed"}.}
#'   \item{\code{deltaStarMedianX}, \code{deltaStarMedianY}}{The median
#'   delta star of the permutations of \code{x} and of \code{y}.}
#'   \item{\code{deltaGridEdge}}{\code{TRUE} when, in either direction, more
#'   than half of the permutations chose the smallest or the largest delta of
#'   the grid: the best bandwidth may lie outside the grid. \code{NA} with
#'   fewer than 3 deltas.}
#'   \item{\code{similarity}}{The share of the kept pixels with
#'   \eqn{|\log_2(y / x)| \le}{|log2(y / x)| <=} \code{foldChange}.}
#'   \item{\code{dissimilarityX}}{The share with
#'   \eqn{\log_2(y / x) < -}{log2(y / x) < -}\code{foldChange} (higher in
#'   \code{x}).}
#'   \item{\code{dissimilarityY}}{The share with
#'   \eqn{\log_2(y / x) >}{log2(y / x) >} \code{foldChange} (higher in
#'   \code{y}).}
#'   \item{\code{nPixelsSimilarity}}{The number of kept pixels.}
#'   \item{\code{thresholdX}, \code{thresholdY}}{The gene's expression
#'   thresholds in \code{x} and \code{y}.}
#'   \item{\code{status}}{\code{"ok"}, \code{"skipped"} (not tested for
#'   correlation: constant or rarely detected in a sample) or
#'   \code{"failed"} (see Details).}
#'   \item{\code{message}}{Why a gene was skipped or failed, or why its
#'   similarity is \code{NA} (negative values); empty otherwise.}
#' }
#' The attributes hold \code{params} (the arguments after defaults, among
#' them the \code{surrogate} mode used and \code{minDetectedPixels}, the
#' number of detected pixels asked for; the sample labels, the assays, the
#' number of pixels and the delta grid actually used, \code{deltaGrid}),
#' \code{call}, \code{runtime} and, with
#' \code{keepNulls = TRUE}, \code{details}: one list per gene with the null
#' correlations (\code{nullX}, \code{nullY}) and delta stars
#' (\code{deltaStarX}, \code{deltaStarY}) of its permutations.
#' \code{print()}, \code{summary()} and \code{as.data.frame()} (which returns
#' a plain data frame) have methods for it; a selection of rows keeps the
#' class, and a selection of columns is a plain data frame.
#'
#' @references Viladomat J, Mazumder R, McInturff A, McCauley DJ, Hastie T
#'   (2014). Assessing the significance of global and local correlations under
#'   spatial autocorrelation: a nonparametric approach. \emph{Biometrics}
#'   70(2):409-418. \doi{10.1111/biom.12139}
#'
#'   Besag J, Clifford P (1991). Sequential Monte Carlo p-values.
#'   \emph{Biometrika} 78(2):301-304. \doi{10.1093/biomet/78.2.301}
#'
#'   Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD (1992). Testing
#'   for nonlinearity in time series: the method of surrogate data.
#'   \emph{Physica D} 58(1-4):77-94. \doi{10.1016/0167-2789(92)90102-S}
#'
#' @seealso \code{\link{spatialCorrelationGeneExp}()} and
#'   \code{\link{spatialSimilarity}()}, which reproduce the published
#'   analyses exactly. The tutorials:
#'   \code{vignette("getting-started-with-STcompare", package = "STcompare")},
#'   \code{vignette("how-STcompare-works", package = "STcompare")} and
#'   \code{vignette("parameters-performance-reproducibility", package = "STcompare")}.
#'
#' @export
#'
#' @examples
#' data(speKidney)
#' # rasterize the three samples together, onto one grid of pixels
#' rast <- SEraster::rasterizeGeneExpression(speKidney, assay_name = "counts",
#'                                           resolution = 0.2, fun = "mean",
#'                                           square = FALSE)
#'
#' # A and B have opposite spatial patterns, A and C similar ones
#' ab <- compareSpatial(list(A = rast$A, B = rast$B), nPermutations = 999,
#'                      nThreads = 2)
#' ab
#' ac <- compareSpatial(list(A = rast$A, C = rast$C), nPermutations = 999,
#'                      nThreads = 2)
#' summary(ac)
#' as.data.frame(ac)[, c("gene", "r", "p", "nPermutations", "similarity")]
compareSpatial <- function(x, y = NULL, assay = 1, genes = NULL,
                           tests = c("correlation", "similarity"),
                           nPermutations = 10000, exceedances = 10,
                           delta = c(0.01, 0.05, seq(0.1, 0.9, 0.1)), maxDistPrctile = 0.25,
                           minDetected = NULL, surrogate = c("remap", "gaussian"),
                           foldChange = 1, minQuantile = 0.05, minPixels = 0.1,
                           adjustMethod = "BH", seed = 0L,
                           nThreads = getOption("STcompare.nThreads", 1L),
                           progress = interactive(), verbose = TRUE, keepNulls = FALSE) {
  call <- match.call()
  t0 <- proc.time()
  .stc_local_rng()

  # ---- arguments (all checked before anything is computed) ----
  if (!is.character(tests) || !length(tests) || anyNA(tests) || !all(tests %in% c("correlation", "similarity"))) {
    stop("tests must be \"correlation\", \"similarity\" or both")
  }
  tests <- intersect(c("correlation", "similarity"), tests)
  do_cor <- "correlation" %in% tests
  do_sim <- "similarity" %in% tests
  surrogate <- match.arg(surrogate)
  nPermutations <- .stc_check_nperm(nPermutations)
  if (nPermutations > .Machine$integer.max - 1L) stop("nPermutations is too large")
  if (!is.numeric(exceedances) || length(exceedances) != 1L || is.na(exceedances) || exceedances < 1 ||
      (is.finite(exceedances) && exceedances != round(exceedances))) {
    stop("exceedances must be a positive whole number or Inf")
  }
  if (!is.numeric(delta) || !length(delta) || !all(is.finite(delta)) || any(delta <= 0 | delta > 1)) {
    stop("delta must be a numeric vector of values in (0, 1]")
  }
  .stc_check_unit_interval(maxDistPrctile, "maxDistPrctile", zero = FALSE)
  if (!is.null(minDetected)) .stc_check_unit_interval(minDetected, "minDetected")
  if (!is.numeric(foldChange) || length(foldChange) != 1L || is.na(foldChange) || foldChange < 0) {
    stop("foldChange must be a single non-negative number")
  }
  .stc_check_unit_interval(minQuantile, "minQuantile")
  .stc_check_unit_interval(minPixels, "minPixels")
  if (!is.character(adjustMethod) || length(adjustMethod) != 1L || !adjustMethod %in% stats::p.adjust.methods) {
    stop("adjustMethod must be one of: ", paste(stats::p.adjust.methods, collapse = ", "))
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != round(seed) ||
      abs(seed) > .Machine$integer.max) {
    stop("seed must be a single integer")
  }
  seed <- as.integer(seed)
  nThreads <- .stc_threads(nThreads)
  for (a in c("progress", "verbose", "keepNulls")) {
    v <- get(a)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) stop(sprintf("%s must be TRUE or FALSE", a))
  }

  # ---- input: pixels, genes and values ----
  d <- .stc_compare_input(x, y, assay, genes, verbose)
  N <- nrow(d$pos)
  G <- length(d$genes)
  lab <- d$labels
  badX <- colSums(!is.finite(d$X)) > 0L
  badY <- colSums(!is.finite(d$Y)) > 0L
  finite <- !badX & !badY
  status <- ifelse(finite, "ok", "failed")
  msg <- rep("", G)
  msg[!finite] <- paste0("missing or infinite values in ", .stc_which_samples(badX, badY, lab)[!finite])

  out <- data.frame(gene = d$genes, nPixels = rep(N, G), stringsAsFactors = FALSE)
  grid <- NULL
  details <- NULL
  min_px <- NA_integer_
  if (do_cor) {
    grid <- .stc_compare_grid(delta, N, verbose)
    # a tested gene is detected at this many pixels of each sample at least: minDetected, or by default none
    # with remapped surrogates (their p-values are calibrated for rarely detected genes) and sqrt(N) with
    # gaussian surrogates (theirs are far too small for such genes)
    min_px <- as.integer(if (!is.null(minDetected)) ceiling(minDetected * N - 1e-9)
                         else if (surrogate == "remap") 0L else ceiling(sqrt(N)))
    cr <- .stc_compare_correlation(d, grid, finite, min_px, nPermutations, exceedances, maxDistPrctile, surrogate,
                                   seed, nThreads, progress, keepNulls)
    grid <- cr$grid
    status <- ifelse(finite, cr$status, status)
    msg <- ifelse(finite, cr$message, msg)
    cr$table$padj <- stats::p.adjust(cr$table$p, method = adjustMethod)
    out <- cbind(out, cr$table[, c("r", "pNaive", "p", "padj", "pX", "pY", "nPermutations", "stop",
                                   "deltaStarMedianX", "deltaStarMedianY", "deltaGridEdge")])
    details <- cr$details
  }
  negative <- integer(0)
  if (do_sim) {
    sim <- .stc_similarity(d$X, d$Y, minQuantile = minQuantile, minPixels = minPixels, foldChange = foldChange)
    # a fold change needs values that are not negative: no similarity for genes with negative values
    negX <- finite & colSums(d$X < 0) > 0L
    negY <- finite & colSums(d$Y < 0) > 0L
    negative <- which(negX | negY)
    sim[!finite | negX | negY, ] <- NA
    note <- paste0("negative values in ", .stc_which_samples(negX, negY, lab)[negative], ": no similarity")
    msg[negative] <- ifelse(nzchar(msg[negative]), paste(msg[negative], note, sep = "; "), note)
    out <- cbind(out, sim)
  }
  out$status <- status
  out$message <- msg
  rownames(out) <- d$genes

  params <- list(samples = lab, assay = d$assay, nPixels = N, genes = d$genes, tests = tests,
                 nPermutations = nPermutations, exceedances = exceedances, delta = delta, deltaGrid = grid,
                 maxDistPrctile = maxDistPrctile, minDetected = minDetected, minDetectedPixels = min_px,
                 surrogate = surrogate, foldChange = foldChange, minQuantile = minQuantile, minPixels = minPixels,
                 adjustMethod = adjustMethod, seed = seed, nThreads = nThreads, progress = progress,
                 verbose = verbose, keepNulls = keepNulls)
  class(out) <- c("STcompareResult", "data.frame")
  attr(out, "params") <- params
  attr(out, "call") <- call
  attr(out, "runtime") <- proc.time() - t0
  if (keepNulls && do_cor) attr(out, "details") <- details

  # ---- messages and warnings ----
  if (verbose && do_cor) {
    skipped <- sum(out$status == "skipped")
    if (skipped) {
      nc <- sum(cr$constant)
      message(sprintf("compareSpatial: %d of %d genes are not tested for correlation (status \"skipped\"): %s", skipped, G,
                      paste(c(if (nc) sprintf("%d constant in a sample", nc),
                              if (skipped > nc) sprintf("%d detected in fewer than %d of the %d shared pixels of a sample (minDetected)",
                                                        skipped - nc, min_px, N)),
                            collapse = " and ")))
    }
    .stc_compare_diagnostics(out, grid, nPermutations, adjustMethod)
  }
  failed <- which(out$status == "failed")
  if (length(failed)) {
    warning(sprintf("compareSpatial: %d of %d genes failed (status \"failed\"; see the message column): %s", length(failed),
                    G, .stc_some(paste0(out$gene[failed], " (", out$message[failed], ")"))), call. = FALSE)
  }
  if (length(negative)) {
    warning(sprintf(paste0("compareSpatial: the similarity compares fold changes, which need values that are not ",
                           "negative (counts or normalized values on a linear scale); it is NA for %d gene%s with ",
                           "negative values: %s"), length(negative), if (length(negative) == 1L) "" else "s",
                    .stc_some(out$gene[negative])), call. = FALSE)
  }
  if (verbose) {
    perms <- if (do_cor) sum(out$nPermutations, na.rm = TRUE) else 0
    message(sprintf("compareSpatial: %d gene%s on %d shared pixels in %s%s", G, if (G == 1L) "" else "s", N,
                    .stc_elapsed(t0),
                    if (do_cor) sprintf(" (%s permutations in total, %d thread%s)", .stc_format_count(perms),
                                        nThreads, if (nThreads == 1L) "" else "s") else ""))
  }
  out
}

# ---------------------------------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------------------------------

.stc_check_unit_interval <- function(v, name, zero = TRUE) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v > 1 || v < 0 || (!zero && v == 0)) {
    stop(sprintf("%s must be a single number in %s, 1]", name, if (zero) "[0" else "(0"))
  }
  invisible(v)
}

# "x", "y" or "x and y" (the sample labels lab) for genes with a property in x (inX) and/or in y (inY).
.stc_which_samples <- function(inX, inY, lab) {
  ifelse(inX & inY, paste(lab, collapse = " and "), ifelse(inX, lab[1], lab[2]))
}

# The first k elements of v for a message, and how many more there are.
.stc_some <- function(v, k = 10L) {
  paste0(paste(utils::head(v, k), collapse = "; "), if (length(v) > k) sprintf("; and %d more", length(v) - k) else "")
}

# The two samples of compareSpatial(), checked: sample labels, resolved assays, genes, the shared pixels
# (in the column order of x, as in the legacy functions, so that runs with more than 1000 pixels draw the
# same variogram subsample) and their coordinates, and pixels x genes matrices of the values.
.stc_compare_input <- function(x, y, assay, genes, verbose) {
  if (is.null(y)) {
    if (!is.list(x) || inherits(x, "SpatialExperiment") || length(x) != 2L) {
      stop("give two SpatialExperiment objects: compareSpatial(x, y), or a list of two as x")
    }
    nm <- names(x)
    labels <- if (!is.null(nm) && all(nzchar(nm)) && !anyNA(nm) && !anyDuplicated(nm)) nm else c("x", "y")
    objs <- unname(x)
  } else {
    objs <- list(x, y)
    labels <- c("x", "y")
  }
  for (i in 1:2) {
    if (!inherits(objs[[i]], "SpatialExperiment")) {
      stop(sprintf("%s must be a SpatialExperiment object (for example rasterized with SEraster)", labels[i]))
    }
  }
  # assays
  if (length(assay) == 1L) assay <- rep(assay, 2L)
  if (length(assay) != 2L || anyNA(assay) || !(is.character(assay) || is.numeric(assay))) {
    stop("assay must be one name or index, or two (one per sample)")
  }
  assay <- as.list(assay)
  assay_label <- character(2)
  for (i in 1:2) {
    an <- SummarizedExperiment::assayNames(objs[[i]])
    na <- length(SummarizedExperiment::assays(objs[[i]]))
    a <- assay[[i]]
    ok <- if (is.character(a)) a %in% an else a == round(a) && a >= 1 && a <= na
    if (!ok) {
      stop(sprintf("%s has no assay %s; its assays are: %s", labels[i],
                   if (is.character(a)) dQuote(a, FALSE) else a,
                   if (length(an)) paste(an, collapse = ", ") else paste0(na, " unnamed")))
    }
    # the assay's name where it has one (recorded in the result), otherwise its index
    assay_label[i] <- if (is.character(a)) a else if (length(an) >= a && nzchar(an[a])) an[a] else as.character(a)
  }
  # genes
  rn <- lapply(objs, rownames)
  for (i in 1:2) {
    if (is.null(rn[[i]]) || anyNA(rn[[i]])) stop(sprintf("the genes (row names) of %s are missing", labels[i]))
    dup <- unique(rn[[i]][duplicated(rn[[i]])])
    if (length(dup)) {
      stop(sprintf("the gene names (row names) of %s must be unique; duplicated: %s", labels[i],
                   paste(utils::head(dup, 5L), collapse = ", ")))
    }
  }
  if (is.null(genes)) {
    genes <- intersect(rn[[1]], rn[[2]])
    if (!length(genes)) stop(sprintf("%s and %s have no gene in common", labels[1], labels[2]))
    only <- c(length(setdiff(rn[[1]], genes)), length(setdiff(rn[[2]], genes)))
    if (verbose && any(only > 0L)) {
      message(sprintf("compareSpatial: comparing the %d genes in both samples (%d genes only in %s and %d only in %s are dropped)",
                      length(genes), only[1], labels[1], only[2], labels[2]))
    }
  } else {
    if (is.factor(genes)) genes <- as.character(genes)
    if (!is.character(genes) || !length(genes) || anyNA(genes)) stop("genes must be a character vector of gene names")
    if (anyDuplicated(genes)) stop("genes must not contain duplicates")
    for (i in 1:2) {
      miss <- setdiff(genes, rn[[i]])
      if (length(miss)) {
        stop(sprintf("%d of the genes are not in %s, for example %s", length(miss), labels[i],
                     paste(utils::head(miss, 5L), collapse = ", ")))
      }
    }
  }
  # shared pixels, matched by name and checked by coordinates
  cn <- lapply(objs, colnames)
  for (i in 1:2) {
    if (is.null(cn[[i]]) || anyNA(cn[[i]])) stop(sprintf("the pixels (column names) of %s are missing", labels[i]))
    if (anyDuplicated(cn[[i]])) stop(sprintf("the pixel names (column names) of %s must be unique", labels[i]))
  }
  shared <- intersect(cn[[1]], cn[[2]])
  if (length(shared) < 3L) {
    stop(sprintf("%s and %s share %d pixel name(s); at least 3 are needed. Rasterize both samples in one call, ",
                 labels[1], labels[2], length(shared)),
         "SEraster::rasterizeGeneExpression(list(x, y), ...), so that they have the same pixels")
  }
  ix <- match(shared, cn[[1]])
  iy <- match(shared, cn[[2]])
  coords <- lapply(1:2, function(i) {
    m <- SpatialExperiment::spatialCoords(objs[[i]])
    if (is.null(m) || ncol(m) < 2L) stop(sprintf("%s needs two spatial coordinates", labels[i]))
    m <- m[if (i == 1L) ix else iy, 1:2, drop = FALSE]
    m <- matrix(as.double(m), ncol = 2L)
    if (!all(is.finite(m))) stop(sprintf("the coordinates of the shared pixels of %s must be finite", labels[i]))
    m
  })
  span <- max(apply(rbind(coords[[1]], coords[[2]]), 2L, function(v) diff(range(v))))
  off <- pmax(abs(coords[[1]][, 1] - coords[[2]][, 1]), abs(coords[[1]][, 2] - coords[[2]][, 2]))
  moved <- which(off > 1e-6 * span)
  if (length(moved)) {
    stop(sprintf(paste0("%d of the %d pixels shared by %s and %s have different coordinates in the two samples ",
                        "(for example %s, which is %.4g apart). Pixel names alone do not identify locations: ",
                        "rasterize both samples in one call, SEraster::rasterizeGeneExpression(list(x, y), ...), ",
                        "so that they are on one grid"),
                 length(moved), length(shared), labels[1], labels[2], shared[moved[1]], off[moved[1]]))
  }
  values <- function(i, cols) {
    m <- SummarizedExperiment::assay(objs[[i]], assay[[i]], withDimnames = TRUE)
    m <- as.matrix(m[genes, cols, drop = FALSE])
    if (!is.numeric(m) && !is.logical(m)) stop(sprintf("the assay of %s must be numeric", labels[i]))
    m <- t(m)
    storage.mode(m) <- "double"
    dimnames(m) <- list(shared, genes)
    m
  }
  list(labels = labels, assay = assay_label, genes = genes, pixels = shared, pos = coords[[1]],
       X = values(1L, ix), Y = values(2L, iy))
}

# The delta grid: sorted distinct values with floor(N * delta + 1e-12) >= 2 (the smoother's window must
# hold at least 2 pixels; the engine computes the window as (int)(N * delta + 1e-12)).
.stc_compare_grid <- function(delta, N, verbose) {
  grid <- sort(unique(as.double(delta)))
  keep <- floor(N * grid + 1e-12) >= 2
  if (!any(keep)) {
    stop(sprintf("no delta is usable with %d shared pixels: the smoothing window must hold at least 2 pixels, so delta must be at least 2 / %d = %.3g",
                 N, N, 2 / N))
  }
  if (verbose && !all(keep)) {
    message(sprintf("compareSpatial: delta %s dropped: with %d shared pixels the smoothing window must hold at least 2 pixels (delta >= %.3g)",
                    paste(format(grid[!keep]), collapse = ", "), N, 2 / N))
  }
  grid[keep]
}

# ---------------------------------------------------------------------------------------------------
# The correlation test
# ---------------------------------------------------------------------------------------------------

# The batch schedule of compareSpatial(): the first batch of permutations, the growth of later batches,
# the largest batch and the permutations per work item. The results do not depend on it (only the work and
# the granularity of the progress display do); the tests change it with option STcompare.compare_schedule.
.stc_compare_schedule <- function() {
  s <- utils::modifyList(list(first_batch = 16L, growth = 2, max_batch = 1024L, chunk = 16L),
                         as.list(getOption("STcompare.compare_schedule", list())))
  s
}

# The correlation columns of compareSpatial() for every gene, the status ("ok", "skipped", "failed") and
# message of every gene, the details (keepNulls), the delta grid used, and which genes are constant. Genes
# that are not finite (finite = FALSE; their status and message are set by the caller), constant in a
# sample, or detected in fewer than min_px pixels of a sample are not given to the engine. A gene is
# detected at a pixel where its value is above its lowest value (a constant gene is detected nowhere).
# surrogate: "gaussian" or "remap" (see .stc_engine_correlate()).
.stc_compare_correlation <- function(d, grid, finite, min_px, n_max, h, maxDistPrctile, surrogate, seed, nThreads,
                                     progress, keepNulls) {
  G <- length(d$genes)
  N <- nrow(d$X)
  lab <- d$labels
  naive <- .stc_cor_tests(d$X, d$Y)
  detected <- function(M) {
    k <- as.integer(colSums(M > rep(apply(M, 2L, min), each = N)))
    k[!finite] <- NA_integer_
    k
  }
  kX <- detected(d$X)
  kY <- detected(d$Y)
  constX <- finite & kX == 0L
  constY <- finite & kY == 0L
  constant <- constX | constY
  rare <- finite & !constant & (kX < min_px | kY < min_px)
  status <- ifelse(constant | rare, "skipped", "ok")
  msg <- rep("", G)
  msg[constant] <- paste0("constant in ", .stc_which_samples(constX, constY, lab)[constant],
                          " (zero variance): not tested")
  msg[rare] <- sprintf("detected in %d (%s) and %d (%s) of the %d shared pixels; minDetected asks for %d in each: not tested",
                       kX[rare], lab[1], kY[rare], lab[2], N, min_px)
  testable <- which(finite & !constant & !rare)
  tab <- data.frame(r = ifelse(finite, naive$r, NA_real_), pNaive = ifelse(finite, naive$p, NA_real_),
                    p = NA_real_, padj = NA_real_, pX = NA_real_, pY = NA_real_,
                    nPermutations = NA_integer_, stop = ifelse(constant | rare, "skipped", "failed"),
                    deltaStarMedianX = NA_real_, deltaStarMedianY = NA_real_, deltaGridEdge = NA,
                    stringsAsFactors = FALSE)
  details <- NULL
  if (keepNulls) {
    empty <- list(nullX = numeric(0), nullY = numeric(0), deltaStarX = numeric(0), deltaStarY = numeric(0))
    details <- stats::setNames(rep(list(empty), G), d$genes)
  }
  result <- function() {
    list(table = tab, status = status, message = msg, details = details, grid = grid, constant = constant)
  }
  if (!length(testable)) return(result())

  sch <- .stc_compare_schedule()
  prepare <- function(grid) {
    .stc_engine_prepare(d$X[, testable, drop = FALSE], d$Y[, testable, drop = FALSE], d$pos, grid, grid, seed,
                        maxDistPrctile, nThreads, FALSE, "pair", sch$chunk, "cpp", streams = "independent",
                        keepNulls = keepNulls, surrogate = surrogate)
  }
  state <- prepare(grid)
  plan <- state$plan
  if (!isTRUE(plan$ok) || plan$nbins < 2L) {
    stop(sprintf("the variograms cannot be computed on the shared pixels with maxDistPrctile = %g: %s",
                 maxDistPrctile,
                 if (isTRUE(plan$ok)) sprintf("%d distance bin(s) with at least 2 pixel pairs; 2 are needed", plan$nbins)
                 else plan$reason),
         ". Try a larger maxDistPrctile.")
  }
  # a smoother can fail for a delta whatever the gene (vertex space exhausted: duplicated coordinates)
  dl <- .stc_engine_deltas(state$session)
  broken <- dl$status != 0L
  if (any(broken)) {
    warning(sprintf("compareSpatial: delta %s dropped: %s", paste(format(dl$delta[broken]), collapse = ", "),
                    paste(unique(dl$message[broken]), collapse = "; ")), call. = FALSE)
    grid <- setdiff(grid, dl$delta[broken])
    if (!length(grid)) stop("no delta is usable: the smoother cannot be built for any delta of the grid")
    state <- prepare(grid)
  }

  K <- length(grid)
  Gt <- length(testable)
  adaptive <- if (is.finite(h)) list(h = h, n_max = n_max, first_batch = sch$first_batch, growth = sch$growth)
  p <- if (progress) .stc_progress(Gt, n_max, h, sch) else NULL
  if (!is.null(p)) on.exit(p$close(), add = TRUE)
  .stc_engine_extend(state, seq_len(Gt), nPermutations = n_max, adaptive = adaptive, batch = sch$max_batch,
                     progress = p)

  # collect
  u <- .stc_engine_units(state$session)
  tr <- .stc_engine_task_results(state$session, seq_len(2L * Gt) - 1L, FALSE)
  st <- u$state
  ok <- st != 3L
  counts <- matrix(unlist(u$counts, use.names = FALSE), nrow = 2L)
  L <- ifelse(ok, u$len, NA_integer_)
  bX <- ifelse(ok, counts[1L, ], NA_integer_)
  bY <- ifelse(ok, counts[2L, ], NA_integer_)
  # delta* per task from the counts per grid position (the grid is sorted): the median as R's median()
  # computes it from the L values, and the share at the two ends of the grid
  dc <- matrix(unlist(lapply(tr, `[[`, "dstar_count"), use.names = FALSE), nrow = K)
  med <- .stc_median_from_counts(dc, grid)
  edge_share <- if (K >= 3L) (dc[1L, ] + dc[K, ]) / colSums(dc) else rep(NA_real_, 2L * Gt)
  edge <- matrix(edge_share > 0.5, nrow = 2L)
  t_idx <- testable
  tab$p[t_idx] <- .stc_p_combined(pmax(bX, bY), L, st == 1L)
  tab$pX[t_idx] <- .stc_p_direction(bX, L)
  tab$pY[t_idx] <- .stc_p_direction(bY, L)
  tab$nPermutations[t_idx] <- L
  tab$stop[t_idx] <- c("active", "exceedances", "limit", "failed")[st + 1L]
  tab$deltaStarMedianX[t_idx] <- ifelse(ok, med[2L * seq_len(Gt) - 1L], NA_real_)
  tab$deltaStarMedianY[t_idx] <- ifelse(ok, med[2L * seq_len(Gt)], NA_real_)
  tab$deltaGridEdge[t_idx] <- ifelse(ok, edge[1L, ] | edge[2L, ], NA)
  fail <- which(!ok)
  if (length(fail)) {
    dir <- ifelse(!is.na(u$fail_task[fail]) & u$fail_task[fail] %% 2L == 1L, lab[2], lab[1])
    msg[t_idx[fail]] <- paste0("permuting ", dir, ": ", u$message[fail])
    status[t_idx[fail]] <- "failed"
  }
  if (keepNulls) {
    for (q in which(ok)) {
      tX <- tr[[2L * q - 1L]]
      tY <- tr[[2L * q]]
      details[[t_idx[q]]] <- list(nullX = tX$nulls[, 1L], nullY = tY$nulls[, 1L],
                                  deltaStarX = grid[tX$dstar], deltaStarY = grid[tY$dstar])
    }
  }
  result()
}

# Medians of delta* per task (columns of counts: permutations per grid position, the grid sorted), equal
# to stats::median() of the delta* values themselves: the middle value, or mean() of the two middle ones
# (looked up in a table of mean(c(grid[i], grid[j])), so that the arithmetic is R's own).
.stc_median_from_counts <- function(counts, grid) {
  K <- length(grid)
  n <- colSums(counts)
  cum <- apply(counts, 2L, cumsum)
  if (K == 1L) cum <- matrix(cum, nrow = 1L)
  half <- (n + 1L) %/% 2L
  lo <- colSums(cum < rep(half, each = K)) + 1L
  hi <- colSums(cum < rep(half + 1L, each = K)) + 1L
  pairs <- outer(seq_len(K), seq_len(K), Vectorize(function(i, j) mean(c(grid[i], grid[j]))))
  med <- ifelse(n %% 2L == 1L, grid[pmin(lo, K)], pairs[cbind(pmin(lo, K), pmin(hi, K))])
  med[n == 0L] <- NA_real_
  med
}

# Messages about the results (verbose): delta* at the edge of the grid, and p-values limited by
# nPermutations.
.stc_compare_diagnostics <- function(out, grid, n_max, adjustMethod) {
  edge <- which(out$deltaGridEdge %in% TRUE)
  if (length(edge)) {
    message(sprintf(paste0("compareSpatial: for %d of %d genes (deltaGridEdge), more than half of the permutations ",
                           "of x or of y chose the smallest (%g) or the largest (%g) delta: the best smoothing may ",
                           "lie outside the grid; consider extending delta"),
                    length(edge), sum(out$status == "ok"), min(grid), max(grid)))
  }
  floor_p <- which(out$stop %in% "limit" & out$p == 1 / (n_max + 1))
  limited <- floor_p[!is.na(out$padj[floor_p]) & out$padj[floor_p] >= 0.05]
  if (length(limited)) {
    message(sprintf(paste0("compareSpatial: %d gene%s reached nPermutations = %d without any exceedance, so p is ",
                           "the smallest possible, 1 / (nPermutations + 1) = %.3g; with %d tested genes their ",
                           "adjusted p-values (%s) are not below 0.05. A larger nPermutations would resolve them"),
                    length(limited), if (length(limited) == 1L) "" else "s", n_max, 1 / (n_max + 1),
                    sum(!is.na(out$p)), adjustMethod))
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------------------------------
# Progress display
# ---------------------------------------------------------------------------------------------------

# A progress reporter for .stc_engine_extend() (dev/compare-spatial-spec.md, section 6). One line, written
# with message() (so suppressMessages() hides it) and rewritten in place ("\r") at most twice a second:
#   compareSpatial:  61% | genes done 640/1046 | 1.9M permutations | elapsed 0:12 | ETA 0:08
# The elapsed time is left out, and then the line cut, where it would be wider than the console
# (getOption("width")). $batch() is called by R before every batch, $tick(done) by the engine's main thread
# while a batch runs (done: the task-permutations of the batch finished so far; R is never called from a
# worker thread), and $finish() at the end, which completes the line with a newline ($close() does that after
# an error or an interrupt). The permutations shown are those the genes kept before the batch (the sum of the
# nPermutations column at the end) plus those computed in the batch so far. In adaptive mode the remaining
# work is a prediction: a gene with b exceedances (the larger direction) after L permutations is expected to
# stop near L * h / b (at nPermutations if b = 0), at the end of the batch that contains it; the percentage is
# kept from going backwards.
.stc_progress <- function(G, n_max, h, schedule, label = "compareSpatial") {
  p <- new.env(parent = emptyenv())
  # the ends of the batches of the schedule (as .stc_engine_extend() forms them)
  ends <- if (is.finite(h)) {
    sizes <- numeric(0)
    total <- 0
    step <- 0
    while (total < n_max) {
      s <- min(schedule$max_batch, max(1, floor(schedule$first_batch * schedule$growth^step)))
      total <- min(n_max, total + s)
      sizes <- c(sizes, total)
      step <- step + 1
    }
    sizes
  } else {
    unique(c(seq_len(n_max %/% schedule$max_batch) * schedule$max_batch, n_max))
  }
  p$t0 <- proc.time()[["elapsed"]]
  p$last <- -Inf
  p$width <- 0L
  p$done <- 0        # task-permutations of the finished batches
  p$batch_n <- 0     # task-permutations of the current batch
  p$remaining <- 0   # predicted task-permutations left at the start of the current batch
  p$kept <- 0        # permutations kept by the genes before the current batch
  p$genes_done <- 0L
  p$pct <- 0L
  p$open <- FALSE
  show <- function(d, final = FALSE) {
    now <- proc.time()[["elapsed"]]
    if (!final && now - p$last < 0.5) return(invisible(NULL))
    p$last <- now
    elapsed <- now - p$t0
    done <- p$done + d
    frac <- if (final) 1 else done / max(1, p$done + p$remaining)
    p$pct <- if (final) 100L else max(p$pct, min(99L, as.integer(floor(100 * frac))))
    eta <- if (final) 0 else if (frac > 0 && elapsed >= 1) elapsed * (1 - frac) / frac else NA
    middle <- if (is.finite(h)) {
      sprintf("genes done %d/%d | %s permutations", p$genes_done, G, .stc_format_count(p$kept + d / 2))
    } else {
      sprintf("%d gene%s x %d permutations", G, if (G == 1L) "" else "s", n_max)
    }
    eta <- if (is.na(eta)) "--:--" else .stc_format_clock(eta)
    width <- max(20L, getOption("width", 80L) - 1L)
    line <- sprintf("%s: %3d%% | %s | elapsed %s | ETA %s", label, p$pct, middle, .stc_format_clock(elapsed), eta)
    if (nchar(line) > width) line <- sprintf("%s: %3d%% | %s | ETA %s", label, p$pct, middle, eta)
    line <- substr(line, 1L, width)
    pad <- strrep(" ", max(0L, p$width - nchar(line)))
    p$width <- nchar(line)
    p$open <- !final
    message("\r", line, pad, appendLF = final)
    invisible(NULL)
  }
  p$batch <- function(u, group, b_from, b_to) {
    p$done <- p$done + p$batch_n
    p$batch_n <- 2 * length(group) * (b_to - b_from + 1)
    p$kept <- sum(u$len[u$state != 3L])  # (failed genes have no permutations in the result)
    active <- u$state == 0L
    p$genes_done <- sum(!active)
    len <- u$len[active]
    if (is.finite(h)) {
      b <- vapply(u$counts[active], max, 0L)
      Lhat <- ifelse(b == 0L, n_max, pmin(n_max, pmax(len + 1, ceiling(len * h / pmax(b, 1L)))))
    } else {
      Lhat <- rep(n_max, length(len))
    }
    end <- ends[findInterval(Lhat - 1, ends) + 1L]
    p$remaining <- 2 * sum(end - len)
    show(0)
  }
  p$tick <- function(d) {
    tryCatch(show(d), error = function(e) NULL)  # a display problem must never stop the computation
  }
  p$finish <- function(u) {
    p$done <- p$done + p$batch_n
    p$batch_n <- 0
    p$kept <- sum(u$len[u$state != 3L])
    p$genes_done <- G
    show(0, final = TRUE)
  }
  p$close <- function() {
    if (p$open) {
      p$open <- FALSE
      message("")
    }
  }
  p
}

.stc_format_clock <- function(s) {
  s <- max(0, round(s))
  if (s >= 3600) sprintf("%d:%02d:%02d", s %/% 3600, (s %% 3600) %/% 60, s %% 60) else sprintf("%d:%02d", s %/% 60, s %% 60)
}

.stc_format_count <- function(n) {
  if (n < 1e4) sprintf("%d", as.integer(round(n))) else if (n < 1e6) sprintf("%.1fk", n / 1e3) else sprintf("%.1fM", n / 1e6)
}

# ---------------------------------------------------------------------------------------------------
# Similarity
# ---------------------------------------------------------------------------------------------------

# The spatial similarity of the genes (columns) of X and Y (pixels in rows, the same pixels in both), with
# the semantics of spatialSimilarity() (R/packageFunction.R, which calls this helper), computed for all genes
# at once: per gene, the thresholds are t1 and t2 or else the minQuantile quantiles (quantile() type 7) of
# its values in X and in Y; pixels are kept when neither value is NA and the value in X is above its
# threshold or the value in Y is above its threshold; zeros are replaced by 1e-4; with l = log2(y / x) on the
# kept pixels, similarity is the share with -foldChange <= l <= foldChange, dissimilarityX the share with
# l < -foldChange and dissimilarityY the share with l > foldChange. A NaN l (values of opposite signs) counts
# in all three, as the original spatialSimilarity() counted it (spatialSimilarity() now stops on negative
# values; compareSpatial() keeps this rule). The three are NA when fewer than minPixels * nrow(X) pixels are
# kept. Returns a data.frame with similarity, dissimilarityX, dissimilarityY, nPixelsSimilarity (the number
# of kept pixels), thresholdX and thresholdY; a column with NA has an NA threshold. With details = TRUE, its
# attribute "details" is a list of kept (per gene, the row indices of the kept pixels), log (per gene, l of
# the kept pixels) and scored (TRUE for the genes with at least minPixels * nrow(X) kept pixels).
.stc_similarity <- function(X, Y, minQuantile = 0.05, minPixels = 0.1, foldChange = 1, t1 = NULL, t2 = NULL,
                            block = 1000L, details = FALSE) {
  N <- nrow(X)
  G <- ncol(X)
  thrX <- if (is.null(t1)) .stc_col_quantile(X, minQuantile) else rep_len(as.double(t1), G)
  thrY <- if (is.null(t2)) .stc_col_quantile(Y, minQuantile) else rep_len(as.double(t2), G)
  sim <- dissX <- dissY <- rep(NA_real_, G)
  kept <- integer(G)
  if (details) {
    kept_rows <- logs <- vector("list", G)
    scored <- logical(G)
  }
  for (cols in split(seq_len(G), (seq_len(G) - 1L) %/% block)) {  # blocks of genes bound the memory
    x <- X[, cols, drop = FALSE]
    y <- Y[, cols, drop = FALSE]
    keep <- !is.na(x) & !is.na(y) & (x > rep(thrX[cols], each = N) | y > rep(thrY[cols], each = N))
    keep[is.na(keep)] <- FALSE
    x[which(x == 0)] <- 1e-4
    y[which(y == 0)] <- 1e-4
    l <- suppressWarnings(log2(y / x))  # NaN for values of opposite signs
    nan <- is.na(l)
    n_keep <- colSums(keep)
    n_sim <- colSums(keep & (nan | (l >= -foldChange & l <= foldChange)))
    n_dx <- colSums(keep & (nan | l < -foldChange))
    n_dy <- colSums(keep & (nan | l > foldChange))
    ok <- !(n_keep < minPixels * N)
    kept[cols] <- as.integer(n_keep)
    sim[cols[ok]] <- n_sim[ok] / n_keep[ok]
    dissX[cols[ok]] <- n_dx[ok] / n_keep[ok]
    dissY[cols[ok]] <- n_dy[ok] / n_keep[ok]
    if (details) {
      scored[cols] <- ok
      for (k in seq_along(cols)) {
        rows <- which(keep[, k], useNames = FALSE)
        kept_rows[[cols[k]]] <- rows
        logs[[cols[k]]] <- as.vector(l[rows, k])
      }
    }
  }
  out <- data.frame(similarity = sim, dissimilarityX = dissX, dissimilarityY = dissY, nPixelsSimilarity = kept,
                    thresholdX = thrX, thresholdY = thrY)
  if (details) attr(out, "details") <- list(kept = kept_rows, log = logs, scored = scored)
  out
}

# quantile(X[, j], prob, names = FALSE) (type 7) of every column, with quantile()'s own arithmetic; NA for
# a column with NA or NaN (where quantile() stops).
.stc_col_quantile <- function(X, prob) {
  N <- nrow(X)
  if (N == 0L) return(rep(NA_real_, ncol(X)))
  prob <- max(0, min(1, prob))
  index <- 1 + max(N - 1, 0) * prob
  lo <- floor(index)
  hi <- ceiling(index)
  h <- index - lo
  vapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    if (anyNA(x)) return(NA_real_)
    x <- sort.int(x, partial = unique(c(lo, hi)))
    qs <- x[lo]
    if (index > lo && x[hi] != qs) qs <- (1 - h) * qs + h * x[hi]
    qs
  }, 0)
}

# ---------------------------------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------------------------------

#' Methods for compareSpatial() results
#'
#' \code{print()} shows a header (genes, pixels, samples and the test
#' settings), the number of genes with an adjusted p-value below 0.05 by the
#' sign of the correlation, how many genes stopped early or were skipped, and
#' the first rows. \code{summary()} returns these counts (and more) as an
#' object with its own \code{print()} method. \code{as.data.frame()} returns
#' the table as a plain data frame, without the class and the attributes.
#' Selecting rows with \code{[} keeps the class and the attributes; selecting
#' columns gives a plain data frame.
#'
#' @param x,object A result of \code{\link{compareSpatial}()}.
#' @param n The number of rows to print.
#' @param alpha The significance level of the counts of \code{summary()}.
#' @param i,j,drop Rows and columns to select, as for data frames.
#' @param row.names,optional,... Passed on or ignored, as for the generics.
#'
#' @return \code{print()} returns \code{x} invisibly; \code{summary()} an
#'   object of class \code{"summary.STcompareResult"}; \code{as.data.frame()} a
#'   data frame.
#'
#' @name STcompareResult
#' @examples
#' data(simRanPatternRasts)
#' # two independent simulated fields, rasterized onto one grid of pixels
#' res <- compareSpatial(simRanPatternRasts[[1]], simRanPatternRasts[[2]],
#'                       nPermutations = 199, nThreads = 2)
#' print(res)
#' summary(res)
#' head(as.data.frame(res))
NULL

# TRUE when df has the columns that print() and summary() count (a result whose columns were removed is
# printed and summarised as a data frame).
.stc_is_result <- function(df) {
  nm <- names(df)
  all(c("gene", "status") %in% nm) && (!"p" %in% nm || all(c("r", "padj", "nPermutations", "stop") %in% nm))
}

#' @rdname STcompareResult
#' @export
print.STcompareResult <- function(x, n = 6L, ...) {
  df <- as.data.frame(x)
  if (!.stc_is_result(df)) {
    print(df, ...)
    return(invisible(x))
  }
  pr <- attr(x, "params")
  G <- nrow(df)
  cat(sprintf("STcompareResult: %d gene%s", G, if (G == 1L) "" else "s"))
  if (!is.null(pr)) {
    cat(sprintf(" on %d shared pixels; %s vs %s (assay %s)", pr$nPixels, pr$samples[1], pr$samples[2],
                paste(unique(pr$assay), collapse = " / ")))
  }
  cat("\n")
  if ("p" %in% names(df)) {
    if (!is.null(pr)) {
      g <- pr$deltaGrid
      cat(if (is.finite(pr$exceedances)) {
        sprintf("Correlation: adaptive permutation p-values (stop at %s exceedances, at most %d permutations)",
                format(pr$exceedances), pr$nPermutations)
      } else {
        sprintf("Correlation: %d permutations per gene", pr$nPermutations)
      }, sprintf("; %d delta%s from %g to %g", length(g), if (length(g) == 1L) "" else "s", min(g), max(g)),
      sprintf("; %s surrogates", if (identical(pr$surrogate, "remap")) "rank-remapped" else "gaussian"), "\n", sep = "")
    } else {
      cat("Correlation:\n")
    }
    sig <- !is.na(df$padj) & df$padj < 0.05
    cat(sprintf("  padj < 0.05: %d gene%s (%d with r > 0, %d with r < 0)\n", sum(sig), if (sum(sig) == 1L) "" else "s",
                sum(sig & df$r > 0, na.rm = TRUE), sum(sig & df$r < 0, na.rm = TRUE)))
    st <- table(factor(df$stop, c("exceedances", "limit", "skipped", "failed")))
    cat(sprintf("  stopped early: %d; reached nPermutations: %d; skipped: %d; failed: %d\n", st[[1]], st[[2]],
                st[[3]], st[[4]]))
  }
  if ("similarity" %in% names(df)) {
    fc <- if (!is.null(pr)) pr$foldChange else NA
    cat(sprintf("Similarity%s: median %s (%d gene%s with NA)\n",
                if (!is.na(fc)) sprintf(" (|log2(y / x)| <= %g)", fc) else "",
                format(stats::median(df$similarity, na.rm = TRUE), digits = 3),
                sum(is.na(df$similarity)), if (sum(is.na(df$similarity)) == 1L) "" else "s"))
  }
  nf <- sum(df$status == "failed", na.rm = TRUE)
  if (nf) cat(sprintf("Failed: %d gene%s (see the message column)\n", nf, if (nf == 1L) "" else "s"))
  if (G) {
    cat("\n")
    print(utils::head(df, n), digits = 3, row.names = FALSE)
    if (G > n) cat(sprintf("... and %d more row%s\n", G - n, if (G - n == 1L) "" else "s"))
  }
  invisible(x)
}

#' @rdname STcompareResult
#' @export
summary.STcompareResult <- function(object, alpha = 0.05, ...) {
  df <- as.data.frame(object)
  if (!.stc_is_result(df)) return(summary(df, ...))
  pr <- attr(object, "params")
  s <- list(genes = nrow(df), nPixels = if (!is.null(pr)) pr$nPixels else NA_integer_,
            samples = if (!is.null(pr)) pr$samples else NULL,
            status = table(factor(df$status, c("ok", "skipped", "failed"))), alpha = alpha,
            runtime = attr(object, "runtime"))
  if ("p" %in% names(df)) {
    sig <- !is.na(df$padj) & df$padj < alpha
    s$significant <- c(total = sum(sig), positive = sum(sig & df$r > 0, na.rm = TRUE),
                       negative = sum(sig & df$r < 0, na.rm = TRUE))
    s$stop <- table(factor(df$stop, c("exceedances", "limit", "skipped", "failed")))
    L <- as.double(df$nPermutations[!is.na(df$nPermutations)])
    s$permutations <- c(total = sum(L), median = if (length(L)) stats::median(L) else NA_real_,
                        max = if (length(L)) max(L) else NA_real_)
    s$deltaGridEdge <- sum(df$deltaGridEdge %in% TRUE)
    s$settings <- if (!is.null(pr)) pr[c("nPermutations", "exceedances", "deltaGrid", "surrogate", "adjustMethod", "seed")]
  }
  if ("similarity" %in% names(df)) s$similarity <- summary(df$similarity)
  class(s) <- "summary.STcompareResult"
  s
}

#' @rdname STcompareResult
#' @export
print.summary.STcompareResult <- function(x, ...) {
  cat(sprintf("compareSpatial: %d gene%s", x$genes, if (x$genes == 1L) "" else "s"))
  if (!is.na(x$nPixels)) cat(sprintf(" on %d shared pixels", x$nPixels))
  if (!is.null(x$samples)) cat(sprintf(" (%s vs %s)", x$samples[1], x$samples[2]))
  cat(sprintf("; %d ok, %d skipped, %d failed\n", x$status[["ok"]], x$status[["skipped"]], x$status[["failed"]]))
  if (!is.null(x$significant)) {
    cat(sprintf("Correlation, padj < %g: %d (%d with r > 0, %d with r < 0)\n", x$alpha, x$significant[["total"]],
                x$significant[["positive"]], x$significant[["negative"]]))
    cat(sprintf("  stopped early: %d; reached nPermutations: %d; skipped: %d; failed: %d\n", x$stop[["exceedances"]],
                x$stop[["limit"]], x$stop[["skipped"]], x$stop[["failed"]]))
    cat(sprintf("  permutations per gene: median %s, max %s; %s in total\n", format(x$permutations[["median"]]),
                format(x$permutations[["max"]]), format(x$permutations[["total"]], big.mark = ",")))
    if (!is.null(x$settings)) {
      g <- x$settings$deltaGrid
      cat(sprintf("  settings: nPermutations = %d, exceedances = %s, %d deltas (%g to %g), %s surrogates, %s adjustment, seed %d\n",
                  x$settings$nPermutations, format(x$settings$exceedances), length(g), min(g), max(g),
                  if (identical(x$settings$surrogate, "remap")) "rank-remapped" else "gaussian",
                  x$settings$adjustMethod, x$settings$seed))
    }
    if (x$deltaGridEdge) cat(sprintf("  delta star at the edge of the grid: %d gene%s\n", x$deltaGridEdge,
                                     if (x$deltaGridEdge == 1L) "" else "s"))
  }
  if (!is.null(x$similarity)) {
    cat("Similarity:\n")
    print(x$similarity, digits = 3)
  }
  if (!is.null(x$runtime)) cat(sprintf("Run time: %.1f s\n", x$runtime[["elapsed"]]))
  invisible(x)
}

#' @rdname STcompareResult
#' @export
`[.STcompareResult` <- function(x, i, j, ..., drop) {
  out <- NextMethod()
  if (!is.data.frame(out)) return(out)
  if (!identical(names(out), names(x))) return(as.data.frame.STcompareResult(out))  # columns selected
  for (a in c("params", "call", "runtime")) attr(out, a) <- attr(x, a)                # rows selected
  out
}

#' @rdname STcompareResult
#' @export
as.data.frame.STcompareResult <- function(x, row.names = NULL, optional = FALSE, ...) {
  for (a in c("params", "call", "runtime", "details")) attr(x, a) <- NULL
  class(x) <- "data.frame"
  if (!is.null(row.names)) rownames(x) <- row.names
  x
}
