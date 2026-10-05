# data-raw/download_data.R
#
# Helper (sourced by the other data-raw scripts) that downloads the public inputs of the STcompare
# test data into a cache directory OUTSIDE the repository and verifies their md5 checksums.
#
#   source("data-raw/download_data.R")     # from the repository root
#   stc_download()                         # every file (130 MB)
#   stc_download(group = "aki")            # only the AKI kidney inputs
#   stc_download_dir()                     # where the downloads live
#
# Running it as a script downloads everything:  Rscript data-raw/download_data.R
#
# Cache layout (the root is $STCOMPARE_DATA_CACHE if set, otherwise tools::R_user_dir("STcompare", "cache")):
#   <root>/data-raw/            public downloads (this file)
#   <root>/data-raw/inputs/     rasterized inputs written by data-raw/build_inputs_*.R
#   <root>/data-raw/extracted/  unpacked tarballs
#
# A file already present with the expected md5 is not downloaded again. A cached file with a different
# md5, or a download whose md5 does not match, stops with an error that names the file, the URL and both
# checksums. Downloads are written to "<file>.part" and renamed only after the checksum has been verified.
# The md5 values below were computed from the copies used to build the committed fixtures; for the Zenodo
# files they equal the checksums published by the Zenodo API.
#
# Licences: every file is CC BY 4.0. The Zenodo licences were read from the Zenodo API. The 10x Genomics
# dataset page (https://www.10xgenomics.com/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0) states "This
# dataset is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license"; this was
# confirmed on 2026-10-03 from the page text as indexed by a web search engine, because the page blocks
# automated downloads. CC BY 4.0 requires attribution: see tests/testthat/fixtures/README.md.

stc_manifest <- local({
  zenodo <- function(record, file) sprintf("https://zenodo.org/records/%s/files/%s?download=1", record, file)
  tenx <- function(file) paste0("https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Mouse_Brain/", file)
  cc_by <- "CC-BY-4.0"
  aki_creators <- "Kalen Clifton, Jean Fan, Hamid Rabb"
  stalign_creators <- "Kalen Clifton, Manjari Anant, Gohta Aihara, Jean Fan"
  tenx_creators <- "10x Genomics"
  rows <- list(
    # group, file, bytes, md5, url, record, doi, licence, creators
    c("aki", "IL3_filtered_feature_bc_matrix.h5", 15103452, "a3ea3a2cb4c3e01403d0116ab580c157",
      zenodo(19074288, "IL3_filtered_feature_bc_matrix.h5"), "Zenodo 19074288", "10.5281/zenodo.19074288", cc_by, aki_creators),
    c("aki", "NL3_filtered_feature_bc_matrix.h5", 13822500, "6793a0f7ceeb082bce383748cfc2805a",
      zenodo(19074288, "NL3_filtered_feature_bc_matrix.h5"), "Zenodo 19074288", "10.5281/zenodo.19074288", cc_by, aki_creators),
    c("aki", "IL3_tissue_positions.csv", 80786, "4e3deabfaa2cd44b3b0231776b3aaf1b",
      zenodo(19074288, "IL3_tissue_positions.csv"), "Zenodo 19074288", "10.5281/zenodo.19074288", cc_by, aki_creators),
    c("aki", "NL3_tissue_positions.csv", 87431, "24b23088b7a22a53086f9b62b21fed29",
      zenodo(19074288, "NL3_tissue_positions.csv"), "Zenodo 19074288", "10.5281/zenodo.19074288", cc_by, aki_creators),
    c("aki", "aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz", 90127, "adede04898833d385eb52fa2032f2434",
      zenodo(19486091, "aki_region_onehot_STalign_to_ctrl_region_onehot_affine_only.csv.gz"),
      "Zenodo 19486091", "10.5281/zenodo.19486091", cc_by, aki_creators),
    c("brain", "STalign_S2R3_to_Visium.csv.gz", 16177162, "c092af84ca82a50539f4fdc0f51e8fdc",
      zenodo(10724029, "STalign_S2R3_to_Visium.csv.gz"), "Zenodo 10724029", "10.5281/zenodo.10724029", cc_by, stalign_creators),
    c("brain", "Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz", 45925294, "017fc078b4aa21e1df39027d0bf3ccad",
      tenx("Visium_FFPE_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"),
      "10x Genomics dataset Adult Mouse Brain (FFPE), Space Ranger 1.3.0 (adult-mouse-brain-ffpe-1-standard-1-3-0)",
      NA, cc_by, tenx_creators),
    c("brain", "Visium_FFPE_Mouse_Brain_spatial.tar.gz", 10632472, "a78bb08009c4defd791354dd5c292bfa",
      tenx("Visium_FFPE_Mouse_Brain_spatial.tar.gz"),
      "10x Genomics dataset Adult Mouse Brain (FFPE), Space Ranger 1.3.0 (adult-mouse-brain-ffpe-1-standard-1-3-0)",
      NA, cc_by, tenx_creators),
    c("celltype", "STalign_cell_type_transcriptional_correlations.csv.gz", 5041, "f2a4dc1997ed8ba2344094d534ff64fd",
      zenodo(19582556, "STalign_cell_type_transcriptional_correlations.csv.gz"), "Zenodo 19582556", "10.5281/zenodo.19582556",
      cc_by, stalign_creators),
    c("celltype", "STalign_S2R3_cell_type_annotations.csv.gz", 1244272, "085e49b9bb8f0693d5854d2013b351dc",
      zenodo(19582556, "STalign_S2R3_cell_type_annotations.csv.gz"), "Zenodo 19582556", "10.5281/zenodo.19582556",
      cc_by, stalign_creators),
    c("celltype", "STalign_Visium_cell_type_annotations.csv.gz", 264260, "897985b732e45a7b4dbee3a64cd357e8",
      zenodo(19582556, "STalign_Visium_cell_type_annotations.csv.gz"), "Zenodo 19582556", "10.5281/zenodo.19582556",
      cc_by, stalign_creators),
    c("merfish", "STalign_S2R2.csv.gz", 12517332, "e104115bb371a92efc06f7a3b908efc2",
      zenodo(10724029, "STalign_S2R2.csv.gz"), "Zenodo 10724029", "10.5281/zenodo.10724029", cc_by, stalign_creators),
    c("merfish", "STalign_S2R3_to_S2R2.csv.gz", 17110794, "73c1969dd2f52c12cbd4870b30f741d1",
      zenodo(10724029, "STalign_S2R3_to_S2R2.csv.gz"), "Zenodo 10724029", "10.5281/zenodo.10724029", cc_by, stalign_creators)
  )
  m <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(m) <- c("group", "file", "bytes", "md5", "url", "record", "doi", "licence", "creators")
  m$bytes <- as.numeric(m$bytes)
  m
})

# Attribution rows (CC BY 4.0: creator, source, licence) for the given downloaded files.
stc_sources <- function(files) {
  m <- stc_manifest[match(files, stc_manifest$file), c("file", "record", "doi", "url", "licence", "creators")]
  if (anyNA(m$file)) stop("Unknown file(s): ", paste(files[is.na(m$file)], collapse = ", "), call. = FALSE)
  rownames(m) <- NULL
  m
}

# Repository root: the data-raw scripts must be run from the directory that holds STcompare's DESCRIPTION.
stc_repo_root <- function() {
  ok <- file.exists("DESCRIPTION") &&
    any(grepl("^Package:[[:space:]]*STcompare[[:space:]]*$", readLines("DESCRIPTION", warn = FALSE)))
  if (!ok) stop("Run the data-raw scripts from the STcompare repository root (the directory containing DESCRIPTION).",
                call. = FALSE)
  normalizePath(".", mustWork = TRUE)
}

stc_cache_root <- function() {
  root <- Sys.getenv("STCOMPARE_DATA_CACHE", unset = "")
  if (!nzchar(root)) root <- tools::R_user_dir("STcompare", which = "cache")
  path.expand(root)
}

# Absolute form of `p` with symbolic links resolved in its deepest existing ancestor (the rest does not exist
# yet, so it cannot be a link).
stc_resolve_path <- function(p) {
  p <- path.expand(p)
  if (!grepl("^(/|[A-Za-z]:)", p)) p <- file.path(getwd(), p)
  rest <- character(0)
  while (!file.exists(p)) {
    parent <- dirname(p)
    if (identical(parent, p)) break
    rest <- c(basename(p), rest)
    p <- parent
  }
  do.call(file.path, as.list(c(normalizePath(p, mustWork = FALSE), rest)))
}

stc_inside_repo <- function(p, repo) {
  a <- paste0(p, "/")
  b <- paste0(repo, "/")
  if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") { a <- tolower(a); b <- tolower(b) }
  startsWith(a, b)
}

# <root>/data-raw/<...>, created if needed. Refuses any location inside the repository, which lives in
# a synced folder and must never receive downloads, caches or rasterized intermediates. The check runs on the
# symlink-resolved path BEFORE anything is created, and again on the created directory.
stc_cache_dir <- function(...) {
  d <- file.path(stc_cache_root(), "data-raw", ...)
  repo <- stc_repo_root()
  target <- stc_resolve_path(d)
  if (stc_inside_repo(target, repo)) {
    stop("Refusing to use a cache directory inside the repository: ", d, if (!identical(target, d)) paste0(" (", target, ")"),
         "\nSet STCOMPARE_DATA_CACHE to a directory outside the repository.", call. = FALSE)
  }
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d <- normalizePath(d, mustWork = TRUE)
  if (stc_inside_repo(d, repo)) stop("Refusing to use a cache directory inside the repository: ", d, call. = FALSE)
  d
}

stc_download_dir <- function() stc_cache_dir()
stc_inputs_dir <- function() stc_cache_dir("inputs")
stc_input_file <- function(name) file.path(stc_inputs_dir(), name)

# The authors' published results (.RData), kept in the repository under bench/published/ (they shipped in
# inst/extdata until version 0.1.0; they are not part of the package).
stc_published_file <- function(...) file.path("bench", "published", ...)

stc_md5 <- function(path) unname(tools::md5sum(path))

# Download (if needed) and verify the requested files; returns their paths, named by file.
stc_download <- function(group = NULL, files = NULL, quiet = FALSE) {
  m <- stc_manifest
  if (!is.null(group)) m <- m[m$group %in% group, , drop = FALSE]
  if (!is.null(files)) m <- m[m$file %in% files, , drop = FALSE]
  if (!nrow(m)) stop("No manifest entries match group = ", deparse(group), ", files = ", deparse(files), call. = FALSE)
  dir <- stc_download_dir()
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  paths <- setNames(file.path(dir, m$file), m$file)
  for (k in seq_len(nrow(m))) {
    f <- m$file[k]; dest <- paths[[f]]; want <- m$md5[k]
    if (file.exists(dest)) {
      got <- stc_md5(dest)
      if (!identical(got, want)) {
        stop(sprintf(paste0("Checksum mismatch for cached file\n  %s\n  expected md5 %s, found %s (%s bytes; expected %s).\n",
                            "The file is corrupt or the upstream file changed. Delete it to download it again from\n  %s"),
                     dest, want, got, format(file.size(dest)), format(m$bytes[k]), m$url[k]), call. = FALSE)
      }
      if (!quiet) message(sprintf("[cached, md5 ok] %s", f))
      next
    }
    if (!quiet) message(sprintf("[download] %s (%.1f MB) from %s", f, m$bytes[k] / 2^20, m$url[k]))
    part <- paste0(dest, ".part")
    status <- tryCatch(utils::download.file(m$url[k], part, mode = "wb", quiet = quiet), error = function(e) e)
    if (inherits(status, "error") || !identical(as.integer(status), 0L)) {
      unlink(part)
      stop(sprintf("Download failed for %s\n  URL: %s\n  %s", f, m$url[k],
                   if (inherits(status, "error")) conditionMessage(status) else paste("status", status)), call. = FALSE)
    }
    got <- stc_md5(part)
    if (!identical(got, want)) {
      bad <- paste0(dest, ".md5-mismatch")
      file.rename(part, bad)
      stop(sprintf(paste0("Checksum mismatch after download of %s\n  URL: %s\n  expected md5 %s, got %s (%s bytes; expected %s).\n",
                          "The download was kept as %s for inspection; the upstream file may have changed."),
                   f, m$url[k], want, got, format(file.size(bad)), format(m$bytes[k]), bad), call. = FALSE)
    }
    if (!file.rename(part, dest)) stop("Could not rename ", part, " to ", dest, call. = FALSE)
    if (!quiet) message(sprintf("[download] %s: md5 ok", f))
  }
  invisible(paths)
}

# Small helpers shared by the build scripts ---------------------------------------------------------

stc_md5_files <- function(paths) {
  paths <- paths[file.exists(paths)]
  setNames(unname(tools::md5sum(paths)), paths)
}

# Git state of the repository, trusted only when the enclosing git work tree IS the repository (not a parent
# repository); NA when git or the repository is unavailable. `status` lists every changed or untracked file
# under the paths that make up the package code and the builders.
stc_git_info <- function(repo = stc_repo_root()) {
  na <- list(available = nzchar(Sys.which("git")), in_repository = NA, commit = NA_character_, dirty = NA,
             status = NA_character_)
  if (!na$available) return(na)
  run <- function(...) {
    out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(repo), "--no-optional-locks", ...),
                                             stdout = TRUE, stderr = FALSE)),
                    error = function(e) structure(character(0), status = 1L))
    st <- attr(out, "status")
    if (!is.null(st) && st != 0) NULL else out
  }
  top <- run("rev-parse", "--show-toplevel")
  same <- function(a, b) {
    a <- normalizePath(a, mustWork = FALSE); b <- normalizePath(b, mustWork = FALSE)
    if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") identical(tolower(a), tolower(b)) else identical(a, b)
  }
  if (length(top) != 1 || !same(top, repo)) { na$in_repository <- FALSE; return(na) }
  commit <- run("rev-parse", "HEAD")
  status <- run("status", "--porcelain", "--untracked-files=all", "--", "R", "NAMESPACE", "DESCRIPTION", "data",
                "data-raw", "tests")
  list(available = TRUE, in_repository = TRUE, commit = if (length(commit)) commit[1] else NA_character_,
       dirty = if (is.null(status)) NA else length(status) > 0, status = if (is.null(status)) NA_character_ else status)
}

# Provenance recorded inside every rasterized input and fixture. `goldens`: published result files
# (bench/published/, formerly inst/extdata/) the output was derived from; `sources`: attribution rows from
# stc_sources().
stc_meta <- function(script, inputs = NULL, sources = NULL, goldens = NULL, extra = list()) {
  repo <- stc_repo_root()
  pv <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  si <- tryCatch(utils::sessionInfo(), error = function(e) NULL)
  c(list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), script = script,
         R = R.version.string, platform = R.version$platform,
         BLAS = if (is.null(si)) NA_character_ else si$BLAS, LAPACK = if (is.null(si)) NA_character_ else si$LAPACK,
         long_double = unname(capabilities("long.double")),
         RNGkind = RNGkind(),
         STcompare_version = tryCatch(as.character(read.dcf("DESCRIPTION", "Version")[1, 1]), error = function(e) NA_character_),
         git = stc_git_info(repo),
         # content hashes identify the code that made the file even when it is not committed
         content_md5 = list(package = stc_md5_files(c("DESCRIPTION", "NAMESPACE", sort(list.files("R", full.names = TRUE)))),
                            data = stc_md5_files(sort(list.files("data", full.names = TRUE))),
                            data_raw = stc_md5_files(sort(list.files("data-raw", pattern = "[.]R$", full.names = TRUE))),
                            test_helpers = stc_md5_files(sort(list.files(file.path("tests", "testthat"),
                                                                         pattern = "^helper.*[.]R$", full.names = TRUE))),
                            goldens = stc_md5_files(goldens)),
         packages = vapply(c("locfit", "geoR", "SEraster", "SpatialExperiment", "SummarizedExperiment", "BiocParallel",
                             "Matrix", "rhdf5", "sf", "rearrr", "testthat"), pv, ""),
         sf_ext_soft = tryCatch(sf::sf_extSoftVersion(), error = function(e) NA_character_),
         inputs = inputs, sources = sources),
    extra)
}

stc_save_rds <- function(x, path, compress = "gzip") {
  base::saveRDS(x, path, compress = compress)
  message(sprintf("wrote %s (%.1f MB)", path, file.size(path) / 2^20))
  invisible(path)
}

if (sys.nframe() == 0L) {
  stc_repo_root()
  invisible(stc_download())
  message("All downloads present in ", stc_download_dir())
}
