#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

  if (length(file_arg) == 0) {
    stop("This script is intended to be run with Rscript.")
  }

  normalizePath(sub("^--file=", "", file_arg[1]))
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

source(file.path(repo_root, "R", "spatialCorrelation.R"))
source(file.path(repo_root, "R", "iterativePermutations.R"))

resolutions <- c(50, 100, 200, 250, 300, 500)
replicate_runs <- 5L
n_permutations <- c(100L, 1000L)
n_threads <- 22L
seed_base <- 20260423L
delta_grid <- c(0.01, 0.05, seq(0.1, 0.9, by = 0.1))
return_permutations <- FALSE

cache_dir <- file.path(script_dir, "cache")
output_dir <- file.path(script_dir, "outputs")
log_dir <- file.path(script_dir, "logs")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

target_url <- "https://zenodo.org/records/10724029/files/STalign_S2R2.csv.gz?download=1"
source_url <- "https://zenodo.org/records/10724029/files/STalign_S2R3_to_S2R2.csv.gz?download=1"

target_path <- file.path(cache_dir, "STalign_S2R2.csv.gz")
source_path <- file.path(cache_dir, "STalign_S2R3_to_S2R2.csv.gz")

download_if_missing <- function(url, destfile) {
  if (!file.exists(destfile)) {
    download.file(url = url, destfile = destfile, mode = "wb", quiet = FALSE)
  }
}

read_csv_gz <- function(path, ...) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.csv(con, ...)
}

get_uptime <- function() {
  tryCatch(
    paste(system2("uptime", stdout = TRUE, stderr = TRUE), collapse = "\n"),
    error = function(e) {
      paste("uptime unavailable:", conditionMessage(e))
    }
  )
}

log_stamp <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", usetz = TRUE)),
    ...,
    "\n",
    sep = ""
  )
}

with_log_capture <- function(log_path, expr) {
  log_con <- file(log_path, open = "wt")
  on.exit(close(log_con), add = TRUE)

  sink(log_con, split = TRUE)
  on.exit(sink(), add = TRUE)

  sink(log_con, type = "message")
  on.exit(sink(type = "message"), add = TRUE)

  force(expr)
}

prepare_merfish_inputs <- function(target_csv, source_csv) {
  target <- read_csv_gz(target_csv, check.names = FALSE)
  source <- read_csv_gz(source_csv, check.names = FALSE)

  pos_target <- target[, c("x", "y"), drop = FALSE]
  rownames(pos_target) <- target$X

  gene_target <- target[, 4:ncol(target), drop = FALSE]
  rownames(gene_target) <- target$X
  gene_target <- gene_target[, !grepl("^Blank\\.", colnames(gene_target)), drop = FALSE]

  target_sparse <- Matrix::Matrix(t(as.matrix(gene_target)), sparse = TRUE)
  spe_target <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = target_sparse),
    spatialCoords = as.matrix(pos_target)
  )

  pos_source <- source[, c("STalign_x", "STalign_y"), drop = FALSE]
  rownames(pos_source) <- source$X
  colnames(pos_source) <- c("x", "y")

  gene_source <- source[, 8:ncol(source), drop = FALSE]
  rownames(gene_source) <- source$X
  gene_source <- gene_source[, !grepl("^Blank\\.", colnames(gene_source)), drop = FALSE]

  source_sparse <- Matrix::Matrix(t(as.matrix(gene_source)), sparse = TRUE)
  spe_source <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = source_sparse),
    spatialCoords = as.matrix(pos_source)
  )

  list(target = spe_target, source = spe_source)
}

infer_n_permutations_used <- function(corr_output) {
  if (!"deltaStarX" %in% colnames(corr_output)) {
    stop("corr_output must contain a deltaStarX column.")
  }

  vapply(corr_output$deltaStarX, function(x) {
    if (is.null(x)) {
      return(NA_integer_)
    }

    if (length(x) == 1L && is.atomic(x) && all(is.na(x))) {
      return(NA_integer_)
    }

    as.integer(length(x))
  }, integer(1))
}

count_genes_by_iteration <- function(final_n_perm, n_perm_vector, total_genes) {
  counts <- vapply(seq_along(n_perm_vector), function(i) {
    if (i == 1L) {
      return(as.integer(total_genes))
    }

    as.integer(sum(!is.na(final_n_perm) & final_n_perm >= n_perm_vector[i]))
  }, integer(1))

  stats::setNames(
    as.list(counts),
    paste0("genes_npermutation_", n_perm_vector)
  )
}

save_named_output <- function(object_name, object_value, output_path) {
  save_env <- new.env(parent = emptyenv())
  assign(object_name, object_value, envir = save_env)
  save(list = object_name, file = output_path, envir = save_env)
}

run_single_benchmark <- function(resolution, run_number, spe_list) {
  object_name <- sprintf("Merfish_res%s_run%s", resolution, run_number)
  log_path <- file.path(log_dir, paste0(object_name, ".log"))
  output_path <- file.path(output_dir, paste0(object_name, ".RData"))

  start_time <- Sys.time()

  summary_row <- with_log_capture(log_path, {
    log_stamp("Starting ", object_name)
    log_stamp("Uptime at start: ", get_uptime())
    log_stamp("Resolution: ", resolution)
    log_stamp("Replicate run: ", run_number)
    log_stamp("nPermutations: ", paste(n_permutations, collapse = ", "))
    log_stamp("nThreads: ", n_threads)

    status <- "success"
    error_message <- NA_character_
    pixel_count <- NA_integer_
    source_pixel_count <- NA_integer_
    target_pixel_count <- NA_integer_
    shared_pixel_count <- NA_integer_
    per_iteration_counts <- stats::setNames(
      as.list(rep(NA_integer_, length(n_permutations))),
      paste0("genes_npermutation_", n_permutations)
    )

    tryCatch({
      raster_output <- SEraster::rasterizeGeneExpression(
        spe_list,
        resolution = resolution,
        BPPARAM = BiocParallel::SerialParam()
      )

      source_pixel_count <- ncol(raster_output$source)
      target_pixel_count <- ncol(raster_output$target)
      shared_pixels <- intersect(
        rownames(SpatialExperiment::spatialCoords(raster_output$source)),
        rownames(SpatialExperiment::spatialCoords(raster_output$target))
      )
      shared_pixel_count <- length(shared_pixels)
      pixel_count <- shared_pixel_count

      log_stamp("Source pixels: ", source_pixel_count)
      log_stamp("Target pixels: ", target_pixel_count)
      log_stamp("Shared pixels: ", shared_pixel_count)

      delta_list <- rep(list(delta_grid), length(rownames(raster_output$target)))
      run_seed <- as.integer(seed_base + (resolution * 100L) + run_number)

      log_stamp("Using seed: ", run_seed)

      corr_output <- spatialCorrelationGeneExpIterPermutations(
        input = raster_output,
        nPermutations = n_permutations,
        deltaX = delta_list,
        deltaY = delta_list,
        returnPermutations = return_permutations,
        nThreads = n_threads,
        verbose = TRUE,
        seed = run_seed
      )

      corr_output$nPermutationsUsed <- infer_n_permutations_used(corr_output)
      per_iteration_counts <- count_genes_by_iteration(
        final_n_perm = corr_output$nPermutationsUsed,
        n_perm_vector = n_permutations,
        total_genes = nrow(corr_output)
      )

      save_named_output(object_name, corr_output, output_path)
      log_stamp("Saved corr_output to: ", output_path)
      log_stamp(
        "Gene counts by permutation level: ",
        paste(
          names(per_iteration_counts),
          unlist(per_iteration_counts, use.names = FALSE),
          sep = "=",
          collapse = ", "
        )
      )
    }, error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      log_stamp("Run failed: ", error_message)
    })

    end_time <- Sys.time()
    elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))

    log_stamp("Uptime at end: ", get_uptime())
    log_stamp("End time: ", format(end_time, "%Y-%m-%d %H:%M:%S %Z", usetz = TRUE))
    log_stamp("Elapsed seconds: ", sprintf("%.3f", elapsed_seconds))
    log_stamp("Status: ", status)

    c(
      list(
        resolution = resolution,
        number_of_pixels = pixel_count,
        source_pixels = source_pixel_count,
        target_pixels = target_pixel_count,
        shared_pixels = shared_pixel_count,
        start_time = format(start_time, "%Y-%m-%d %H:%M:%S %Z", usetz = TRUE),
        end_time = format(end_time, "%Y-%m-%d %H:%M:%S %Z", usetz = TRUE),
        total_time_seconds = elapsed_seconds,
        replicate_run = run_number,
        log_file = normalizePath(log_path, winslash = "/", mustWork = FALSE),
        corr_output_file = normalizePath(output_path, winslash = "/", mustWork = FALSE),
        corr_output_object = object_name,
        status = status,
        error_message = error_message
      ),
      per_iteration_counts
    )
  })

  as.data.frame(summary_row, stringsAsFactors = FALSE, check.names = FALSE)
}

download_if_missing(target_url, target_path)
download_if_missing(source_url, source_path)

spe_list <- prepare_merfish_inputs(target_path, source_path)

summary_rows <- vector("list", length(resolutions) * replicate_runs)
row_index <- 1L

for (resolution in resolutions) {
  for (run_number in seq_len(replicate_runs)) {
    summary_rows[[row_index]] <- run_single_benchmark(
      resolution = resolution,
      run_number = run_number,
      spe_list = spe_list
    )
    row_index <- row_index + 1L
  }
}

benchmark_results <- do.call(rbind, summary_rows)
rownames(benchmark_results) <- NULL

summary_csv <- file.path(output_dir, "Merfish_resolution_benchmark_summary.csv")
summary_rds <- file.path(output_dir, "Merfish_resolution_benchmark_summary.rds")

utils::write.csv(benchmark_results, summary_csv, row.names = FALSE)
saveRDS(benchmark_results, summary_rds)

print(benchmark_results)
