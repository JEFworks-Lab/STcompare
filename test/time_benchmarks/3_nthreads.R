# Test at different threads ------------------------------------------------------------------


devtools::load_all()

load("test/time_benchmarks/inputs/spe_list.Rdata")

nThreads <- seq(2, 22, by = 2)

# interesting to see how this changes as a result of change in resolution 
resolutions <- c(600, 400, 200, 100)

n_runs <- 5
max_threads <- max(nThreads)

dir.create("test/time_benchmarks/outputs/3_nthreads/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("test/time_benchmarks/outputs/3_nthreads/corr_output", recursive = TRUE, showWarnings = FALSE)

summary_file <- "test/time_benchmarks/outputs/3_nthreads/nthreads_benchmark_summary.csv"

get_uptime <- function() {
  paste(system("uptime", intern = TRUE), collapse = " ")
}

append_csv <- function(df, file) {
  if (!file.exists(file)) {
    write.csv(df, file, row.names = FALSE)
  } else {
    write.table(
      df,
      file = file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

count_genes_by_nperm <- function(corr_output) {
  if (is.null(corr_output) || nrow(corr_output) == 0) {
    return(list(
      n_genes_perm100 = NA_integer_,
      n_genes_perm1000 = NA_integer_
    ))
  }

  n_perms <- vapply(
    corr_output$deltaStarX,
    function(x) {
      if (length(x) == 1 && is.atomic(x) && all(is.na(x))) return(NA_integer_)
      length(x)
    },
    integer(1)
  )

  list(
    n_genes_perm100 = sum(n_perms == 100, na.rm = TRUE),
    n_genes_perm1000 = sum(n_perms == 1000, na.rm = TRUE)
  )
}

get_n_pixels <- function(rast) {
  source_spe <- rast[[1]]
  target_spe <- rast[[2]]

  shared_pixels <- intersect(
    rownames(SpatialExperiment::spatialCoords(source_spe)),
    rownames(SpatialExperiment::spatialCoords(target_spe))
  )

  length(shared_pixels)
}

run_benchmark <- function(rast, resolution, n_thread, run_num, seed = 0) {
  var_name <- paste0(
    "merfish_res", resolution,
    "_nthreads", n_thread,
    "_run_", run_num
  )
  log_file <- file.path("test/time_benchmarks/outputs/3_nthreads/logs", paste0(var_name, ".log"))
  out_file <- file.path("test/time_benchmarks/outputs/3_nthreads/corr_output", paste0(var_name, ".RData"))

  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  status <- "success"
  error_message <- ""
  n_pixels <- get_n_pixels(rast)
  start_time <- Sys.time()
  end_time <- start_time

  cat("=====================================================\n")
  cat("Benchmark run:", var_name, "\n")
  cat("Resolution:", resolution, "\n")
  cat("nThreads:", n_thread, "\n")
  cat("Run number:", run_num, "\n")
  cat("Start time:", as.character(start_time), "\n")
  cat("Uptime at start:", get_uptime(), "\n")
  cat("Number of shared pixels:", n_pixels, "\n")
  cat("=====================================================\n\n")

  result <- tryCatch({
    corr_output <- spatialCorrelationGeneExpIterPermutations(
      rast,
      nThreads = n_thread,
      returnPermutations = FALSE,
      seed = seed,
      BPPARAM = BiocParallel::MulticoreParam(workers = n_thread)
    )
    end_time <<- Sys.time()

    gene_counts <- count_genes_by_nperm(corr_output)

    assign(var_name, corr_output, envir = .GlobalEnv)
    save(list = var_name, file = out_file, envir = .GlobalEnv)

    cat("Saved corr_output as:", out_file, "\n")
    cat("Object name:", var_name, "\n")

    list(
      n_genes_perm100 = gene_counts$n_genes_perm100,
      n_genes_perm1000 = gene_counts$n_genes_perm1000
    )
  }, error = function(e) {
    status <<- "error"
    error_message <<- conditionMessage(e)
    end_time <<- Sys.time()

    cat("\nERROR:\n")
    cat(error_message, "\n")

    list(
      n_genes_perm100 = NA_integer_,
      n_genes_perm1000 = NA_integer_
    )
  })

  total_time_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

  cat("\n=====================================================\n")
  cat("End time:", as.character(end_time), "\n")
  cat("Total time (clock time, sec):", total_time_sec, "\n")
  cat("Uptime at end:", get_uptime(), "\n")
  cat("Status:", status, "\n")
  if (error_message != "") {
    cat("Error message:", error_message, "\n")
  }
  cat("=====================================================\n")

  data.frame(
    resolution = resolution,
    n_threads = n_thread,
    n_pixels = n_pixels,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    total_time_sec = total_time_sec,
    n_genes_perm100 = result$n_genes_perm100,
    n_genes_perm1000 = result$n_genes_perm1000,
    replicate_run_number = run_num,
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

all_results <- list()
idx <- 1

for (res in resolutions) {
  cat("Rasterizing resolution", res, "with", max_threads, "threads\n")
  cat("Rasterizing...\n")
  rast_start_time <- Sys.time()
  rast <- SEraster::rasterizeGeneExpression(
    spe_list,
    resolution = res,
    BPPARAM = BiocParallel::MulticoreParam(workers = max_threads)
  )
  rast_end_time <- Sys.time()
  cat("Rast time diff: ", (rast_end_time - rast_start_time), "\n")

  for (run_num in seq_len(n_runs)) {
    for (n_thread in nThreads) {
      row <- run_benchmark(
        rast = rast,
        resolution = res,
        n_thread = n_thread,
        run_num = run_num,
        seed = run_num
      )
      all_results[[idx]] <- row
      idx <- idx + 1
      append_csv(row, summary_file)
    }
  }
}

benchmark_summary <- do.call(rbind, all_results)
benchmark_summary




