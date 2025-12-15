setup_parallelisation <- function(reserve_cores = 4, verbose = TRUE) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel package is required.")
  }
  
  n_total <- parallel::detectCores()
  n_workers <- max(1, n_total - reserve_cores)
  
  if (verbose) {
    message("Detected cores: ", n_total)
    message("Using workers:  ", n_workers)
  }
  
  if (.Platform$OS.type == "unix") {
    param <- BiocParallel::MulticoreParam(workers = n_workers)
  } else {
    param <- BiocParallel::SnowParam(workers = n_workers, type = "SOCK")
  }
  
  BiocParallel::register(param)
  
  invisible(param)
}
