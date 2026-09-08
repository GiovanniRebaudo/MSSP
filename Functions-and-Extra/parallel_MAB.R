# Parallelize complete replicas, never the MCMC steps within a replica.
# Each replica function must retain its original explicit seeds and call order.
.mab_replica_with_rng <- function(seed, replica_fun, ...){
  value = replica_fun(seed, ...)
  list(value = value, rng_state = get(".Random.seed", envir = globalenv()))
}

run_mab_replicas <- function(seeds, replica_fun, ..., workers = 4L,
                             code_dir = getwd()){
  stopifnot(length(workers) == 1L, is.finite(workers),
            workers >= 1, workers == as.integer(workers))
  if(length(seeds) == 0L) return(list())
  workers = min(as.integer(workers), length(seeds))
  if(workers == 1L) return(lapply(seeds, replica_fun, ...))

  code_dir = normalizePath(code_dir, mustWork = TRUE)
  cluster = parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  # Use the same R libraries and RNG kind; do not replace the existing seeds
  # with parallel RNG streams, which would change the numerical results.
  parallel::clusterCall(cluster, function(code_dir, library_paths, rng_kind){
    .libPaths(library_paths)
    do.call(RNGkind, as.list(rng_kind))
    setwd(code_dir)
    source("mSSPmab.R")
    NULL
  }, code_dir, .libPaths(), RNGkind())

  message("Running ", length(seeds), " replicas on ", workers, " workers...")
  results = parallel::parLapply(cluster, seeds, .mab_replica_with_rng,
                               replica_fun = replica_fun, ...)
  # parLapply returns results in input order, regardless of completion order.
  # Also leave the caller's RNG in the state of the last sequential replica.
  assign(".Random.seed", results[[length(results)]]$rng_state, envir = globalenv())
  message("Completed ", length(seeds), " replicas.")
  lapply(results, `[[`, "value")
}
