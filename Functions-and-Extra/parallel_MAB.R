# Parallelize complete replicas, never the MCMC steps within a replica.
# Each replica function must retain its original explicit seeds and call order.
.mab_replica_with_rng <- function(seed, replica_fun, ...){
  value = replica_fun(seed, ...)
  list(value = value, rng_state = get(".Random.seed", envir = globalenv()))
}

run_mab_replicas <- function(seeds, replica_fun, ..., workers = 4L,
                             code_dir = getwd(), progress = TRUE){
  stopifnot(length(workers) == 1L, is.finite(workers),
            workers >= 1, workers == as.integer(workers),
            is.logical(progress), length(progress) == 1L, !is.na(progress))
  if(length(seeds) == 0L) return(list())
  workers = min(as.integer(workers), length(seeds))
  started = function(i){
    if(progress) message("Started replica ", i, "/", length(seeds),
                         " (seed ", seeds[[i]], ").")
  }
  finished = function(i, count){
    if(progress) message("Completed replica ", i, "/", length(seeds),
                         " (seed ", seeds[[i]], "); finished ", count,
                         "/", length(seeds), ".")
  }
  if(workers == 1L){
    return(lapply(seq_along(seeds), function(i){
      started(i)
      value = replica_fun(seeds[[i]], ...)
      finished(i, i)
      value
    }))
  }

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

  if(progress) message("Running ", length(seeds), " replicas on ", workers, " workers...")
  # Use the same dispatch/receive primitives as parallel's clusterApplyLB, so
  # progress is printed by the main R session (also visible in RStudio).
  # Each free worker receives one whole replica; its internal call order is unchanged.
  args = list(...)
  submit = function(node, job){
    parallel:::sendCall(cluster[[node]], .mab_replica_with_rng,
                       c(list(seed = seeds[[job]], replica_fun = replica_fun), args),
                       tag = job)
    started(job)
  }
  for(i in seq_len(workers)) submit(i, i)
  results = vector("list", length(seeds))
  for(count in seq_along(seeds)){
    completed = parallel:::recvOneResult(cluster)
    if(inherits(completed$value, "try-error"))
      stop("Replica ", completed$tag, " (seed ", seeds[[completed$tag]],
           ") failed: ", completed$value, call. = FALSE)
    results[[completed$tag]] = completed$value
    finished(completed$tag, count)
    next_job = workers + count
    if(next_job <= length(seeds)) submit(completed$node, next_job)
  }
  # Results are restored to input order, regardless of completion order.
  # Also leave the caller's RNG in the state of the last sequential replica.
  assign(".Random.seed", results[[length(results)]]$rng_state, envir = globalenv())
  if(progress) message("Completed ", length(seeds), " replicas.")
  lapply(results, `[[`, "value")
}
