# On the sensitivities to the modifiable areal unit problem
# Raster sensitivity calculation functions
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# yexiang@nnu.edu.cn
# 2026-11-12

# Please cite the following reference when part or all of the code in this file
# is reused under the license of CC-BY-4.0:
# Ye, X., & Chen, J. (2026). On the sensitivities to the modifiable areal unit
# problem. Big Earth Data, 10(3), 1539-1574.
# https://doi.org/10.1080/20964471.2026.2692263

# -- Shared Monte Carlo execution core ---------------------------------------
# Keep this block identical in the dated vector and raster calculation modules.
# Simulation validity and sample collection remain the responsibility of callers.

default_parallel_workers <- function(reserve_cores = 10L) {
  # Choose a worker count from the detected logical cores, leaving the
  # requested number free. Fall back to one worker if detection is unavailable.

  if (!is.numeric(reserve_cores) || length(reserve_cores) != 1L ||
      !is.finite(reserve_cores) || reserve_cores < 0 ||
      reserve_cores != floor(reserve_cores)) {
    stop("reserve_cores must be one non-negative integer")
  }
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || cores < 2L) return(1L)
  as.integer(max(1, cores - reserve_cores))
}

restore_monte_carlo_rng <- function(kind, seed, had_seed) {
  # Restore the calling process's random-number settings after seeded work.
  # kind stores RNGkind(); seed is the previous .Random.seed. Restore the
  # generator first, since changing it can also change the random state.
  do.call(RNGkind, as.list(kind))
  if (had_seed) {
    assign(".Random.seed", seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    # No state existed before the operation: remove the newly created one.
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

# -- Monte Carlo execution lifecycle -----------------------------------------
# Create a record -> submit batches with run_monte_carlo -> close the record.
# Workers start only when needed. The sensitivity algorithm, not these
# helpers, decides which paths are valid and whether to submit more tasks.

create_monte_carlo_executor <- function(FUN,
                                        parallel = FALSE,
                                        n_workers = default_parallel_workers(),
                                        random_seed = NULL,
                                        export_names = character(),
                                        export_env = .GlobalEnv,
                                        packages = character(),
                                        progress_label = "Monte Carlo") {
  # Prepare a reusable record containing the task and its execution settings.
  # FUN performs one complete simulation; this function does not run it or
  # start worker processes. The returned environment is reused across batches.

  # -- Validate the execution configuration --------------------------------
  if (!is.function(FUN)) stop("FUN must be a function")
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE")
  }
  if (!is.numeric(n_workers) || length(n_workers) != 1L ||
      !is.finite(n_workers) || n_workers < 1 || n_workers != floor(n_workers) ||
      n_workers > .Machine$integer.max) {
    stop("n_workers must be one positive integer")
  }
  if (!is.null(random_seed) &&
      (!is.numeric(random_seed) || length(random_seed) != 1L ||
       !is.finite(random_seed) || random_seed < 0 ||
       random_seed != floor(random_seed) || random_seed > .Machine$integer.max)) {
    stop("random_seed must be NULL or one non-negative integer")
  }
  for (name in c("export_names", "packages")) {
    value <- get(name)
    if (!is.character(value) || anyNA(value) || any(!nzchar(value))) {
      stop(sprintf("%s must be a character vector of non-empty names", name))
    }
  }
  if (!is.environment(export_env)) stop("export_env must be an environment")
  if (!is.character(progress_label) || length(progress_label) != 1L ||
      is.na(progress_label) || !nzchar(progress_label)) {
    stop("progress_label must be one non-empty string")
  }

  # -- Resolve dependencies before starting workers ------------------------
  # The caller supplies package dependencies; this shared core does not
  # assume a vector or raster representation. No packages are loaded by default.
  packages <- unique(packages)
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop(sprintf("Missing requested package(s): %s", paste(missing_packages, collapse = ", ")))
  }
  # Look up the named helper functions or data in export_env and its parents.
  # Workers will receive these values under the same names.
  export_names <- unique(export_names)
  values <- setNames(lapply(export_names, function(name) {
    if (!exists(name, envir = export_env, inherits = TRUE)) {
      stop(sprintf("Export object '%s' was not found", name))
    }
    get(name, envir = export_env, inherits = TRUE)
  }), export_names)

  # -- Capture the task before adding live worker connections ---------------
  # Serializing and reading back copies local function environments and
  # exported objects into a snapshot for worker initialisation.
  # Freeze local closure environments before the caller stores this executor
  # beside FUN. Only this snapshot is sent to workers, once per pool.
  frozen <- unserialize(serialize(list(FUN = FUN, exports = values), NULL))
  # An environment is updated by reference, so later batches see the same
  # connections, submitted-task count and next random stream.
  executor <- new.env(parent = emptyenv())
  class(executor) <- "monte_carlo_executor"
  # Keep the original settings to detect conflicting arguments on reuse.
  executor$config <- list(FUN = FUN, parallel = parallel, n_workers = n_workers,
                          random_seed = random_seed, export_names = export_names,
                          export_env = export_env, packages = packages,
                          progress_label = progress_label)
  executor$FUN <- frozen$FUN
  executor$exports <- frozen$exports
  # Mutable state: cluster stores worker connections; closed blocks reuse
  # after shutdown; n_submitted counts tasks across batches; next_seed holds
  # the next task's RNG state (NULL means no explicit stream is assigned).
  executor$cluster <- NULL
  executor$closed <- FALSE
  executor$n_submitted <- 0
  executor$next_seed <- NULL

  # -- Prepare the first task's random stream ------------------------------
  # L'Ecuyer-CMRG supports separate streams for successive tasks. Preserve
  # the caller's RNG settings with on.exit(), including on an early error.
  if (!is.null(random_seed)) {
    kind <- RNGkind()
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
    on.exit(restore_monte_carlo_rng(kind, seed, had_seed), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(as.integer(random_seed))
    executor$next_seed <- get(".Random.seed", .GlobalEnv)
  }
  executor
}

close_monte_carlo_executor <- function(executor) {
  # Release the worker processes. The creator registers this with on.exit()
  # so normal completion and errors both trigger cleanup.

  if (!inherits(executor, "monte_carlo_executor") || !is.environment(executor)) {
    stop("executor must be a Monte Carlo executor")
  }
  # Repeated cleanup is harmless. Mark the executor closed before shutdown
  # so later submissions cannot reuse these connections.
  if (executor$closed) return(invisible(NULL))
  executor$closed <- TRUE
  cl <- executor$cluster
  executor$cluster <- NULL
  # A broken connection can also fail during cleanup. Do not let that
  # secondary error hide the error that originally ended the calculation.
  if (!is.null(cl)) tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
  invisible(NULL)
}

start_monte_carlo_workers <- function(executor, n_tasks) {
  # Start a PSOCK pool when a batch needs parallel execution. Reuse an
  # existing pool without resending its task function and dependencies.

  # -- Reuse the existing pool or leave execution serial -------------------
  if (!is.null(executor$cluster) || !executor$config$parallel ||
      executor$config$n_workers <= 1L || n_tasks <= 1L) return(invisible(NULL))

  # -- Launch independent R processes -------------------------------------
  # Cap the initial process count at the number of submitted tasks.
  cl <- parallel::makePSOCKcluster(min(executor$config$n_workers, n_tasks))
  # Until initialisation succeeds, this function owns cleanup of the new
  # pool. Afterwards, the executor's creator is responsible for closing it.
  initialized <- FALSE
  on.exit(if (!initialized) {
    tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
  }, add = TRUE)

  # -- Initialise each worker once ----------------------------------------
  # This function runs inside each worker: .GlobalEnv below belongs to that
  # worker, where the packages, exported objects and task function are kept.
  initialize <- function(FUN, exports, packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Missing requested worker package: %s", pkg))
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
    list2env(exports, envir = .GlobalEnv)
    assign(".maup_monte_carlo_FUN", FUN, envir = .GlobalEnv)
    NULL
  }
  # Only base functions and explicit arguments are needed. Using baseenv()
  # avoids sending this function's enclosing coordinator environment.
  environment(initialize) <- baseenv()
  # clusterCall runs the same initialisation on every worker in the pool.
  parallel::clusterCall(cl, initialize, executor$FUN, executor$exports,
                        executor$config$packages)
  executor$cluster <- cl
  initialized <- TRUE
  invisible(NULL)
}

dispatch_monte_carlo_task <- function(task, FUN = NULL) {
  # Execute one complete simulation. task contains value (the argument for
  # FUN), ordinal (its submission number across batches), and seed (its
  # assigned random stream, or NULL). Serial calls supply FUN directly;
  # parallel calls retrieve the function previously installed on the worker.

  # This function's environment is reset to baseenv before dispatch: no
  # coordinator, executor, or cluster is serialized with individual tasks.
  if (is.null(FUN)) FUN <- get(".maup_monte_carlo_FUN", envir = .GlobalEnv)
  # Select this task's stream independently of which worker receives it.
  if (!is.null(task$seed)) assign(".Random.seed", task$seed, envir = .GlobalEnv)
  # ok indicates whether FUN returned normally, not whether a simulated
  # path is valid. A returned list(error = ...) is kept as a normal value;
  # only an uncaught R error is wrapped with ok = FALSE and the task number.
  tryCatch(list(ok = TRUE, value = FUN(task$value)),
           error = function(e) list(ok = FALSE, ordinal = task$ordinal,
                                    message = conditionMessage(e)))
}

run_monte_carlo <- function(X,
                            FUN,
                            parallel = FALSE,
                            n_workers = default_parallel_workers(),
                            random_seed = NULL,
                            export_names = character(0),
                            export_env = .GlobalEnv,
                            progress_label = "Monte Carlo",
                            executor = NULL,
                            progress = TRUE,
                            chunk_size = NULL,
                            packages = character()) {
  # Execute one submission of independent tasks and return a list in input
  # order. Each element of X is passed to FUN; no results are filtered here.

  # -- Validate controls for this submission -------------------------------
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE")
  }
  if (!is.null(chunk_size) &&
      (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
       !is.finite(chunk_size) || chunk_size < 1 ||
       chunk_size != floor(chunk_size) || chunk_size > .Machine$integer.max)) {
    stop("chunk_size must be NULL or one positive integer")
  }

  # -- Create a temporary executor or borrow an existing one ---------------
  # Only an executor created here is closed on return. When reusing one
  # across retry batches, the calling sensitivity function owns cleanup.
  owned <- is.null(executor)
  if (owned) {
    if (missing(FUN)) stop("FUN must be provided when executor is NULL")
    executor <- create_monte_carlo_executor(
      FUN, parallel, n_workers, random_seed, export_names, export_env,
      packages, progress_label
    )
    on.exit(close_monte_carlo_executor(executor), add = TRUE)
  } else {
    if (!inherits(executor, "monte_carlo_executor") || !is.environment(executor)) {
      stop("executor must be a Monte Carlo executor")
    }
    # Explicit configuration is allowed only when it matches this session.
    # Inspect explicitly supplied arguments only: omitted defaults must not
    # overwrite the settings recorded when this executor was created.
    supplied <- as.list(match.call(expand.dots = FALSE))[-1L]
    config_names <- intersect(names(supplied), names(executor$config))
    for (name in config_names) {
      value <- get(name)
      expected <- executor$config[[name]]
      same <- if (name %in% c("n_workers", "random_seed")) {
        isTRUE(all.equal(value, expected))
      } else identical(value, expected)
      if (!same) stop(sprintf("'%s' conflicts with the supplied executor", name))
    }
  }
  if (executor$closed) stop("Cannot submit tasks to a closed executor")
  if (!length(X)) return(setNames(vector("list", 0L), names(X)))

  # -- Preserve the caller's random-number settings ------------------------
  # Explicit seeds never alter the caller's RNG state, including on errors.
  if (!is.null(executor$config$random_seed)) {
    kind <- RNGkind()
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
    on.exit(restore_monte_carlo_rng(kind, seed, had_seed), add = TRUE)
  }

  # -- Start workers if needed, then prepare tasks and random streams -------
  start_monte_carlo_workers(executor, length(X))
  # Number tasks and assign streams before dispatch. next_seed persists
  # across submissions, so a retry batch continues instead of restarting.
  tasks <- lapply(seq_along(X), function(i) {
    task <- list(value = X[[i]], ordinal = executor$n_submitted + i,
                 seed = executor$next_seed)
    if (!is.null(executor$next_seed)) {
      executor$next_seed <- parallel::nextRNGStream(executor$next_seed)
    }
    task
  })
  # This counts submitted tasks, including simulations that may later fail.
  executor$n_submitted <- executor$n_submitted + length(X)

  # -- Group tasks for progress reporting ---------------------------------
  # Aim for roughly ten updates, with at least one task per worker in each
  # full batch by default. The last batch may contain fewer tasks.
  if (is.null(chunk_size)) {
    workers <- if (is.null(executor$cluster)) 1L else length(executor$cluster)
    chunk_size <- max(workers, ceiling(length(X) / 10))
  }
  batches <- split(seq_along(X), ceiling(seq_along(X) / chunk_size))
  results <- setNames(vector("list", length(X)), names(X))
  dispatch <- dispatch_monte_carlo_task
  environment(dispatch) <- baseenv()

  # -- Execute batches and collect results --------------------------------
  # parLapply distributes one batch across workers and waits for all replies
  # before the next batch. lapply executes in the current R process.
  for (batch in batches) {
    replies <- if (is.null(executor$cluster)) {
      lapply(tasks[batch], function(task) dispatch(task, FUN = executor$FUN))
    } else {
      parallel::parLapply(executor$cluster, tasks[batch], dispatch)
    }
    for (j in seq_along(batch)) {
      reply <- replies[[j]]
      if (!reply$ok) {
        stop(sprintf("Monte Carlo task %s failed: %s", reply$ordinal, reply$message),
             call. = FALSE)
      }
      # Preserve NULL as a list entry instead of deleting its position.
      results[batch[j]] <- list(reply$value)
    }
    if (progress) message(sprintf("%s: completed %d/%d iterations",
                                  executor$config$progress_label, max(batch), length(X)))
  }
  results
}


# -- Raster compatibility entry points --------------------------------------
default_parallel_workers_raster <- function(reserve_cores = 1L) {
  # Preserve the raster worker-count policy while sharing validation and core
  # detection with the common helper. Vector callers still reserve ten cores.
  default_parallel_workers(reserve_cores = reserve_cores)
}


run_monte_carlo_raster <- function(X,
                                   FUN,
                                   parallel = FALSE,
                                   n_workers = default_parallel_workers_raster(),
                                   random_seed = NULL,
                                   export_names = character(0),
                                   export_env = .GlobalEnv,
                                   progress_label = "Raster Monte Carlo",
                                   progress = TRUE,
                                   chunk_size = NULL,
                                   packages = character()) {
  # Preserve the original eight positional arguments and delegate execution to
  # the shared core. Return every attempted path, including failures and NULLs;
  # raster callers retain their fixed-attempt sample-collection rule.
  run_monte_carlo(
    X = X, FUN = FUN, parallel = parallel, n_workers = n_workers,
    random_seed = random_seed, export_names = export_names,
    export_env = export_env, progress_label = progress_label,
    progress = progress, chunk_size = chunk_size, packages = packages
  )
}


summary_90th_percentile_raster <- function(x) {
  # Compute the 90th percentile of a numeric vector while ignoring missing values.
  as.numeric(stats::quantile(x, probs = 0.9, na.rm = TRUE))
}


call_raster_summary_function <- function(summary_function, zone_values) {
  # Apply a user-supplied summary function to county-level raster values.
  value <- summary_function(zone_values)

  # Downstream sensitivity calculations require exactly one numeric statistic.
  if (!is.numeric(value) || length(value) != 1L) {
    stop("summary_function must return one numeric value.")
  }

  as.numeric(value)
}


build_raster_neighbor_index <- function(r_zones, valid_idx, directions = 4L) {
  # Build a compact neighbor list for valid raster cells.
  # Cell IDs are mapped from full raster cell numbers to dense valid-cell indices.

  if (!is.numeric(directions) || length(directions) != 1L ||
      !is.finite(directions) || !directions %in% c(4L, 8L)) {
    stop("directions must be 4 (rook) or 8 (queen).")
  }

  # ── Extract adjacency using the original raster cell numbers ─────────────
  adj_matrix <- terra::adjacent(
    r_zones,
    cells      = valid_idx,
    directions = directions,
    pairs      = TRUE
  )

  # Keep only pairs whose two cells belong to the zone geometry. Missing data
  # values do not remove zone cells or change their spatial connectivity.
  valid_adj <- adj_matrix[adj_matrix[, 2L] %in% valid_idx, , drop = FALSE]

  # ── Map full raster indices to compact indices 1, ..., n_valid ────────────
  # All Monte Carlo engines work on compact vectors to reduce memory use.
  map_idx <- integer(terra::ncell(r_zones))
  map_idx[valid_idx] <- seq_along(valid_idx)

  adj_pairs <- cbind(
    map_idx[valid_adj[, 1L]],
    map_idx[valid_adj[, 2L]]
  )

  # Store both directions so every cell can retrieve all of its neighbors.
  adj_pairs <- adj_pairs[adj_pairs[, 1L] > 0L & adj_pairs[, 2L] > 0L, , drop = FALSE]
  adj_pairs <- unique(rbind(adj_pairs, adj_pairs[, c(2L, 1L), drop = FALSE]))

  # ── Convert the pair matrix into a per-cell neighbor list ─────────────────
  neighbor_list <- vector("list", length(valid_idx))

  if (nrow(adj_pairs) > 0L) {
    for (i in seq_len(nrow(adj_pairs))) {
      from_cell <- adj_pairs[i, 1L]
      neighbor_list[[from_cell]] <- c(neighbor_list[[from_cell]], adj_pairs[i, 2L])
    }
  }

  list(
    adj_pairs     = adj_pairs,
    neighbor_list = neighbor_list
  )
}


check_zone_connected_after_removal_raster <- function(zone_vals,
                                                      pixel_idx,
                                                      zone_id,
                                                      neighbor_list,
                                                      bfs_env) {
  # Check whether removing one pixel preserves connectivity within its original
  # zone under the selected neighbor rule.

  # Only same-zone neighbors can become disconnected when pixel_idx is removed.
  nbrs_same <- neighbor_list[[pixel_idx]]
  nbrs_same <- nbrs_same[zone_vals[nbrs_same] == zone_id]
  n_target_nbrs <- length(nbrs_same)

  # Zero or one remaining same-zone neighbor is connected by definition.
  if (n_target_nbrs <= 1L) {
    return(TRUE)
  }

  # ── Reuse preallocated BFS workspace ─────────────────────────────────────
  # `visited` and `queue` live in bfs_env so they are not repeatedly allocated
  # for every candidate pixel tested during a Monte Carlo path.
  touched <- integer(1024L)
  n_touched <- 2L
  touched[1L] <- pixel_idx
  touched[2L] <- nbrs_same[1L]

  bfs_env$visited[pixel_idx] <- TRUE
  bfs_env$visited[nbrs_same[1L]] <- TRUE

  found_nbrs <- 1L
  queue_head <- 1L
  queue_tail <- 2L
  bfs_env$queue[queue_head] <- nbrs_same[1L]

  is_connected <- FALSE

  # Search within the original zone while treating pixel_idx as already visited
  # (and therefore unavailable). Connectivity is proved once every same-zone
  # neighbor of the removed pixel is reachable from the first neighbor.
  while (queue_head < queue_tail) {
    curr <- bfs_env$queue[queue_head]
    queue_head <- queue_head + 1L

    new_n <- neighbor_list[[curr]]
    new_n <- new_n[!bfs_env$visited[new_n] & zone_vals[new_n] == zone_id]

    if (length(new_n) > 0L) {
      for (nn in new_n) {
        bfs_env$visited[nn] <- TRUE

        n_touched <- n_touched + 1L
        if (n_touched > length(touched)) {
          length(touched) <- length(touched) * 2L
        }
        touched[n_touched] <- nn

        if (nn %in% nbrs_same) {
          found_nbrs <- found_nbrs + 1L
          if (found_nbrs == n_target_nbrs) {
            is_connected <- TRUE
            break
          }
        }
      }

      if (is_connected) {
        break
      }

      bfs_env$queue[queue_tail:(queue_tail + length(new_n) - 1L)] <- new_n
      queue_tail <- queue_tail + length(new_n)
    }
  }

  # Reset only entries touched by this search before the workspace is reused.
  bfs_env$visited[touched[1L:n_touched]] <- FALSE

  is_connected
}


reassign_one_raster_pixel <- function(zone_vals,
                                      pixel_idx,
                                      zone_counts,
                                      neighbor_list,
                                      bfs_env) {
  # Try to move one boundary pixel from its current zone to an adjacent zone.

  # ── Identify candidate destination zones along the pixel boundary ────────
  own_zone <- zone_vals[pixel_idx]
  own_key <- as.character(own_zone)
  nbrs <- neighbor_list[[pixel_idx]]

  target_cands <- unique(zone_vals[nbrs[zone_vals[nbrs] != own_zone]])

  # An interior pixel has no neighboring zone and therefore cannot be moved.
  if (length(target_cands) == 0L) {
    return(list(result_status = "NoTarget"))
  }

  # Never remove the final pixel of a zone.
  if (zone_counts[own_key] <= 1L) {
    return(list(result_status = "ConnectivityFail"))
  }

  # Reject moves that would fragment the source zone.
  if (!check_zone_connected_after_removal_raster(
    zone_vals      = zone_vals,
    pixel_idx      = pixel_idx,
    zone_id        = own_zone,
    neighbor_list  = neighbor_list,
    bfs_env        = bfs_env
  )) {
    return(list(result_status = "ConnectivityFail"))
  }

  # All remaining destination zones are adjacent; choose one at random to
  # avoid systematically favoring any boundary direction.
  target_zone <- target_cands[sample.int(length(target_cands), 1L)]

  list(
    result_status = "Success",
    target_zone   = target_zone
  )
}


reassign_raster_pixels_by_count <- function(zone_vals,
                                            zone_counts,
                                            k_pixels,
                                            neighbor_list,
                                            adj_pairs,
                                            max_failure = 1000L) {
  # Reassign exactly k boundary pixels while preserving source-zone connectivity.

  # ── Initialize reassignment state ────────────────────────────────────────
  n_cells <- length(zone_vals)
  remaining <- as.integer(k_pixels)
  consecutive_failure <- 0L
  total_failure <- 0L
  n_reassigned <- 0L

  is_candidate <- logical(n_cells)
  is_reassigned <- logical(n_cells)

  # Candidate pixels lie on a boundary between two different zone IDs.
  # A pixel can be moved at most once within a single simulation path.
  diff_mask <- zone_vals[adj_pairs[, 1L]] != zone_vals[adj_pairs[, 2L]]
  boundary_cells <- unique(adj_pairs[diff_mask, 1L])
  is_candidate[boundary_cells] <- TRUE

  # Allocate one BFS workspace for all connectivity checks in this path.
  bfs_env <- new.env(parent = emptyenv())
  bfs_env$visited <- logical(n_cells)
  bfs_env$queue <- integer(n_cells)

  # ── Reassign pixels until the target is reached or failures accumulate ───
  while (consecutive_failure < max_failure && remaining > 0L) {
    avail <- which(is_candidate & !is_reassigned)

    if (length(avail) == 0L) {
      break
    }

    # Remove the selected cell from the current candidate pool before testing;
    # successful boundary changes may add nearby cells back below.
    sel <- if (length(avail) == 1L) avail else sample(avail, 1L)
    is_candidate[sel] <- FALSE

    res <- reassign_one_raster_pixel(
      zone_vals     = zone_vals,
      pixel_idx     = sel,
      zone_counts   = zone_counts,
      neighbor_list = neighbor_list,
      bfs_env       = bfs_env
    )

    if (res$result_status == "Success") {
      # Apply the zone-label change and keep zone sizes synchronized.
      own_key <- as.character(zone_vals[sel])
      target_key <- as.character(res$target_zone)

      zone_vals[sel] <- res$target_zone
      zone_counts[own_key] <- zone_counts[own_key] - 1L
      zone_counts[target_key] <- zone_counts[target_key] + 1L

      remaining <- remaining - 1L
      n_reassigned <- n_reassigned + 1L
      is_reassigned[sel] <- TRUE
      consecutive_failure <- 0L

      # Reinspect neighboring cells because moving sel changes the local zone
      # boundary and can create new eligible reassignment candidates.
      for (nid in neighbor_list[[sel]]) {
        nn <- neighbor_list[[nid]]
        if (length(nn) > 0L && any(zone_vals[nn] != zone_vals[nid])) {
          is_candidate[nid] <- TRUE
        }
      }
    } else {
      # Only consecutive failures control termination; a successful move resets
      # that counter, while total_failure remains available for diagnostics.
      consecutive_failure <- consecutive_failure + 1L
      total_failure <- total_failure + 1L
    }
  }

  # A path is valid only when it completes the requested number of moves.
  status <- if (n_reassigned == k_pixels) "Success" else "Failed"

  list(
    result_status       = status,
    n_reassigned        = n_reassigned,
    total_failure       = total_failure,
    consecutive_failure = consecutive_failure,
    zone_vals           = zone_vals
  )
}


exact_pixel_reassignment_sensitivity_raster <- function(
    r_data,
    r_zones,
    k_pixels,
    summary_function = summary_90th_percentile_raster,
    n_iterations     = 100L,
    random_seed      = NULL,
    max_failure      = 1000L,
    directions       = 4L,
    parallel         = FALSE,
    n_workers        = default_parallel_workers_raster(),
    save_path        = NULL
) {
  # Internal engine for pixel reassignment sensitivity.
  # Each iteration reassigns exactly k_pixels boundary pixels to adjacent zones.
  # The public paper-facing wrapper is pixel_reassignment_sensitivity_raster(),
  # where k_pixels is derived from alpha * s.

  # ── Input validation ─────────────────────────────────────────────────────
  if (!inherits(r_data, "SpatRaster")) {
    stop("r_data must be a terra SpatRaster.")
  }

  if (!inherits(r_zones, "SpatRaster")) {
    stop("r_zones must be a terra SpatRaster.")
  }

  if (terra::nlyr(r_data) != 1L || terra::nlyr(r_zones) != 1L) {
    stop("r_data and r_zones must each have exactly one raster layer.")
  }

  if (!terra::compareGeom(r_data, r_zones, stopOnError = FALSE)) {
    stop("r_data and r_zones must have matching geometry.")
  }

  if (!is.numeric(k_pixels) || length(k_pixels) != 1L || k_pixels < 1L) {
    stop("k_pixels must be one positive integer.")
  }

  if (!is.numeric(n_iterations) || length(n_iterations) != 1L || n_iterations < 1L) {
    stop("n_iterations must be one positive integer.")
  }

  k_pixels <- as.integer(k_pixels)
  n_iterations <- as.integer(n_iterations)

  # ── Extract zone geometry independently of observed data ──────────────────
  message("Preparing raster values and adjacency index...")

  data_values <- terra::values(r_data, mat = FALSE)
  zone_values <- terra::values(r_zones, mat = FALSE)

  valid_mask <- !is.na(zone_values)
  valid_idx <- which(valid_mask)

  if (length(valid_idx) == 0L) {
    stop("No valid zone cells were found in r_zones.")
  }

  zone_vals_init <- zone_values[valid_idx]
  data_vals_init <- data_values[valid_idx]
  zone_counts_init <- table(zone_vals_init)

  # Build cell adjacency once and reuse it in every Monte Carlo iteration.
  neighbor_index <- build_raster_neighbor_index(
    r_zones     = r_zones,
    valid_idx   = valid_idx,
    directions  = directions
  )

  if (nrow(neighbor_index$adj_pairs) == 0L) {
    stop("No adjacent valid raster cells were found.")
  }

  # ── Compute the baseline statistic before any reassignment ───────────────
  zone_means_base <- calculate_raster_zone_values(
    data_values = data_vals_init,
    zone_values = zone_vals_init,
    merge_method = "mean"
  )

  origin_value <- call_raster_summary_function(
    summary_function = summary_function,
    zone_values      = zone_means_base
  )

  message(sprintf("Baseline summary value: %.6f", origin_value))
  message(sprintf("Running raster pixel reassignment engine with k_pixels = %d", k_pixels))

  # Copy values into local bindings captured by the worker closure. This keeps
  # each PSOCK task self-contained when run_one_iteration is serialized.
  adj_pairs_local <- neighbor_index$adj_pairs
  neighbor_list_local <- neighbor_index$neighbor_list
  zone_vals_init_local <- zone_vals_init
  data_vals_init_local <- data_vals_init
  zone_counts_init_local <- zone_counts_init
  summary_function_local <- summary_function
  max_failure_local <- max_failure
  k_pixels_local <- k_pixels

  run_one_iteration <- function(iter) {
    # Each iteration starts from the original zoning, producing independent
    # Monte Carlo paths rather than continuing from the preceding iteration.
    zone_vals <- zone_vals_init_local
    zone_counts <- zone_counts_init_local
    failure_reason <- NA_character_

    res_run <- tryCatch(
      reassign_raster_pixels_by_count(
        zone_vals     = zone_vals,
        zone_counts   = zone_counts,
        k_pixels      = k_pixels_local,
        neighbor_list = neighbor_list_local,
        adj_pairs     = adj_pairs_local,
        max_failure   = max_failure_local
      ),
      error = function(e) {
        failure_reason <<- conditionMessage(e)
        message("Iteration error: ", conditionMessage(e))
        NULL
      }
    )

    # Preserve failed or errored paths in the distribution with NA statistics
    # so success rates and failure diagnostics remain auditable.
    if (is.null(res_run) || res_run$result_status != "Success") {
      return(data.frame(
        iteration           = iter,
        summary_value       = NA_real_,
        n_reassigned        = if (is.null(res_run)) NA_integer_ else res_run$n_reassigned,
        total_failure       = if (is.null(res_run)) NA_integer_ else res_run$total_failure,
        consecutive_failure = if (is.null(res_run)) NA_integer_ else res_run$consecutive_failure,
        status              = if (is.null(res_run)) "Error" else "Failed",
        failure_reason      = if (is.null(res_run)) failure_reason else "Requested pixel count was not reached",
        stringsAsFactors    = FALSE
      ))
    }

    # The data raster never changes; only its grouping by zone label changes.
    # A move can leave a zone without observations. Retain that path as an
    # error diagnostic instead of interrupting all remaining simulations.
    summary_value <- tryCatch({
      zone_means_new <- calculate_raster_zone_values(
        data_values = data_vals_init_local,
        zone_values = res_run$zone_vals,
        merge_method = "mean"
      )
      value <- call_raster_summary_function(
        summary_function = summary_function_local,
        zone_values      = zone_means_new
      )
      if (!is.finite(value)) stop("summary_function returned a non-finite value.")
      value
    }, error = function(e) {
      failure_reason <<- conditionMessage(e)
      message(sprintf("Iteration %d statistic error: %s", iter, failure_reason))
      NA_real_
    })

    data.frame(
      iteration           = iter,
      summary_value       = summary_value,
      n_reassigned        = res_run$n_reassigned,
      total_failure       = res_run$total_failure,
      consecutive_failure = res_run$consecutive_failure,
      status              = if (is.finite(summary_value)) "Success" else "Error",
      failure_reason      = failure_reason,
      stringsAsFactors    = FALSE
    )
  }

  # ── Run independent reassignment paths ───────────────────────────────────
  result_list <- run_monte_carlo_raster(
    X              = seq_len(n_iterations),
    FUN            = run_one_iteration,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = c(
      "check_zone_connected_after_removal_raster",
      "reassign_one_raster_pixel",
      "reassign_raster_pixels_by_count",
      "calculate_raster_zone_values",
      "call_raster_summary_function"
    ),
    export_env     = environment(),
    progress_label = "Raster pixel reassignment"
  )

  distribution <- do.call(rbind, result_list)
  rownames(distribution) <- NULL

  # ── Assemble the result bundle ───────────────────────────────────────────
  out <- list(
    distribution      = distribution,
    origin_value      = origin_value,
    k_pixels          = k_pixels,
    n_iterations      = n_iterations,
    directions        = directions,
    n_valid_cells     = length(valid_idx),
    n_observed_cells  = sum(is.finite(data_vals_init)),
    summary_function  = deparse(substitute(summary_function)),
    sensitivity_type  = "pixel_reassignment_exact_count_raster",
    random_seed       = random_seed,
    parallel          = parallel,
    n_workers         = if (isTRUE(parallel)) n_workers else 1L,
    call              = match.call()
  )

  # ── Persist to disk (optional) ───────────────────────────────────────────
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
    out$save_path <- save_path
  }

  message(sprintf(
    "Raster pixel reassignment engine finished: %d / %d successful iterations",
    sum(distribution$status == "Success", na.rm = TRUE),
    n_iterations
  ))

  out
}


pixel_reassignment_sensitivity_raster <- function(
    r_data,
    r_zones,
    alpha_multiplier,
    summary_function = summary_90th_percentile_raster,
    n_iterations     = 100L,
    random_seed      = NULL,
    max_failure      = 1000L,
    directions       = 4L,
    parallel         = FALSE,
    n_workers        = default_parallel_workers_raster(),
    save_path        = NULL
) {
  # Calculate raster pixel reassignment sensitivity.
  # alpha_multiplier is multiplied by s, where s is the number of pixels
  # covering the average county size.

  # ── Validate the zone raster and estimate the reference county size ──────
  if (!inherits(r_zones, "SpatRaster")) {
    stop("r_zones must be a terra SpatRaster.")
  }

  zone_values <- terra::values(r_zones, mat = FALSE)
  zone_values <- zone_values[!is.na(zone_values)]

  if (length(zone_values) == 0L) {
    stop("No valid zone pixels were found in r_zones.")
  }

  county_pixel_counts <- table(zone_values)
  # alpha_multiplier is the reassigned-area scale relative to the average zone.
  # s_pixels is the mean number of valid raster pixels per zone.
  # k_pixels is the integer number of pixels reassigned in each simulation:
  # k_pixels = round(alpha_multiplier * s_pixels), with a minimum of one pixel.
  s_pixels <- mean(as.numeric(county_pixel_counts), na.rm = TRUE)
  k_pixels <- max(1L, as.integer(round(alpha_multiplier * s_pixels)))

  message(sprintf(
    "Pixel reassignment alpha = %.2f s | s = %.2f pixels | target pixels = %d",
    alpha_multiplier,
    s_pixels,
    k_pixels
  ))

  # ── Delegate to the exact-count Monte Carlo engine ───────────────────────
  out <- exact_pixel_reassignment_sensitivity_raster(
    r_data           = r_data,
    r_zones          = r_zones,
    k_pixels         = k_pixels,
    summary_function = summary_function,
    n_iterations     = n_iterations,
    random_seed      = random_seed,
    max_failure      = max_failure,
    directions       = directions,
    parallel         = parallel,
    n_workers        = n_workers,
    save_path        = save_path
  )

  # Retain the area multiplier and average-zone pixel count alongside k_pixels
  # so the conversion from the area scale to a pixel count remains explicit.
  out$s_pixels <- s_pixels
  out$alpha_multiplier <- alpha_multiplier
  out$sensitivity_type <- "pixel_reassignment_raster"

  if (!is.null(save_path)) {
    # The internal engine saved before wrapper metadata was appended; overwrite
    # the same file once so the on-disk object matches the returned object.
    saveRDS(out, file = save_path)
  }

  out
}


calculate_raster_zone_values <- function(data_values,
                                         zone_values,
                                         area_values = NULL,
                                         merge_method = c("mean", "area_weighted")) {
  # Calculate one value for each raster zone from vectorized raster cells.
  # For PRISM precipitation, area_weighted returns area-weighted mean
  # precipitation for each zone.

  # ── Select the zone aggregation rule ─────────────────────────────────────
  merge_method <- match.arg(merge_method)

  if (merge_method == "area_weighted" && is.null(area_values)) {
    stop("area_values must be supplied when merge_method = 'area_weighted'.")
  }

  # Zone labels define geometry; only finite observations contribute to a
  # statistic. Exclude a missing observation's area from the denominator too.
  observed <- !is.na(zone_values) & is.finite(data_values)
  if (merge_method == "area_weighted") {
    observed <- observed & is.finite(area_values) & area_values > 0
  }

  # Do not silently drop an entire zone or turn its missing values into zero.
  # Report unobserved zone IDs so the input or a proposed split can be checked.
  zone_ids <- sort(unique(zone_values[!is.na(zone_values)]))
  observed_counts <- table(factor(zone_values[observed], levels = zone_ids))
  missing_zones <- zone_ids[observed_counts == 0L]
  if (length(missing_zones) > 0L) {
    stop(sprintf("No valid data values in zone(s): %s.",
                 paste(missing_zones, collapse = ", ")))
  }

  # Equal-cell averaging is appropriate when cells have the same effective
  # area or when an explicitly unweighted statistic has been requested.
  if (merge_method == "mean") {
    zone_means <- tapply(data_values[observed], zone_values[observed], mean)
    return(as.numeric(zone_means))
  }

  # Geographic cells may have unequal areas. Aggregate the weighted numerator
  # and observed-area denominator over exactly the same set of cells.
  weighted_sum <- tapply(data_values[observed] * area_values[observed],
                         zone_values[observed], sum)
  area_sum <- tapply(area_values[observed], zone_values[observed], sum)

  as.numeric(weighted_sum / area_sum)
}


build_raster_zone_adjacency <- function(zone_values,
                                        adj_pairs) {
  # Build a zone-level adjacency list from cell-level adjacency pairs.

  # ── Translate cell boundaries into pairs of neighboring zone IDs ─────────
  zone_pair_matrix <- cbind(
    zone_values[adj_pairs[, 1L]],
    zone_values[adj_pairs[, 2L]]
  )

  # Discard missing cells and within-zone pairs; neither represents a boundary
  # between two distinct zones.
  zone_pair_matrix <- zone_pair_matrix[
    !is.na(zone_pair_matrix[, 1L]) &
      !is.na(zone_pair_matrix[, 2L]) &
      zone_pair_matrix[, 1L] != zone_pair_matrix[, 2L],
    ,
    drop = FALSE
  ]

  # Initialize every zone, including isolated zones with no neighbors.
  zone_ids <- sort(unique(zone_values[!is.na(zone_values)]))
  adjacency_list <- vector("list", length(zone_ids))
  names(adjacency_list) <- as.character(zone_ids)

  if (nrow(zone_pair_matrix) == 0L) {
    return(adjacency_list)
  }

  # Make adjacency symmetric before constructing the named neighbor list.
  zone_pair_matrix <- unique(rbind(
    zone_pair_matrix,
    zone_pair_matrix[, c(2L, 1L), drop = FALSE]
  ))

  for (zone_id in zone_ids) {
    neighbors <- unique(zone_pair_matrix[zone_pair_matrix[, 1L] == zone_id, 2L])
    neighbors <- neighbors[neighbors != zone_id & !is.na(neighbors)]
    adjacency_list[[as.character(zone_id)]] <- neighbors
  }

  adjacency_list
}


build_raster_zone_merge_pairs <- function(adj_list) {
  # Keep each unordered adjacent zone pair once, regardless of zone degree or
  # shared boundary length, so every available pair has equal sampling weight.
  pairs <- lapply(names(adj_list), function(key) {
    zone_id <- as.numeric(key)
    neighbors <- unique(adj_list[[key]])
    neighbors <- neighbors[!is.na(neighbors) & neighbors > zone_id]
    cbind(rep(zone_id, length(neighbors)), neighbors)
  })
  result <- do.call(rbind, pairs)
  if (is.null(result)) result <- matrix(numeric(), ncol = 2L)
  unique(result)
}


update_raster_zone_adjacency_after_merge <- function(adj_list,
                                                     target_id,
                                                     source_id) {
  # Update a zone-level adjacency list after source_id is merged into target_id.

  # ── Combine the two zones' neighborhoods ─────────────────────────────────
  source_key <- as.character(source_id)
  target_key <- as.character(target_id)

  source_neighbors <- adj_list[[source_key]]
  target_neighbors <- adj_list[[target_key]]

  new_neighbors <- unique(c(target_neighbors, source_neighbors))
  new_neighbors <- new_neighbors[new_neighbors != target_id & new_neighbors != source_id]

  # The absorbed source zone disappears, while the target inherits every
  # external neighbor of either original zone.
  adj_list[[target_key]] <- new_neighbors
  adj_list[[source_key]] <- NULL

  # ── Redirect all references from source_id to target_id ──────────────────
  for (zone_key in names(adj_list)) {
    neighbors <- adj_list[[zone_key]]

    if (source_id %in% neighbors) {
      neighbors <- neighbors[neighbors != source_id]

      if (!(target_id %in% neighbors) && zone_key != target_key) {
        neighbors <- c(neighbors, target_id)
      }

      adj_list[[zone_key]] <- unique(neighbors)
    }
  }

  adj_list
}


merging_sensitivity_raster <- function(
    r_data,
    r_zones,
    summary_function = summary_90th_percentile_raster,
    k_order          = 3L,
    n_iterations     = 100L,
    merge_method     = c("mean", "area_weighted"),
    random_seed      = NULL,
    directions       = 4L,
    parallel         = FALSE,
    n_workers        = default_parallel_workers_raster(),
    save_path        = NULL
) {
  # Calculate raster-based first- through kth-order merging sensitivity.
  # Each Monte Carlo simulation starts from the original spatial partitioning
  # scheme and performs up to k consecutive merges. After each merge, the summary
  # statistic is recalculated and recorded.

  # ── Input validation ─────────────────────────────────────────────────────
  merge_method <- match.arg(merge_method)

  if (!inherits(r_data, "SpatRaster")) {
    stop("r_data must be a terra SpatRaster.")
  }

  if (!inherits(r_zones, "SpatRaster")) {
    stop("r_zones must be a terra SpatRaster.")
  }

  if (terra::nlyr(r_data) != 1L || terra::nlyr(r_zones) != 1L) {
    stop("r_data and r_zones must each have exactly one raster layer.")
  }

  if (!terra::compareGeom(r_data, r_zones, stopOnError = FALSE)) {
    stop("r_data and r_zones must have matching geometry.")
  }

  if (!is.numeric(k_order) || length(k_order) != 1L || k_order < 1L) {
    stop("k_order must be one positive integer.")
  }

  if (!is.numeric(n_iterations) || length(n_iterations) != 1L || n_iterations < 1L) {
    stop("n_iterations must be one positive integer.")
  }

  k_order <- as.integer(k_order)
  n_iterations <- as.integer(n_iterations)

  # ── Extract zone geometry and precompute topology ─────────────────────────
  message("Preparing raster merging sensitivity inputs...")

  data_values_all <- terra::values(r_data, mat = FALSE)
  zone_values_all <- terra::values(r_zones, mat = FALSE)

  valid_mask <- !is.na(zone_values_all)
  valid_idx <- which(valid_mask)

  if (length(valid_idx) == 0L) {
    stop("No valid zone cells were found in r_zones.")
  }

  data_values <- data_values_all[valid_idx]
  zone_values_init <- zone_values_all[valid_idx]

  area_values <- NULL
  if (merge_method == "area_weighted") {
    # Cell areas are computed only when requested because cellSize can be
    # expensive for large rasters.
    area_raster <- terra::cellSize(r_data)
    area_values <- terra::values(area_raster, mat = FALSE)[valid_idx]
  }

  neighbor_index <- build_raster_neighbor_index(
    r_zones    = r_zones,
    valid_idx  = valid_idx,
    directions = directions
  )

  if (nrow(neighbor_index$adj_pairs) == 0L) {
    stop("No adjacent valid raster cells were found.")
  }

  # Cell adjacency is converted to a zone graph once, then updated cheaply
  # after each simulated merge rather than recalculated from the raster.
  initial_adj_list <- build_raster_zone_adjacency(
    zone_values = zone_values_init,
    adj_pairs   = neighbor_index$adj_pairs
  )

  # ── Compute the baseline statistic ────────────────────────────────────────
  baseline_values <- calculate_raster_zone_values(
    data_values  = data_values,
    zone_values  = zone_values_init,
    area_values  = area_values,
    merge_method = merge_method
  )

  origin_value <- call_raster_summary_function(
    summary_function = summary_function,
    zone_values      = baseline_values
  )

  message(sprintf("Baseline summary value: %.6f", origin_value))
  message(sprintf("Running raster merging sensitivity for k = 1 to %d", k_order))

  # Bind all read-only inputs locally for safe serialization to PSOCK workers.
  data_values_local <- data_values
  zone_values_init_local <- zone_values_init
  area_values_local <- area_values
  initial_adj_list_local <- initial_adj_list
  summary_function_local <- summary_function
  merge_method_local <- merge_method
  k_order_local <- k_order

  run_one_merge_path <- function(iter) {
    # Every path starts from the unchanged zone labels and adjacency graph.
    current_zone_values <- zone_values_init_local
    current_adj_list <- initial_adj_list_local

    rows <- vector("list", k_order_local)

    # ── Sequentially extend one path from first- to kth-order merging ───────
    for (merge_level in seq_len(k_order_local)) {
      # Rebuild the unique pair table after each graph contraction. Sampling
      # zones first would give pairs involving low-degree zones greater weight.
      merge_pairs <- build_raster_zone_merge_pairs(current_adj_list)

      # If the graph has no remaining edge, higher-order merges are impossible;
      # retain a failed row so every iteration has the same k-level schema.
      if (nrow(merge_pairs) == 0L) {
        rows[[merge_level]] <- data.frame(
          iteration        = iter,
          k                = merge_level,
          target_id        = NA_real_,
          source_id        = NA_real_,
          summary_value    = NA_real_,
          status           = "Failed",
          stringsAsFactors = FALSE
        )
        next
      }

      # Sample one unordered pair uniformly, then randomize which ID is kept.
      # Index-based sampling also handles singleton numeric ID vectors safely.
      selected_pair <- merge_pairs[sample.int(nrow(merge_pairs), 1L), ]
      retained_index <- sample.int(2L, 1L)
      target_id <- selected_pair[retained_index]
      source_id <- selected_pair[3L - retained_index]

      current_zone_values[current_zone_values == source_id] <- target_id

      current_adj_list <- update_raster_zone_adjacency_after_merge(
        adj_list   = current_adj_list,
        target_id  = target_id,
        source_id  = source_id
      )

      # Reaggregate the unchanged data values under the updated zone labels.
      current_values <- calculate_raster_zone_values(
        data_values  = data_values_local,
        zone_values  = current_zone_values,
        area_values  = area_values_local,
        merge_method = merge_method_local
      )

      summary_value <- call_raster_summary_function(
        summary_function = summary_function_local,
        zone_values      = current_values
      )

      rows[[merge_level]] <- data.frame(
        iteration        = iter,
        k                = merge_level,
        target_id        = target_id,
        source_id        = source_id,
        summary_value    = summary_value,
        status           = "Success",
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, rows)
  }

  # ── Run independent merge paths ──────────────────────────────────────────
  result_list <- run_monte_carlo_raster(
    X              = seq_len(n_iterations),
    FUN            = run_one_merge_path,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = c(
      "calculate_raster_zone_values",
      "build_raster_zone_merge_pairs",
      "update_raster_zone_adjacency_after_merge",
      "call_raster_summary_function"
    ),
    export_env     = environment(),
    progress_label = "Raster merging sensitivity"
  )

  distribution <- do.call(rbind, result_list)
  rownames(distribution) <- NULL

  # ── Assemble the result bundle ───────────────────────────────────────────
  out <- list(
    distribution      = distribution,
    origin_value      = origin_value,
    merge_method      = merge_method,
    merge_sampling    = "uniform_adjacent_pairs",
    k_order           = k_order,
    n_iterations      = n_iterations,
    directions        = directions,
    n_valid_cells     = length(valid_idx),
    n_observed_cells  = sum(is.finite(data_values)),
    summary_function  = deparse(substitute(summary_function)),
    sensitivity_type  = "merging_raster",
    random_seed       = random_seed,
    parallel          = parallel,
    n_workers         = if (isTRUE(parallel)) n_workers else 1L,
    call              = match.call()
  )

  # ── Persist to disk (optional) ───────────────────────────────────────────
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
    out$save_path <- save_path
  }

  message(sprintf(
    "Raster merging sensitivity finished: %d successful records",
    sum(distribution$status == "Success", na.rm = TRUE)
  ))

  out
}


split_raster_zone_values_once <- function(zone_values,
                                          target_id,
                                          neighbor_list,
                                          next_zone_id,
                                          n_attempts = 10L) {
  # Randomly split one raster zone using a two-seed region-growing algorithm.
  # Cells grow from two randomly selected seed cells to form two connected regions.
  # One region keeps the original zone ID, and the other receives a new zone ID.

  # ── Extract the target zone and build its internal adjacency graph ────────
  zone_cells <- which(zone_values == target_id)

  # Both child regions must contain at least one cell.
  if (length(zone_cells) < 2L) {
    return(NULL)
  }

  n_nodes <- length(zone_cells)
  cell_to_node <- integer(length(zone_values))
  cell_to_node[zone_cells] <- seq_len(n_nodes)

  adj_list <- vector("list", n_nodes)

  for (node in seq_len(n_nodes)) {
    global_cell <- zone_cells[node]
    nbrs <- neighbor_list[[global_cell]]
    nbrs <- nbrs[zone_values[nbrs] == target_id]
    adj_list[[node]] <- unique(cell_to_node[nbrs])
    adj_list[[node]] <- adj_list[[node]][adj_list[[node]] > 0L]
  }

  # ── Grow two connected regions from randomly selected starting cells ─────
  for (attempt in seq_len(as.integer(n_attempts))) {
    # Each seed anchors one child region. Alternating queue expansion gives
    # both regions an opportunity to grow through the target-zone graph.
    seeds <- sample(seq_len(n_nodes), 2L, replace = FALSE)

    assigned <- integer(n_nodes)
    assigned[seeds[1L]] <- 1L
    assigned[seeds[2L]] <- 2L

    queue1 <- seeds[1L]
    queue2 <- seeds[2L]

    while (length(queue1) > 0L || length(queue2) > 0L) {
      if (length(queue1) > 0L) {
        curr <- queue1[1L]
        queue1 <- queue1[-1L]
        nbrs <- adj_list[[curr]]
        new_n <- nbrs[assigned[nbrs] == 0L]

        if (length(new_n) > 0L) {
          assigned[new_n] <- 1L
          queue1 <- c(queue1, new_n)
        }
      }

      if (length(queue2) > 0L) {
        curr <- queue2[1L]
        queue2 <- queue2[-1L]
        nbrs <- adj_list[[curr]]
        new_n <- nbrs[assigned[nbrs] == 0L]

        if (length(new_n) > 0L) {
          assigned[new_n] <- 2L
          queue2 <- c(queue2, new_n)
        }
      }
    }

    # Unreached nodes belong to components without a seed. Retry instead of
    # assigning disconnected cells arbitrarily to one of the child regions.
    if (any(assigned == 0L)) next

    region1_cells <- zone_cells[assigned == 1L]
    region2_cells <- zone_cells[assigned == 2L]

    # A valid split must leave both child regions non-empty.
    if (length(region1_cells) == 0L || length(region2_cells) == 0L) {
      next
    }

    # Preserve target_id for the first child and assign a fresh ID to the
    # second child.
    new_zone_values <- zone_values
    new_zone_values[region2_cells] <- next_zone_id

    return(list(
      zone_values   = new_zone_values,
      split_zone_id = target_id,
      new_zone_id   = next_zone_id,
      region1_cells = region1_cells,
      region2_cells = region2_cells
    ))
  }

  NULL
}


splitting_sensitivity_raster <- function(
    r_data,
    r_zones,
    summary_function = summary_90th_percentile_raster,
    k_order          = 3L,
    n_iterations     = 100L,
    split_method     = c("mean", "area_weighted"),
    random_seed      = NULL,
    n_attempts       = 10L,
    directions       = 4L,
    parallel         = FALSE,
    n_workers        = default_parallel_workers_raster(),
    save_path        = NULL,
    max_split_attempts = 100L
) {
  # Calculate raster-based first- through kth-order splitting sensitivity.
  # Each Monte Carlo simulation starts from the original zoning and performs
  # up to k consecutive splits. After each split, the summary statistic is
  # recalculated and recorded.

  # ── Input validation ─────────────────────────────────────────────────────
  split_method <- match.arg(split_method)

  if (!is.numeric(max_split_attempts) || length(max_split_attempts) != 1L ||
      !is.finite(max_split_attempts) || max_split_attempts < 1L ||
      max_split_attempts != floor(max_split_attempts) ||
      max_split_attempts > .Machine$integer.max) {
    stop("max_split_attempts must be one positive integer.")
  }
  if (!is.numeric(n_attempts) || length(n_attempts) != 1L ||
      !is.finite(n_attempts) || n_attempts < 1L ||
      n_attempts != floor(n_attempts) || n_attempts > .Machine$integer.max) {
    stop("n_attempts must be one positive integer.")
  }

  if (!inherits(r_data, "SpatRaster")) {
    stop("r_data must be a terra SpatRaster.")
  }

  if (!inherits(r_zones, "SpatRaster")) {
    stop("r_zones must be a terra SpatRaster.")
  }

  if (terra::nlyr(r_data) != 1L || terra::nlyr(r_zones) != 1L) {
    stop("r_data and r_zones must each have exactly one raster layer.")
  }

  if (!terra::compareGeom(r_data, r_zones, stopOnError = FALSE)) {
    stop("r_data and r_zones must have matching geometry.")
  }

  if (!is.numeric(k_order) || length(k_order) != 1L || k_order < 1L) {
    stop("k_order must be one positive integer.")
  }

  if (!is.numeric(n_iterations) || length(n_iterations) != 1L || n_iterations < 1L) {
    stop("n_iterations must be one positive integer.")
  }

  k_order <- as.integer(k_order)
  n_iterations <- as.integer(n_iterations)

  # ── Extract zone geometry and precompute cell adjacency ────────────────────
  message("Preparing raster splitting sensitivity inputs...")

  data_values_all <- terra::values(r_data, mat = FALSE)
  zone_values_all <- terra::values(r_zones, mat = FALSE)

  valid_mask <- !is.na(zone_values_all)
  valid_idx <- which(valid_mask)

  if (length(valid_idx) == 0L) {
    stop("No valid zone cells were found in r_zones.")
  }

  data_values <- data_values_all[valid_idx]
  zone_values_init <- zone_values_all[valid_idx]

  area_values <- NULL
  if (split_method == "area_weighted") {
    # Compute cell areas only for the area-weighted aggregation option.
    area_raster <- terra::cellSize(r_data)
    area_values <- terra::values(area_raster, mat = FALSE)[valid_idx]
  }

  neighbor_index <- build_raster_neighbor_index(
    r_zones    = r_zones,
    valid_idx  = valid_idx,
    directions = directions
  )

  if (length(neighbor_index$neighbor_list) == 0L) {
    stop("No valid raster-cell neighbors were found.")
  }

  # ── Compute the baseline statistic and initial zone state ─────────────────
  baseline_values <- calculate_raster_zone_values(
    data_values  = data_values,
    zone_values  = zone_values_init,
    area_values  = area_values,
    merge_method = split_method
  )

  origin_value <- call_raster_summary_function(
    summary_function = summary_function,
    zone_values      = baseline_values
  )

  # Zone sizes determine split eligibility; new IDs increase monotonically so
  # they cannot collide with any original zone identifier.
  initial_cell_counts <- table(zone_values_init)
  next_zone_id_init <- max(zone_values_init, na.rm = TRUE) + 1L

  message(sprintf("Baseline summary value: %.6f", origin_value))
  message(sprintf("Running raster splitting sensitivity for k = 1 to %d", k_order))

  # Bind all read-only inputs locally for safe serialization to PSOCK workers.
  data_values_local <- data_values
  zone_values_init_local <- zone_values_init
  area_values_local <- area_values
  neighbor_list_local <- neighbor_index$neighbor_list
  summary_function_local <- summary_function
  split_method_local <- split_method
  k_order_local <- k_order
  n_attempts_local <- as.integer(n_attempts)
  max_split_attempts_local <- as.integer(max_split_attempts)
  initial_cell_counts_local <- initial_cell_counts
  next_zone_id_init_local <- next_zone_id_init

  run_one_split_path <- function(iter) {
    # Every iteration begins with the original zoning and an unused new-zone ID.
    current_zone_values <- zone_values_init_local
    cell_counts <- initial_cell_counts_local
    next_zone_id <- next_zone_id_init_local

    # Failed rows retain the requested k-level schema, but never contribute a
    # statistic. Only completed splits are counted in successful_splits.
    rows <- lapply(seq_len(k_order_local), function(level) data.frame(
      iteration = iter, k = level, split_id = NA_real_, new_id = NA_real_,
      summary_value = NA_real_, status = "Failed", step_attempts = 0L,
      successful_splits = 0L, failure_reason = "Previous split order was not completed",
      stringsAsFactors = FALSE
    ))

    # ── Sequentially extend one path from first- to kth-order splitting ─────
    for (split_level in seq_len(k_order_local)) {
      # A zone needs at least two cells to form two non-empty child regions.
      eligible_ids <- as.numeric(names(cell_counts[cell_counts >= 2L]))

      target_id <- NA_real_
      split_result <- NULL
      step_attempts <- 0L
      failure_reason <- "No zone contains enough cells for another split"

      # Retry the same order from the unchanged state. n_attempts bounds seed
      # trials per target; max_split_attempts bounds target trials per order.
      if (length(eligible_ids) > 0L) {
        for (split_attempt in seq_len(max_split_attempts_local)) {
          step_attempts <- split_attempt
          target_id <- eligible_ids[sample.int(length(eligible_ids), 1L)]
          split_result <- split_raster_zone_values_once(
            zone_values   = current_zone_values,
            target_id     = target_id,
            neighbor_list = neighbor_list_local,
            next_zone_id  = next_zone_id,
            n_attempts    = n_attempts_local
          )
          if (is.null(split_result)) {
            failure_reason <- "No connected split within the retry limit"
            next
          }

          # Missing-data cells remain part of both child geometries, but each
          # child must contain observations to support a zone statistic.
          child1_observed <- any(is.finite(data_values_local[split_result$region1_cells]))
          child2_observed <- any(is.finite(data_values_local[split_result$region2_cells]))
          if (!child1_observed || !child2_observed) {
            failure_reason <- "A split child has no valid data values"
            split_result <- NULL
            next
          }
          break
        }
      }

      # A failed order ends this path; later rows cannot be labeled successful
      # without the required preceding splits. Keep the failure for diagnostics.
      if (is.null(split_result)) {
        rows[[split_level]]$split_id <- target_id
        rows[[split_level]]$step_attempts <- step_attempts
        rows[[split_level]]$failure_reason <- failure_reason
        for (level in seq.int(split_level, k_order_local)) {
          rows[[level]]$successful_splits <- split_level - 1L
        }
        break
      }

      # Update zone labels, child sizes, and the next unused ID in lockstep.
      current_zone_values <- split_result$zone_values

      cell_counts[as.character(target_id)] <- length(split_result$region1_cells)
      cell_counts[as.character(split_result$new_zone_id)] <- length(split_result$region2_cells)
      next_zone_id <- next_zone_id + 1L

      # Reaggregate the unchanged raster data under the new child-zone labels.
      current_values <- calculate_raster_zone_values(
        data_values  = data_values_local,
        zone_values  = current_zone_values,
        area_values  = area_values_local,
        merge_method = split_method_local
      )

      summary_value <- call_raster_summary_function(
        summary_function = summary_function_local,
        zone_values      = current_values
      )

      rows[[split_level]] <- data.frame(
        iteration        = iter,
        k                = split_level,
        split_id         = target_id,
        new_id           = split_result$new_zone_id,
        summary_value    = summary_value,
        status           = "Success",
        step_attempts    = step_attempts,
        successful_splits = split_level,
        failure_reason   = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, rows)
  }

  # ── Run independent split paths ──────────────────────────────────────────
  result_list <- run_monte_carlo_raster(
    X              = seq_len(n_iterations),
    FUN            = run_one_split_path,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = c(
      "split_raster_zone_values_once",
      "calculate_raster_zone_values",
      "call_raster_summary_function"
    ),
    export_env     = environment(),
    progress_label = "Raster splitting sensitivity"
  )

  distribution <- do.call(rbind, result_list)
  rownames(distribution) <- NULL

  # ── Assemble the result bundle ───────────────────────────────────────────
  out <- list(
    distribution      = distribution,
    origin_value      = origin_value,
    split_method      = split_method,
    k_order           = k_order,
    n_iterations      = n_iterations,
    n_attempts        = as.integer(n_attempts),
    max_split_attempts = as.integer(max_split_attempts),
    directions        = directions,
    n_valid_cells     = length(valid_idx),
    n_observed_cells  = sum(is.finite(data_values)),
    summary_function  = deparse(substitute(summary_function)),
    sensitivity_type  = "splitting_raster",
    random_seed       = random_seed,
    parallel          = parallel,
    n_workers         = if (isTRUE(parallel)) n_workers else 1L,
    call              = match.call()
  )

  # ── Persist to disk (optional) ───────────────────────────────────────────
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
    out$save_path <- save_path
  }

  message(sprintf(
    "Raster splitting sensitivity finished: %d successful records",
    sum(distribution$status == "Success", na.rm = TRUE)
  ))

  out
}


