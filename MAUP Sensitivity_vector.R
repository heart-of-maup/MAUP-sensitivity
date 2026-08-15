# On the sensitivities to the modifiable areal unit problem
# Vector sensitivity calculation functions
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# yexiang@nnu.edu.cn
# 2026-08-16

# Please cite the following reference when part or all of the code in this file
# is reused under the license of CC-BY-4.0:
# Ye, X., & Chen, J. (2026). On the sensitivities to the modifiable areal unit
# problem. Big Earth Data, 1–36. https://doi.org/10.1080/20964471.2026.2692263

# Required packages are checked when the module is sourced. Data import and
# plotting packages belong to the walkthrough and plotting module, respectively.
required_vector_packages <- c("sf", "spdep")
missing_vector_packages <- required_vector_packages[
  !vapply(required_vector_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_vector_packages) > 0L) {
  stop(sprintf(
    "Install the required vector package(s): %s",
    paste(missing_vector_packages, collapse = ", ")
  ))
}

# Attach sf and spdep because the calculation code uses their spatial verbs
# repeatedly. The optional splitting dependency lwgeom is checked only when a
# polygon split is requested, so merging and reassignment remain usable alone.
suppressPackageStartupMessages({
  library(sf)
  library(spdep)
})

validate_summary_fields <- function(sf, summary_field, object_name = "sf") {
  # Validate that the requested summary field(s) exist and are numeric.
  # The returned field vector is safe to pass into summary functions.
  
  if (missing(summary_field) || is.null(summary_field)) {
    stop("summary_field must be provided")
  }
  
  # Coerce to character so both bare names and character vectors are handled
  # consistently by the checks below.
  summary_field <- as.character(summary_field)
  
  # Reject empty, missing, or blank field names before looking them up.
  if (length(summary_field) < 1L || any(is.na(summary_field)) || any(!nzchar(summary_field))) {
    stop("summary_field must contain at least one valid field name")
  }
  
  # Check that every requested field is present in the supplied data object.
  missing_fields <- setdiff(summary_field, names(sf))
  if (length(missing_fields) > 0L) {
    stop(sprintf(
      "Field(s) %s do not exist in %s",
      paste(sprintf("'%s'", missing_fields), collapse = ", "),
      object_name
    ))
  }
  
  # Sensitivity statistics are computed from numeric vectors only.
  non_numeric_fields <- summary_field[!vapply(
    summary_field,
    function(field_name) is.numeric(sf[[field_name]]),
    logical(1)
  )]
  if (length(non_numeric_fields) > 0L) {
    stop(sprintf(
      "Field(s) %s must be numeric",
      paste(sprintf("'%s'", non_numeric_fields), collapse = ", ")
    ))
  }
  
  summary_field
}

call_summary_function <- function(sf, summary_field, summary_function) {
  # Extract one or more fields from an sf object and pass them to the selected
  # summary function as ordinary vectors.
  
  args <- lapply(summary_field, function(field_name) sf[[field_name]])
  do.call(summary_function, unname(args))
}

default_parallel_workers <- function(reserve_cores = 1L) {
  # Choose a conservative number of parallel workers while leaving cores free
  # for the operating system and other tasks.
  
  n_cores <- parallel::detectCores(logical = TRUE)
  
  if (is.na(n_cores) || n_cores < 2L) {
    return(1L)
  }
  
  max(1L, n_cores - reserve_cores)
}

run_monte_carlo <- function(X,
                            FUN,
                            parallel  = FALSE,
                            n_workers = default_parallel_workers(),
                            random_seed = NULL,
                            export_names = character(0),
                            export_env = .GlobalEnv,
                            progress_label = "Monte Carlo") {
  # Run a Monte-Carlo job either serially or with a temporary PSOCK cluster.
  # This wrapper keeps random seeding and worker setup consistent.

  if (!is.function(FUN)) {
    stop("FUN must be a function")
  }

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE")
  }

  if (!is.null(random_seed)) {
    if (!is.numeric(random_seed) || length(random_seed) != 1L ||
        !is.finite(random_seed) || random_seed < 0 ||
        random_seed != floor(random_seed) ||
        random_seed > .Machine$integer.max) {
      stop("random_seed must be NULL or one non-negative integer")
    }
    random_seed <- as.integer(random_seed)
  }

  if (length(X) == 0L) {
    return(list())
  }

  run_serial <- function() {
    if (!is.null(random_seed)) {
      set.seed(random_seed)
    }
    lapply(X, FUN)
  }

  if (!isTRUE(parallel)) {
    return(run_serial())
  }

  if (!is.numeric(n_workers) || length(n_workers) != 1L ||
      !is.finite(n_workers) || n_workers < 1L ||
      n_workers != floor(n_workers)) {
    stop("n_workers must be one positive integer")
  }

  n_workers <- as.integer(n_workers)
  n_workers <- max(1L, min(n_workers, length(X)))

  if (n_workers <= 1L || length(X) <= 1L) {
    message(sprintf(
      "%s: one worker/task available; using serial execution",
      progress_label
    ))
    return(run_serial())
  }
  
  message(sprintf(
    "%s: using %d parallel workers for %d iterations",
    progress_label, n_workers, length(X)
  ))
  
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  if (!is.null(random_seed)) {
    parallel::clusterSetRNGStream(cl, random_seed)
  }
  
  parallel::clusterEvalQ(cl, {
    for (pkg in c("sf", "lwgeom", "spdep")) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        suppressPackageStartupMessages(
          library(pkg, character.only = TRUE)
        )
      }
    }
    NULL
  })
  
  export_names <- unique(as.character(export_names))
  export_names <- export_names[nzchar(export_names)]
  if (length(export_names) > 0L) {
    parallel::clusterExport(cl, export_names, envir = export_env)
  }
  
  # Run parallel jobs in batches so the main process can report progress
  # instead of staying silent until every worker finishes.
  batch_size <- max(n_workers, ceiling(length(X) / 10))
  batch_ids <- split(seq_along(X), ceiling(seq_along(X) / batch_size))
  results <- vector("list", length(X))
  
  for (batch in batch_ids) {
    batch_results <- parallel::parLapply(cl, X[batch], FUN)
    results[batch] <- batch_results
    
    n_done <- max(batch)
    message(sprintf(
      "%s: completed %d/%d iterations",
      progress_label, n_done, length(X)
    ))
  }
  
  results
}

summary_function_export_names <- function(summary_function) {
  # Identify external variables used inside a user-supplied summary function
  # so PSOCK workers can see objects such as spatial weights.
  
  if (!is.function(summary_function) ||
      !requireNamespace("codetools", quietly = TRUE)) {
    return(character(0))
  }
  
  globals <- codetools::findGlobals(summary_function, merge = FALSE)
  variables <- globals$variables
  if (is.null(variables)) {
    return(character(0))
  }
  
  variables <- unique(as.character(variables))
  variables[nzchar(variables)]
}

merging_operation <- function(sf, i, j, field_rules) {
  # Merge two areal units and recalculate their attributes according to
  # the user-supplied field rules.
  
  # ── Locate target row indices ─────────────────────────────────────────────
  idx_i <- which(sf$ID == i)
  idx_j <- which(sf$ID == j)
  
  if (length(idx_i) == 0) stop(sprintf("ID '%s' does not exist in the sf object", i))
  if (length(idx_j) == 0) stop(sprintf("ID '%s' does not exist in the sf object", j))
  
  # ── Compute areas of the two units (used for the area_weighted rule) ──────
  area_i     <- as.numeric(st_area(sf[idx_i, ]))
  area_j     <- as.numeric(st_area(sf[idx_j, ]))
  total_area <- area_i + area_j
  
  # ── Merge geometries; new ID is current max + 1 ───────────────────────────
  merged_geom <- st_union(sf$geometry[c(idx_i, idx_j)])
  new_id      <- max(sf$ID, na.rm = TRUE) + 1L
  
  # Use row i as the template for the new row, then replace geometry and ID
  new_row    <- sf[idx_i, ]
  new_row$ID <- new_id
  st_geometry(new_row) <- merged_geom
  
  if (!is.null(field_rules)) {
    
    # ── Step 1: handle all basic rules (sum / area_weighted / reserve) ──────
    for (field_name in names(field_rules)) {
      rule <- field_rules[[field_name]]$merge
      if (is.list(rule)) next
      if (!field_name %in% names(sf)) next
      
      val_i <- sf[[field_name]][idx_i]
      val_j <- sf[[field_name]][idx_j]
      new_row[[field_name]] <- switch(rule,
                                      sum           = val_i + val_j,
                                      area_weighted = (val_i * area_i + val_j * area_j) / total_area,
                                      reserve       = val_i,
                                      stop(sprintf("Unknown merge rule: '%s'", rule))
      )
    }
    
    # ── Step 2: recompute derived fields using the updated count fields ─────
    for (field_name in names(field_rules)) {
      rule <- field_rules[[field_name]]$merge
      if (!is.list(rule)) next
      
      num_field   <- rule$numerator
      denom_field <- rule$denominator
      
      if (!num_field %in% names(new_row))   stop(sprintf("Numerator field '%s' does not exist", num_field))
      if (!denom_field %in% names(new_row)) stop(sprintf("Denominator field '%s' does not exist", denom_field))
      
      denom_val <- new_row[[denom_field]]
      new_row[[field_name]] <- if (!is.na(denom_val) && denom_val != 0) {
        new_row[[num_field]] / denom_val
      } else {
        NA_real_
      }
    }
  }
  
  # ── Remove the two original rows and append the merged new row ────────────
  rbind(sf[-c(idx_i, idx_j), ], new_row)
}

calculate_merge_complexity <- function(sf_object, n_simulations = 1000) {
  # Estimate the number of possible k-step merge paths using Monte-Carlo
  # simulations on the rook-contiguity graph of the input sf object.
  
  # ---- Build the initial topology graph -------------------------------
  message("Building topology graph...")
  nb <- spdep::poly2nb(sf_object, queen = FALSE)
  
  if (any(spdep::card(nb) == 0)) {
    warning("Some isolated regions cannot be merged.")
  }
  
  adj_mat <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  g_initial <- igraph::graph_from_adjacency_matrix(
    adj_mat,
    mode = "undirected"
  )
  g_initial <- igraph::simplify(
    g_initial,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )

  max_k <- igraph::vcount(g_initial) - 1L
  if (max_k < 1L) {
    return(data.frame(
      k_order = integer(),
      avg_available_merges = numeric(),
      estimated_total_combinations_log10 = numeric(),
      estimated_total_combinations = numeric()
    ))
  }
  
  available_choices_matrix <- matrix(
    NA_real_,
    nrow = n_simulations,
    ncol = max_k
  )
  
  message(sprintf(
    "Starting %d topology merge simulations...",
    n_simulations
  ))
  
  # ---- Monte-Carlo topology merge simulations -------------------------
  for (sim in seq_len(n_simulations)) {
    g_current <- g_initial
    
    for (k in seq_len(max_k)) {
      n_edges <- igraph::ecount(g_current)
      available_choices_matrix[sim, k] <- n_edges
      
      if (n_edges == 0L) break
      
      # Randomly choose one adjacent pair to merge.
      edge_idx  <- sample.int(n_edges, 1L)
      edge_ends <- igraph::ends(g_current, edge_idx, names = FALSE)
      v1 <- edge_ends[1]
      v2 <- edge_ends[2]
      
      # Contract v2 into v1 by reconnecting v2's neighbors to v1.
      neighbors_v2 <- igraph::neighbors(g_current, v2)
      neighbors_v2 <- neighbors_v2[neighbors_v2 != v1]
      
      if (length(neighbors_v2) > 0L) {
        new_edges <- as.vector(
          rbind(rep(v1, length(neighbors_v2)), neighbors_v2)
        )
        g_current <- igraph::add_edges(g_current, new_edges)
      }
      
      # Remove v2 and simplify the topology graph.
      g_current <- igraph::delete_vertices(g_current, v2)
      g_current <- igraph::simplify(
        g_current,
        remove.multiple = TRUE,
        remove.loops = TRUE
      )
    }
    
    if (sim %% 100L == 0L) {
      message(sprintf(
        "Completed simulations: %d/%d",
        sim,
        n_simulations
      ))
    }
  }
  
  # ---- Summarize simulated branching factors --------------------------
  message("Calculating complexity summary...")
  
  avg_choices <- colMeans(available_choices_matrix, na.rm = TRUE)
  avg_choices[!is.finite(avg_choices)] <- 0
  
  log10_choices <- ifelse(avg_choices > 0, log10(avg_choices), -Inf)
  log10_cumulative <- cumsum(log10_choices)
  
  data.frame(
    k_order = seq_len(max_k),
    avg_available_merges = round(avg_choices, 2),
    estimated_total_combinations_log10 = round(log10_cumulative, 2),
    estimated_total_combinations = 10^log10_cumulative
  )
}

merging_sensitivity <- function(sf,
                                field_rules,
                                summary_field,
                                summary_function,
                                k_order              = 1,
                                exhaustive_threshold = 1e3,
                                n_iterations         = 100,
                                random_seed          = NULL,
                                n_simulations        = 200,
                                parallel             = FALSE,
                                n_workers            = default_parallel_workers(),
                                save_path            = NULL) {
  # Estimate how much a summary statistic changes after k sequential merges.
  # Small candidate spaces are enumerated; larger ones are sampled randomly.
  
  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.function(summary_function))         stop("summary_function must be a function")
  summary_field <- validate_summary_fields(sf, summary_field)
  if (k_order < 1)                            stop("k_order must be >= 1")
  
  if (!is.null(random_seed)) set.seed(random_seed)
  
  # ── Compute baseline statistic ────────────────────────────────────────────
  origin_value <- tryCatch({
    call_summary_function(sf, summary_field, summary_function)
  }, error = function(e) {
    stop(paste0("summary_function failed on the initial data: ", e$message))
  })
  
  # ── Estimate the number of combinations to decide exhaustive vs sampling ──
  message("Estimating the complexity of k-order merging...")
  complexity_results     <- calculate_merge_complexity(sf, n_simulations = n_simulations)
  target_complexity      <- complexity_results[complexity_results$k_order == k_order, ]
  estimated_combinations <- target_complexity$estimated_total_combinations
  
  message(sprintf("Estimated number of combinations for %d-order merging: %.2e", k_order, estimated_combinations))
  
  use_exhaustive  <- estimated_combinations <= exhaustive_threshold
  analysis_method <- if (use_exhaustive) "exhaustive" else "random_sampling"
  
  # ── Exhaustively enumerate all possible k-step merge paths ────────────────
  exhaustive_merging <- function() {
    # Enumerate every valid k-step merge path and store each step's statistic.
    
    message(sprintf("Starting exhaustive enumeration of all %d-order merge paths...", k_order))
    
    all_results <- data.frame(
      path_id = integer(),
      k       = integer(),
      i       = character(),
      j       = character(),
      new     = integer(),
      summary_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    path_counter <- 0L
    
    explore_merge_paths <- function(current_sf, current_k, path_history) {
      # Recursively extend one partial merge path until the requested order is
      # reached or no adjacent merge candidates remain.
      
      if (current_k > k_order) {
        # Target order reached; record the complete path
        path_counter <<- path_counter + 1L
        for (step in path_history) {
          all_results <<- rbind(all_results, data.frame(
            path_id = path_counter,
            k       = step$k,
            i       = step$i,
            j       = step$j,
            new     = step$new,
            summary_value = step$summary_value,
            stringsAsFactors = FALSE
          ))
        }
        if (path_counter %% 100 == 0) {
          message(sprintf("Explored %d complete paths so far", path_counter))
        }
        return()
      }
      
      G          <- poly2nb(current_sf, queen = FALSE)
      candidates <- which(sapply(G, length) > 0)
      
      if (length(candidates) == 0) {
        warning(sprintf("Path %d has no mergeable regions at order %d", path_counter + 1L, current_k))
        return()
      }
      
      # Iterate over all adjacent pairs; require i < j to avoid duplicates
      for (i_index in candidates) {
        for (j_index in G[[i_index]]) {
          if (i_index >= j_index) next
          
          id_i    <- current_sf$ID[i_index]
          id_j    <- current_sf$ID[j_index]
          
          temp_sf <- tryCatch(
            merging_operation(current_sf, id_i, id_j, field_rules),
            error = function(e) {
              warning(sprintf("Exhaustive merge failed; skipping this branch: %s", e$message))
              NULL
            }
          )
          if (is.null(temp_sf)) next
          
          new_id <- max(temp_sf$ID, na.rm = TRUE)
          
          current_value <- tryCatch({
            call_summary_function(temp_sf, summary_field, summary_function)
          }, error = function(e) {
            warning(paste0("summary_function failed (path ", path_counter + 1L,
                           ", k=", current_k, "): ", e$message))
            NA
          })
          
          explore_merge_paths(
            temp_sf,
            current_k + 1L,
            c(path_history, list(list(
              k     = current_k,
              i     = as.character(id_i),
              j     = as.character(id_j),
              new   = new_id,
              summary_value = current_value
            )))
          )
        }
      }
    }
    
    explore_merge_paths(sf, 1L, list())
    message(sprintf("Exhaustive enumeration finished; %d paths explored", path_counter))
    return(all_results)
  }
  
  # ── Randomly sample n_iterations k-step merge paths ───────────────────────
  random_merging <- function() {
    # Sample k-step merge paths when exhaustive enumeration is too expensive.
    
    message(sprintf("Starting random sampling of %d-order merging, %d iterations...", k_order, n_iterations))
    
    sf_local               <- sf
    field_rules_local      <- field_rules
    summary_field_local    <- summary_field
    summary_function_local <- summary_function
    merge_op               <- get("merging_operation", mode = "function", inherits = TRUE)
    call_summary_fn        <- get("call_summary_function", mode = "function", inherits = TRUE)
    
    run_one_merge_path <- function(iter) {
      # Run one independent random merge path and return its step log.
      temp_sf       <- sf_local
      iter_rows     <- vector("list", k_order)
      iter_valid    <- TRUE
      attempt_error <- NULL
      
      for (merge_level in seq_len(k_order)) {
        
        G          <- poly2nb(temp_sf, queen = FALSE)
        candidates <- which(sapply(G, length) > 0)
        
        if (length(candidates) == 0) {
          iter_valid <- FALSE
          break
        }
        
        # Randomly pick a pair of adjacent units
        i_index <- sample(candidates, 1)
        j_index <- G[[i_index]][sample(length(G[[i_index]]), 1)]
        id_i    <- temp_sf$ID[i_index]
        id_j    <- temp_sf$ID[j_index]
        
        temp_sf <- tryCatch(
          merge_op(temp_sf, id_i, id_j, field_rules_local),
          error = function(e) {
            attempt_error <<- sprintf(
              "merge failed at iteration %d, k=%d, IDs=%s/%s: %s",
              iter, merge_level, id_i, id_j, e$message
            )
            NULL
          }
        )
        
        if (is.null(temp_sf)) {
          iter_valid <- FALSE
          break
        }
        
        new_id <- max(temp_sf$ID, na.rm = TRUE)
        
        current_value <- tryCatch({
          call_summary_fn(temp_sf, summary_field_local, summary_function_local)
        }, error = function(e) {
          attempt_error <<- paste0("summary_function failed (iteration ", iter,
                                   ", k=", merge_level, "): ", e$message)
          NA
        })
        
        iter_rows[[merge_level]] <- data.frame(
          path_id = iter,
          k       = merge_level,
          i       = as.character(id_i),
          j       = as.character(id_j),
          new     = new_id,
          summary_value = current_value,
          stringsAsFactors = FALSE
        )
      }
      
      if (!iter_valid) {
        if (is.null(attempt_error)) {
          attempt_error <- sprintf(
            "no mergeable candidates at iteration %d before reaching k=%d",
            iter, k_order
          )
        }
        return(list(error = attempt_error))
      }
      do.call(rbind, iter_rows)
    }
    
    iter_results <- run_monte_carlo(
      X              = seq_len(n_iterations),
      FUN            = run_one_merge_path,
      parallel       = parallel,
      n_workers      = n_workers,
      random_seed    = random_seed,
      export_names   = summary_function_export_names(summary_function),
      export_env     = environment(summary_function),
      progress_label = "Merging sensitivity"
    )
    
    error_messages <- vapply(
      Filter(function(x) is.list(x) && !is.null(x$error), iter_results),
      function(x) x$error,
      character(1)
    )
    valid_results <- Filter(is.data.frame, iter_results)
    skipped       <- n_iterations - length(valid_results)
    results       <- if (length(valid_results) > 0L) {
      do.call(rbind, valid_results)
    } else {
      data.frame(
        path_id = integer(),
        k       = integer(),
        i       = character(),
        j       = character(),
        new     = integer(),
        summary_value = numeric(),
        stringsAsFactors = FALSE
      )
    }
    
    message(sprintf("Random sampling finished: %d successful, %d discarded", n_iterations - skipped, skipped))
    if (length(error_messages) > 0L) {
      warning(sprintf(
        "First merging failure messages:\n%s",
        paste(utils::head(unique(error_messages), 5L), collapse = "\n")
      ))
    }
    return(results)
  }
  
  # ── Run the analysis ──────────────────────────────────────────────────────
  if (use_exhaustive) {
    message(sprintf("Combination count (%.2e) is below the threshold (%.2e); using exhaustive method",
                    estimated_combinations, exhaustive_threshold))
    if (isTRUE(parallel)) {
      message("Exhaustive merge enumeration is recursive and runs serially; the parallel setting applies only to random sampling.")
    }
    distribution <- exhaustive_merging()
  } else {
    message(sprintf("Combination count (%.2e) exceeds the threshold (%.2e); using random sampling (n=%d)",
                    estimated_combinations, exhaustive_threshold, n_iterations))
    distribution <- random_merging()
  }
  
  # ── Assemble the final result object ──────────────────────────────────────
  result <- list(
    origin_value    = origin_value,
    k_order         = k_order,
    analysis_method = analysis_method,
    distribution    = distribution
  )
  
  # ── Save the result as an .rds file ───────────────────────────────────────
  if (is.null(save_path)) {
    save_path <- sprintf("merging_sensitivity_k%d_%s_%s.rds",
                         k_order,
                         analysis_method,
                         format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Make sure the target directory exists
  save_dir <- dirname(save_path)
  if (nzchar(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  tryCatch({
    saveRDS(result, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
  }, error = function(e) {
    warning(sprintf("Failed to save .rds file: %s", e$message))
  })
  
  result$save_path <- save_path
  return(result)
}

generate_smooth_guided_path <- function(start, end) {
  # Create a random smooth polyline between two boundary points.
  # The line is used as a stochastic cutting path for polygon splitting.
  
  # Number of intermediate steps follows a Poisson distribution
  n_steps <- rpois(1, lambda = 10)
  
  # Direction vector and its length from start to end
  direction <- end - start
  total_d <- sqrt(sum(direction^2))
  direction_unit <- direction / total_d
  
  # Perpendicular unit vector, used to add lateral offsets
  perp <- c(-direction_unit[2], direction_unit[1])
  
  # Sample positions along the segment (in [0, 1]) and sort them
  t_vec <- sort(c(0, runif(n_steps), 1))
  
  # Maximum allowed offset at each position; tapers to 0 at the two endpoints
  # so the path remains anchored at start and end
  r_max_vec <- total_d * pmin(t_vec, 1 - t_vec)
  r_vec <- runif(length(t_vec), min = -r_max_vec, max = r_max_vec)
  
  # Base points along the straight line plus the lateral offsets
  base_points <- matrix(start, nrow = length(t_vec), ncol = 2, byrow = TRUE) + 
    t_vec %*% t(direction)
  offsets <- r_vec %*% t(perp)
  path <- base_points + offsets
  return(path)
}

splitting_operation <- function(
    sf,
    target_id,
    field_rules = NULL,
    max_iter    = 100
) {
  # Split one selected areal unit into two parts and allocate attributes to
  # the new parts according to the supplied field rules.

  if (!requireNamespace("lwgeom", quietly = TRUE)) {
    stop("Package 'lwgeom' is required for vector splitting operations.")
  }
  
  # ---- Input validation and target lookup -----------------------------
  if (!inherits(sf, "sf"))       stop("Input must be an sf object")
  if (!target_id %in% sf$ID)     stop(paste0("Target ID does not exist: ", target_id))
  
  target_index <- which(sf$ID == target_id)
  target       <- sf[target_index, ] |> st_make_valid()
  
  attempts <- 0
  parts    <- NULL
  
  # ---- Generate a valid two-part split -------------------------------
  while (is.null(parts) || length(parts) != 2) {
    attempts <- attempts + 1
    if (attempts > max_iter) {
      stop("splitting_operation: exceeded maximum number of attempts; unable to split")
    }
    
    # Randomly pick two points on the boundary and generate a smooth guided
    # path between them to be used as the cutting line
    boundary <- st_cast(st_boundary(st_geometry(target)), "LINESTRING")
    coords   <- st_coordinates(boundary)
    
    if (nrow(coords) < 4) {
      stop("Insufficient boundary points for a valid split; check whether the polygon is degenerate")
    }
    
    idx        <- sample(1:nrow(coords), 2)
    p1         <- coords[idx[1], 1:2]
    p2         <- coords[idx[2], 1:2]
    path       <- generate_smooth_guided_path(p1, p2)
    split_line <- st_sfc(st_linestring(path), crs = st_crs(sf))
    
    parts <- tryCatch({
      lwgeom::st_split(st_geometry(target), split_line) |>
        st_collection_extract("POLYGON") |>
        st_make_valid()
    }, error = function(e) NULL)
  }
  
  # ---- Compute area shares for attribute allocation -------------------
  # Compute area ratios, used by the area_weighted rule
  areas      <- st_area(parts)
  ratios     <- as.numeric(areas / sum(areas))
  base_attrs <- sf[target_index, , drop = FALSE] |> st_drop_geometry()
  
  # ---- Build the two new rows and update their attributes -------------
  new_parts <- lapply(1:2, function(i) {
    attrs <- base_attrs
    
    if (!is.null(field_rules)) {
      
      # ── Step 1: handle simple-rule fields ─────────────────────────────
      # Simple rules: split value is a string ("area_weighted" or "keep").
      # These must be processed before derived fields, because derived
      # fields rely on the values produced in this step.
      for (field_name in names(field_rules)) {
        if (!field_name %in% names(sf)) next
        rule <- field_rules[[field_name]]$split
        if (is.list(rule)) next  # derived fields are deferred to step 2
        
        orig_val <- sf[[field_name]][target_index]
        
        attrs[[field_name]] <- switch(
          rule,
          area_weighted = orig_val * ratios[i],  # allocate by area share
          keep          = orig_val,              # keep the original value
          orig_val  # fall back to the original value for unknown rules
        )
      }
      
      # ── Step 2: handle derived fields ─────────────────────────────────
      # Derived fields: split value is a list with `numerator` and
      # `denominator` sub-fields. Ratio-type fields (e.g. percentages)
      # cannot be split directly by area; they must be recomputed from
      # numerator / denominator after the split, so the resulting values
      # remain internally consistent (e.g. insurance coverage =
      # insured population / total population).
      for (field_name in names(field_rules)) {
        if (!field_name %in% names(sf)) next
        rule <- field_rules[[field_name]]$split
        if (!is.list(rule)) next  # only handle derived fields here
        
        num_field <- rule$numerator
        den_field <- rule$denominator
        
        # Prefer values already updated in step 1 (post-split values)
        num_val <- if (num_field %in% names(attrs)) attrs[[num_field]] else sf[[num_field]][target_index]
        den_val <- if (den_field %in% names(attrs)) attrs[[den_field]] else sf[[den_field]][target_index]
        
        attrs[[field_name]] <- if (!is.na(den_val) && den_val != 0) {
          num_val / den_val
        } else {
          NA_real_
        }
      }
    }
    
    st_as_sf(cbind(attrs, geometry = parts[i]))
  })
  
  # ---- Assign IDs and rebuild the sf object ---------------------------
  # Assign new IDs; this matches the setdiff(after$ID, before$ID) tracking
  # logic used in splitting_sensitivity
  new_ids           <- max(sf$ID, na.rm = TRUE) + 1:2
  new_parts[[1]]$ID <- new_ids[1]
  new_parts[[2]]$ID <- new_ids[2]
  
  out_sf <- rbind(sf[-target_index, ], new_parts[[1]], new_parts[[2]])
  return(out_sf)
}

splitting_sensitivity <- function(sf,
                                  field_rules,
                                  summary_field,
                                  summary_function,
                                  k_order      = 1,
                                  n_iterations = 100,
                                  random_seed  = NULL,
                                  parallel     = FALSE,
                                  n_workers    = default_parallel_workers(),
                                  save_path    = NULL) {
  # Estimate how much a summary statistic changes after k sequential splits.
  # Split paths are sampled randomly because the search space is very large.
  
  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.function(summary_function))   stop("summary_function must be a function")
  summary_field <- validate_summary_fields(sf, summary_field)
  if (k_order < 1)                      stop("k_order must be >= 1")
  
  # ── Compute baseline statistic ────────────────────────────────────────────
  # Record the statistic on the original data as a reference baseline before
  # any split operation is performed
  origin_value <- tryCatch({
    call_summary_function(sf, summary_field, summary_function)
  }, error = function(e) {
    stop(paste0("summary_function failed on the initial data: ", e$message))
  })
  
  # ── Random sampling ───────────────────────────────────────────────────────
  # The number of possible split paths is huge (the candidate count grows
  # rapidly with k), so we do not enumerate exhaustively and always rely on
  # random sampling. n_iterations is the target number of *successful* paths,
  # not the total number of attempts.
  message(sprintf("Starting random sampling of %d-order splitting; target successes: %d...",
                  k_order, n_iterations))
  
  # Initialise the result table; each row corresponds to one step of one path
  results <- data.frame(
    path_id = integer(),   # path index (incremented by success order)
    k       = integer(),   # which split step within the path
    i       = character(), # ID of the region being split
    new1    = integer(),   # ID of the first new region produced
    new2    = integer(),   # ID of the second new region produced
    summary_value = numeric(), # summary statistic after the split
    stringsAsFactors = FALSE
  )
  
  successful   <- 0L
  attempts     <- 0L
  error_messages <- character(0)
  # Cap the total number of attempts to prevent infinite loops when split
  # conditions are very restrictive
  max_attempts <- max(n_iterations * 10L, n_iterations + 1000L)
  
  field_rules_local      <- field_rules
  summary_field_local    <- summary_field
  summary_function_local <- summary_function
  guided_path_fn         <- get("generate_smooth_guided_path", mode = "function", inherits = TRUE)
  split_op               <- get("splitting_operation", mode = "function", inherits = TRUE)
  call_summary_fn        <- get("call_summary_function", mode = "function", inherits = TRUE)

  # Make the splitter self-contained for PSOCK workers. Its body calls
  # generate_smooth_guided_path(), so bind that name in the splitter's
  # enclosing environment before run_monte_carlo() serializes the closure.
  generate_smooth_guided_path <- guided_path_fn
  environment(split_op) <- environment()

  run_one_split_attempt <- function(attempt_id) {
    # Try to generate one successful k-step split path from the original sf.
    temp_sf       <- sf
    iter_rows     <- vector("list", k_order)
    iter_valid    <- TRUE
    attempt_error <- NULL

    for (split_level in seq_len(k_order)) {
      if (nrow(temp_sf) == 0) {
        iter_valid    <- FALSE
        attempt_error <- "temporary sf became empty"
        break
      }

      i_index    <- sample(seq_len(nrow(temp_sf)), 1)
      id_i       <- temp_sf$ID[i_index]
      ids_before <- temp_sf$ID

      temp_sf <- tryCatch(
        split_op(temp_sf, id_i, field_rules_local),
        error = function(e) {
          attempt_error <<- paste0(
            "split failed at attempt ", attempt_id,
            ", k=", split_level,
            ", ID=", id_i,
            ": ", e$message
          )
          NULL
        }
      )

      if (is.null(temp_sf)) {
        iter_valid <- FALSE
        break
      }

      new_ids <- sort(setdiff(temp_sf$ID, ids_before))
      new1    <- if (length(new_ids) >= 1) new_ids[1] else NA_integer_
      new2    <- if (length(new_ids) >= 2) new_ids[2] else NA_integer_

      current_value <- tryCatch(
        call_summary_fn(temp_sf, summary_field_local, summary_function_local),
        error = function(e) {
          attempt_error <<- paste0(
            "summary_function failed at attempt ", attempt_id,
            ", k=", split_level,
            ": ", e$message
          )
          NA_real_
        }
      )

      iter_rows[[split_level]] <- data.frame(
        path_id = attempt_id,
        k       = split_level,
        i       = as.character(id_i),
        new1    = new1,
        new2    = new2,
        summary_value = current_value,
        stringsAsFactors = FALSE
      )
    }

    if (!iter_valid) {
      if (is.null(attempt_error)) attempt_error <- "unknown split failure"
      return(list(error = attempt_error))
    }

    do.call(rbind, iter_rows)
  }

  valid_results <- list()

  while (successful < n_iterations && attempts < max_attempts) {
    remaining_successes <- n_iterations - successful

    # Preserve the existing execution strategies: serial mode submits one
    # attempt at a time, while parallel mode oversamples candidate attempts
    # to compensate for failed paths.
    batch_size <- if (isTRUE(parallel)) {
      min(
        max(remaining_successes * 2L, n_workers),
        max_attempts - attempts
      )
    } else {
      1L
    }

    batch_ids  <- attempts + seq_len(batch_size)
    batch_seed <- if (is.null(random_seed)) NULL else random_seed + attempts

    attempt_results <- run_monte_carlo(
      X              = batch_ids,
      FUN            = run_one_split_attempt,
      parallel       = parallel,
      n_workers      = n_workers,
      random_seed    = batch_seed,
      export_names   = summary_function_export_names(summary_function),
      export_env     = environment(summary_function),
      progress_label = "Splitting sensitivity"
    )

    batch_errors <- vapply(
      Filter(function(x) is.list(x) && !is.null(x$error), attempt_results),
      function(x) x$error,
      character(1)
    )
    if (length(batch_errors) > 0L) {
      error_messages <- c(error_messages, batch_errors)
    }

    batch_valid <- Filter(is.data.frame, attempt_results)
    if (length(batch_valid) > 0L) {
      valid_results <- c(valid_results, batch_valid)
    }

    attempts   <- attempts + batch_size
    successful <- min(length(valid_results), n_iterations)

    if (!isTRUE(parallel) && length(batch_valid) > 0L &&
        successful > 0L && successful %% 10L == 0L) {
      message(sprintf(
        "Successfully completed %d/%d (total attempts: %d)",
        successful, n_iterations, attempts
      ))
    }
  }

  if (successful > 0L) {
    selected <- valid_results[seq_len(successful)]
    selected <- Map(function(df, path_id) {
      df$path_id <- path_id
      df
    }, selected, seq_len(successful))
    results <- do.call(rbind, selected)
  }
  
  # If the attempt cap is reached without meeting the target, raise a clear
  # warning rather than silently returning a deficient result
  if (successful < n_iterations) {
    warning(sprintf(
      "Reached the maximum attempt cap (%d); only %d/%d paths succeeded. Consider relaxing the split conditions or lowering k_order",
      max_attempts, successful, n_iterations
    ))
    
    if (length(error_messages) > 0L) {
      warning(sprintf(
        "First splitting failure messages:\n%s",
        paste(utils::head(unique(error_messages), 5L), collapse = "\n")
      ))
    }
  }
  
  message(sprintf("Random sampling finished: %d successes, %d total attempts",
                  successful, attempts))
  
  # ── Assemble the final result object ──────────────────────────────────────
  result <- list(
    origin_value = origin_value,  # baseline statistic on the original data
    k_order      = k_order,       # split order of this analysis
    distribution = results        # step-by-step record of all successful paths
  )
  
  # ── Save the result as an .rds file ───────────────────────────────────────
  if (is.null(save_path)) {
    save_path <- sprintf("splitting_sensitivity_k%d_%s.rds",
                         k_order,
                         format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Make sure the target directory exists
  save_dir <- dirname(save_path)
  if (nzchar(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  tryCatch({
    saveRDS(result, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
  }, error = function(e) {
    warning(sprintf("Failed to save .rds file: %s", e$message))
  })
  
  result$save_path <- save_path
  return(result)
}

continuous_reassignment_operation <- function(
    sf,
    alpha,
    field_rules,
    summary_field,
    summary_function,
    max_iter  = 1000,
    tol_ratio = 0.1
) {
  # Reassign a continuous amount of area by splitting a source region and
  # merging one split part into a neighboring target region.
  
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(sf, "sf"))              stop("Input must be an sf object")
  summary_field <- validate_summary_fields(sf, summary_field)
  if (!is.function(summary_function))   stop("summary_function must be a function")
  
  # ── Initialize state variables ──────────────────────────────────────────
  tol             <- tol_ratio * alpha   # absolute tolerance derived from alpha
  remaining_alpha <- alpha               # area budget still to be reassigned
  failure         <- 0L                  # cumulative failed attempts
  total_iter      <- 0L                  # total while-loop iterations
  step_count      <- 0L                  # number of successful reassignment steps
  source_blocked  <- c()                 # IDs that should no longer be picked as source/target
  current_sf      <- sf                  # working copy of the sf object
  step_log_list   <- list()              # per-step log entries (assembled at the end)
  
  # ── Compute baseline summary statistic ──────────────────────────────────
  stat_before <- tryCatch(
    call_summary_function(sf, summary_field, summary_function),
    error = function(e) stop(paste0("summary_function failed: ", e$message))
  )
  current_stat       <- stat_before
  termination_reason <- NA_character_
  
  # ── Main reassignment loop ──────────────────────────────────────────────
  # The loop terminates when any of the following holds:
  #   (a) too many cumulative failures (`failure  >= max_iter`),
  #   (b) too many overall iterations  (`total_iter >= max_iter`),
  #   (c) the remaining area budget is within tolerance of zero, i.e.
  #       enough area has been reassigned to satisfy `alpha` up to `tol`.
  while (failure < max_iter && total_iter < max_iter && abs(remaining_alpha) > tol) {
    
    total_iter <- total_iter + 1L
    
    # ---- Pick a source polygon from the unblocked candidate pool ------
    candidates <- setdiff(current_sf$ID, source_blocked)
    if (length(candidates) == 0L) {
      termination_reason <- "candidate set is empty"
      break
    }
    
    source_id  <- sample(candidates, 1)
    source_row <- which(current_sf$ID == source_id)
    
    # ---- Identify topological neighbors of the source -----------------
    # Only direct neighbors (sharing a boundary) are valid merge targets;
    # this preserves spatial contiguity in the resulting zoning.
    pre_neighbor_idx <- st_touches(
      current_sf[source_row, ], current_sf, sparse = TRUE
    )[[1]]
    pre_neighbor_ids <- current_sf$ID[pre_neighbor_idx]
    
    # An isolated source polygon can never produce a valid reassignment,
    # so block it permanently to avoid resampling it.
    if (length(pre_neighbor_ids) == 0L) {
      source_blocked <- union(source_blocked, source_id)
      next
    }
    
    # ---- Split the source polygon into sub-parts ----------------------
    split_result <- tryCatch(
      splitting_operation(current_sf, source_id, field_rules),
      error = function(e) NULL
    )
    if (is.null(split_result)) {
      failure <- failure + 1L
      next
    }
    
    # Newly created sub-parts are exactly the IDs that did not exist
    # before the split.
    new_ids   <- setdiff(split_result$ID, current_sf$ID)
    new_rows  <- which(split_result$ID %in% new_ids)
    new_geoms <- split_result[new_rows, ]
    
    # Randomize the order in which sub-parts are tried for merging, so
    # the algorithm does not systematically prefer the first fragment
    # produced by the splitter (which would bias the area distribution).
    try_order <- sample(seq_len(nrow(new_geoms)))
    
    # ---- Try to merge one sub-part into a neighboring polygon ---------
    success <- FALSE
    for (i in try_order) {
      part_to_merge <- new_geoms[i, ]
      part_area     <- as.numeric(st_area(part_to_merge))
      
      # Skip parts that would overshoot the remaining area budget.
      if (part_area > remaining_alpha + tol) next
      
      # Valid merge targets must (1) be a pre-existing neighbor of the
      # source, and (2) still exist in the post-split sf (a neighbor
      # might have been replaced if it was also affected by the split).
      valid_target_ids <- intersect(pre_neighbor_ids, split_result$ID)
      if (length(valid_target_ids) == 0L) next
      
      target_id <- sample(valid_target_ids, 1)
      part_id   <- part_to_merge$ID
      
      # Attempt the merge; on failure, try the next sub-part.
      merged_result <- tryCatch(
        merging_operation(split_result, part_id, target_id, field_rules),
        error = function(e) NULL
      )
      if (is.null(merged_result)) next
      
      # Recompute the summary statistic on the updated sf.
      new_stat <- tryCatch(
        call_summary_function(merged_result, summary_field, summary_function),
        error = function(e) NA_real_
      )
      
      # ---- Log this successful step ----------------------------------
      step_count <- step_count + 1L
      step_log_list[[step_count]] <- data.frame(
        step        = step_count,
        source_id   = as.character(source_id),
        target_id   = as.character(target_id),
        part_area   = part_area,
        stat_before = current_stat,
        stat_after  = new_stat,
        delta_stat  = new_stat - current_stat,
        stringsAsFactors = FALSE
      )
      
      # ---- Update state for the next iteration -----------------------
      new_merged_id   <- max(merged_result$ID, na.rm = TRUE)
      current_stat    <- new_stat
      remaining_alpha <- remaining_alpha - part_area
      # Block the newly merged unit so it cannot be picked as a source in
      # the next iteration. The old target_id no longer exists after
      # merging_operation(), which replaces both input rows with a new ID.
      source_blocked  <- union(source_blocked, new_merged_id)
      current_sf      <- merged_result
      success         <- TRUE
      break
    }
    
    # If no sub-part could be merged this round, count it as a failure.
    if (!success) failure <- failure + 1L
  }
  
  # ── Determine termination reason ────────────────────────────────────────
  if (is.na(termination_reason)) {
    termination_reason <- if (abs(remaining_alpha) <= tol) "success" else "max_iter reached"
  }
  
  # ── Final statistics on the resulting sf ────────────────────────────────
  converged  <- abs(remaining_alpha) <= tol
  stat_after <- tryCatch(
    call_summary_function(current_sf, summary_field, summary_function),
    error = function(e) NA_real_
  )
  
  # ── Assemble the step log ───────────────────────────────────────────────
  # When no successful step was recorded, return an empty data.frame whose
  # schema matches the populated case exactly. This lets downstream code
  # rbind results across runs without hitting column-mismatch errors.
  step_log <- if (length(step_log_list) > 0) {
    do.call(rbind, step_log_list)
  } else {
    data.frame(
      step = integer(0), source_id = character(0), target_id = character(0),
      part_area = numeric(0), stat_before = numeric(0), stat_after = numeric(0),
      delta_stat = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  # ── Return result bundle ────────────────────────────────────────────────
  return(list(
    converged          = converged,
    termination_reason = termination_reason,
    stat_before        = stat_before,
    stat_after         = stat_after,
    stat_diff          = stat_after - stat_before,
    remaining_alpha    = remaining_alpha,
    total_failure      = failure,
    updated_sf         = current_sf,
    step_log           = step_log
  ))
}

continuous_reassignment_sensitivity <- function(
    sf,
    field_rules,
    summary_field,
    summary_function,
    alpha,
    n_iterations = 100,
    max_iter     = 1000,
    tol_ratio    = 0.1,
    random_seed  = NULL,
    parallel     = FALSE,
    n_workers    = default_parallel_workers(),
    save_path    = NULL
) {
  # Run repeated continuous reassignment simulations and summarize the
  # resulting distribution of the selected statistic.
  
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(sf, "sf"))
    stop("Input must be an sf object")
  if (!is.function(summary_function))
    stop("summary_function must be a function")
  summary_field <- validate_summary_fields(sf, summary_field)
  if (!is.numeric(alpha) || alpha <= 0)
    stop("alpha must be a positive number")
  if (!is.numeric(n_iterations) || n_iterations < 1)
    stop("n_iterations must be a positive integer")
  if (!is.numeric(max_iter) || max_iter < 1)
    stop("max_iter must be a positive integer")
  if (!is.numeric(tol_ratio) || tol_ratio <= 0 || tol_ratio >= 1)
    stop("tol_ratio must be in the open interval (0, 1)")
  
  # Force lazily evaluated arguments before they are captured by PSOCK
  # workers. This avoids expressions such as
  # continuous_reassignment_plan[[fig_id]]$alpha being evaluated later inside
  # workers where the plan object does not exist.
  alpha        <- as.numeric(alpha)
  n_iterations <- as.integer(n_iterations)
  max_iter     <- as.integer(max_iter)
  tol_ratio    <- as.numeric(tol_ratio)
  
  # ── Compute baseline statistic ──────────────────────────────────────────
  # Recorded once on the untouched data so every iteration's result can be
  # compared against a single, stable reference value.
  origin_value <- tryCatch(
    call_summary_function(sf, summary_field, summary_function),
    error = function(e) stop(paste0("summary_function failed: ", e$message))
  )
  message(sprintf("Baseline statistic: %.6f", origin_value))
  
  # ── Initialize result container ─────────────────────────────────────────
  attempt_results <- data.frame(
    iteration       = integer(0),
    summary_value   = numeric(0),
    converged       = logical(0),
    termination_reason = character(0),
    total_failure   = integer(0),
    remaining_alpha = numeric(0),
    successful_steps = integer(0),
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "Starting continuous reassignment sensitivity analysis: %d iterations",
    n_iterations
  ))
  
  # ── Main loop ───────────────────────────────────────────────────────────
  # Each iteration starts from the original `sf` (not the previous result),
  # so the n_iterations runs are mutually independent Monte-Carlo samples.
  sf_local               <- sf
  field_rules_local      <- field_rules
  summary_field_local    <- summary_field
  summary_function_local <- summary_function
  continuous_op          <- get("continuous_reassignment_operation", mode = "function", inherits = TRUE)
  guided_path_fn         <- get("generate_smooth_guided_path", mode = "function", inherits = TRUE)
  split_op               <- get("splitting_operation", mode = "function", inherits = TRUE)
  merging_operation      <- get("merging_operation", mode = "function", inherits = TRUE)
  call_summary_function  <- get("call_summary_function", mode = "function", inherits = TRUE)
  validate_summary_fields <- get("validate_summary_fields", mode = "function", inherits = TRUE)

  # Make the splitter self-contained before the worker closure is serialized.
  # splitting_operation() calls generate_smooth_guided_path(), which is not
  # otherwise available in a fresh PSOCK worker's global environment.
  generate_smooth_guided_path <- guided_path_fn
  environment(split_op) <- environment()
  splitting_operation <- split_op
  environment(continuous_op) <- environment()
  
  run_one_continuous <- function(run) {
    # Run one independent continuous reassignment simulation.
    
    res <- tryCatch(
      continuous_op(
        sf               = sf_local,
        alpha            = alpha,
        field_rules      = field_rules_local,
        summary_field    = summary_field_local,
        summary_function = summary_function_local,
        max_iter         = max_iter,
        tol_ratio        = tol_ratio
      ),
      error = function(e) {
        return(list(error = sprintf("Run %d failed: %s", run, e$message)))
      }
    )
    
    if (is.list(res) && !is.null(res$error)) return(res)
    
    data.frame(
      iteration       = run,
      summary_value   = res$stat_after,
      converged       = isTRUE(res$converged),
      termination_reason = as.character(res$termination_reason),
      total_failure   = res$total_failure,
      remaining_alpha = res$remaining_alpha,
      successful_steps = nrow(res$step_log),
      stringsAsFactors = FALSE
    )
  }
  
  run_results <- run_monte_carlo(
    X              = seq_len(n_iterations),
    FUN            = run_one_continuous,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = summary_function_export_names(summary_function),
    export_env     = environment(summary_function),
    progress_label = "Continuous reassignment sensitivity"
  )
  
  error_messages <- vapply(
    Filter(function(x) is.list(x) && !is.null(x$error), run_results),
    function(x) x$error,
    character(1)
  )
  returned_attempts <- Filter(is.data.frame, run_results)
  attempt_results <- if (length(returned_attempts) > 0L) {
    do.call(rbind, returned_attempts)
  } else {
    attempt_results
  }

  # Only completed area-budget paths define the sensitivity distribution.
  # Non-converged attempts remain available in diagnostics but must not be
  # plotted as if their unchanged or partially changed values were valid.
  converged_mask <- attempt_results$converged & is.finite(attempt_results$summary_value)
  results <- attempt_results[
    converged_mask,
    c("iteration", "summary_value", "total_failure", "remaining_alpha"),
    drop = FALSE
  ]

  n_converged     <- nrow(results)
  n_nonconverged  <- sum(!attempt_results$converged)
  n_worker_errors <- length(error_messages)

  message(sprintf(
    "Done: %d/%d converged runs collected (%d non-converged, %d worker errors)",
    n_converged,
    n_iterations,
    n_nonconverged,
    n_worker_errors
  ))

  if (n_nonconverged > 0L) {
    failed_reasons <- attempt_results$termination_reason[!attempt_results$converged]
    reason_counts <- table(failed_reasons, useNA = "ifany")
    reason_summary <- paste0(
      names(reason_counts), "=", as.integer(reason_counts),
      collapse = ", "
    )
    warning(sprintf(
      "%d continuous reassignment runs did not converge and were excluded from the distribution. Reasons: %s",
      n_nonconverged,
      reason_summary
    ))
  }

  if (length(error_messages) > 0L) {
    warning(sprintf(
      "First continuous reassignment failure messages:\n%s",
      paste(utils::head(unique(error_messages), 5L), collapse = "\n")
    ))
  }
  
  out <- list(
    origin_value = origin_value,
    distribution = results,
    summary_field = summary_field,
    diagnostics  = list(
      requested_iterations = n_iterations,
      converged_runs       = n_converged,
      nonconverged_runs    = n_nonconverged,
      worker_errors        = n_worker_errors,
      attempts             = attempt_results[
        , c(
          "iteration",
          "converged",
          "termination_reason",
          "total_failure",
          "remaining_alpha",
          "successful_steps"
        ),
        drop = FALSE
      ],
      error_messages       = unname(error_messages)
    )
  )
  
  # ── Persist to disk (optional) ──────────────────────────────────────────
  if (!is.null(save_path)) {
    save_dir <- dirname(save_path)
    if (nzchar(save_dir) && !dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
  }
  
  return(out)
}

discrete_reassignment_operation <- function(
    sf_fine,
    sf_coarse,
    i,
    neighbor_list,
    fine_id,
    coarse_id,
    field_rules
) {
  # Move one fine unit from its current coarse unit to a neighboring coarse
  # unit, then update coarse-level attributes and geometries.
  
  # ── Locate source coarse polygon ────────────────────────────────────────
  # `coarse_id` is the polygon that currently owns the fine unit `i`.
  # If it cannot be found in sf_coarse, abort by returning NULL so the
  # caller can decide how to handle the missed reassignment.
  coarse_row <- which(sf_coarse$ID == coarse_id)
  if (length(coarse_row) == 0L) return(NULL)
  
  # ── Pick a neighboring coarse polygon as the reassignment target ────────
  # Only topological neighbors (precomputed via st_touches) are valid
  # targets — this preserves spatial contiguity of the resulting zoning.
  neighbor_idx <- neighbor_list[[coarse_row]]
  valid_nb_ids <- sf_coarse$ID[neighbor_idx]
  if (length(valid_nb_ids) == 0L) return(NULL)   # isolated coarse polygon
  
  target_id  <- sample(valid_nb_ids, 1)
  target_row <- which(sf_coarse$ID == target_id)
  
  # ── Reassign the fine unit at the lookup table level ────────────────────
  # The fine layer is unchanged geometrically; only its parent ID flips.
  sf_fine$coarse_ID[i] <- target_id
  
  # Recompute fine-row memberships *after* the flip, so subsequent
  # aggregations reflect the new parent assignment.
  source_fine_rows <- which(sf_fine$coarse_ID == coarse_id)
  target_fine_rows <- which(sf_fine$coarse_ID == target_id)
  
  # ── Recompute coarse-level fields per merge rule ────────────────────────
  # Two rule shapes are supported:
  #   (a) "sum"  → re-aggregate from the fine layer (extensive variable).
  #   (b) list(numerator, denominator) → recompute as a ratio from
  #       already-updated coarse fields (intensive variable / rate).
  # Note: ratio rules must be evaluated *after* their numerator and
  # denominator have been re-summed — the caller is expected to order
  # `field_rules` so sums precede ratios that depend on them.
  for (fn in names(field_rules)) {
    rule       <- field_rules[[fn]]
    merge_rule <- rule$merge
    
    if (is.character(merge_rule) && merge_rule == "sum") {
      if (!fn %in% names(sf_fine)) next
      
      # Source may now be empty (all its fine units flipped away);
      # in that case its aggregate becomes NA rather than 0 to
      # signal "no underlying data" instead of "true zero".
      sf_coarse[[fn]][coarse_row] <- if (length(source_fine_rows) > 0L)
        sum(sf_fine[[fn]][source_fine_rows], na.rm = TRUE) else NA_real_
      
      sf_coarse[[fn]][target_row] <- if (length(target_fine_rows) > 0L)
        sum(sf_fine[[fn]][target_fine_rows], na.rm = TRUE) else NA_real_
      
    } else if (is.list(merge_rule)) {
      num_fn <- merge_rule$numerator
      den_fn <- merge_rule$denominator
      if (!all(c(num_fn, den_fn) %in% names(sf_coarse))) next
      
      # Guard against division by zero / NA denominators.
      den_s <- sf_coarse[[den_fn]][coarse_row]
      sf_coarse[[fn]][coarse_row] <-
        if (!is.na(den_s) && den_s != 0)
          sf_coarse[[num_fn]][coarse_row] / den_s else NA_real_
      
      den_t <- sf_coarse[[den_fn]][target_row]
      sf_coarse[[fn]][target_row] <-
        if (!is.na(den_t) && den_t != 0)
          sf_coarse[[num_fn]][target_row] / den_t else NA_real_
    }
  }
  
  # ── Rebuild coarse geometries from their fine-unit members ──────────────
  # Source geometry is only rebuilt if it still owns at least one fine
  # unit; otherwise its old geometry is left in place and the caller
  # should treat it as "logically empty" (its summary fields are NA).
  if (length(source_fine_rows) > 0L) {
    st_geometry(sf_coarse)[coarse_row] <-
      st_geometry(st_union(sf_fine[source_fine_rows, ]))
  }
  # The target always gains a unit, so its geometry is always rebuilt.
  st_geometry(sf_coarse)[target_row] <-
    st_geometry(st_union(sf_fine[target_fine_rows, ]))
  
  return(list(sf_fine = sf_fine, sf_coarse = sf_coarse))
}

discrete_reassignment_sensitivity_area <- function(
    finer_sf,
    coarser_sf,
    field_rules,
    n_iterations,
    alpha,
    summary_field,
    summary_function,
    tol_ratio  = 0.1,
    keep_maps  = FALSE,
    random_seed = NULL,
    parallel   = FALSE,
    n_workers  = default_parallel_workers(),
    save_path  = NULL        # If non-NULL, results are saved as a .rds file
) {
  # Run discrete reassignment simulations constrained by a total reassigned
  # area budget and summarize the statistic after each simulation.
  
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(finer_sf,   "sf")) stop("finer_sf must be an sf object")
  if (!inherits(coarser_sf, "sf")) stop("coarser_sf must be an sf object")
  summary_field <- validate_summary_fields(coarser_sf, summary_field, "coarser_sf")
  if (!is.function(summary_function))
    stop("summary_function must be a function")
  if (!is.numeric(alpha) || alpha <= 0)
    stop("alpha must be a positive number")
  if (!is.numeric(tol_ratio) || tol_ratio <= 0 || tol_ratio >= 1)
    stop("tol_ratio must be in the open interval (0, 1)")
  if (!is.numeric(n_iterations) || n_iterations < 1)
    stop("n_iterations must be a positive integer")
  if (!is.null(save_path) && !is.character(save_path))
    stop("save_path must be a character string path or NULL")
  
  n_iterations <- as.integer(n_iterations)
  delta_alpha  <- tol_ratio * alpha
  # ── Baseline statistic & precomputation ─────────────────────────────────
  origin_value  <- call_summary_function(coarser_sf, summary_field, summary_function)
  neighbor_list <- st_touches(coarser_sf, sparse = TRUE)
  fine_areas    <- as.numeric(st_area(finer_sf))
  error_messages <- character(0)
  
  field_rules_local      <- field_rules
  summary_field_local    <- summary_field
  summary_function_local <- summary_function
  keep_maps_local        <- keep_maps
  neighbor_list_local    <- neighbor_list
  fine_areas_local       <- fine_areas
  reassignment_op        <- get("discrete_reassignment_operation", mode = "function", inherits = TRUE)
  call_summary_fn        <- get("call_summary_function", mode = "function", inherits = TRUE)

  run_one_discrete_area <- function(iter) {
    # Run one area-budgeted discrete reassignment simulation.
    budget <- runif(
      1,
      min = alpha - delta_alpha,
      max = alpha + delta_alpha
    )

    current_fine   <- finer_sf
    current_coarse <- coarser_sf
    remaining      <- budget
    available_ids  <- current_fine$ID
    any_reassigned <- FALSE
    attempt_error   <- NULL

    while (remaining > 0 && length(available_ids) > 0L) {
      fine_id  <- sample(available_ids, 1)
      fine_row <- which(current_fine$ID == fine_id)
      area_i   <- fine_areas_local[fine_row]

      if (area_i > remaining) {
        available_ids <- setdiff(available_ids, fine_id)
        next
      }

      coarse_id <- current_fine$coarse_ID[fine_row]

      result <- tryCatch(
        reassignment_op(
          sf_fine       = current_fine,
          sf_coarse     = current_coarse,
          i             = fine_row,
          neighbor_list = neighbor_list_local,
          fine_id       = fine_id,
          coarse_id     = coarse_id,
          field_rules   = field_rules_local
        ),
        error = function(e) {
          attempt_error <<- sprintf(
            "area reassignment failed at iteration %d, fine ID=%s: %s",
            iter, fine_id, e$message
          )
          NULL
        }
      )

      if (is.null(result)) {
        available_ids <- setdiff(available_ids, fine_id)
        next
      }

      remaining      <- remaining - area_i
      current_fine   <- result$sf_fine
      current_coarse <- result$sf_coarse
      available_ids  <- setdiff(available_ids, fine_id)
      any_reassigned <- TRUE
    }

    if (!any_reassigned) {
      if (is.null(attempt_error)) {
        attempt_error <- sprintf(
          "no eligible fine unit could be reassigned within the area budget in iteration %d",
          iter
        )
      }
      return(list(error = attempt_error))
    }

    list(
      stat_after = call_summary_fn(
        current_coarse,
        summary_field_local,
        summary_function_local
      ),
      updated_sf = if (keep_maps_local) current_coarse else NULL
    )
  }

  results <- run_monte_carlo(
    X              = seq_len(n_iterations),
    FUN            = run_one_discrete_area,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = summary_function_export_names(summary_function),
    export_env     = environment(summary_function),
    progress_label = "Discrete reassignment sensitivity by area"
  )
  error_messages <- vapply(
    Filter(function(x) is.list(x) && !is.null(x$error), results),
    function(x) x$error,
    character(1)
  )
  results <- Filter(function(x) is.list(x) && !is.null(x$stat_after), results)
  success_count <- length(results)

  # ── Assemble distribution & deviation summary ───────────────────────────
  valid <- Filter(function(x) is.list(x) && !is.null(x$stat_after), results)
  distribution <- if (length(valid) > 0L) {
    data.frame(summary_value = vapply(valid, function(x) x$stat_after, numeric(1)))
  } else {
    data.frame(summary_value = numeric(0))
  }
  devs         <- distribution$summary_value / origin_value
  
  if (success_count == 0L) {
    warning("No successful discrete area reassignment simulations were collected.")
  }
  if (length(error_messages) > 0L) {
    warning(sprintf(
      "First discrete area failure messages:\n%s",
      paste(utils::head(unique(error_messages), 5L), collapse = "\n")
    ))
  }
  
  out <- list(
    origin_value  = origin_value,
    distribution  = distribution,
    alpha         = alpha,
    tol_ratio     = tol_ratio,
    n_iterations  = n_iterations,
    n_success     = success_count,
    deviations    = devs,
    sd_deviation  = sd(devs),
    maps          = if (keep_maps) lapply(valid, `[[`, "updated_sf") else NULL
  )
  
  # ── Persist to disk (optional) ──────────────────────────────────────────
  # The save step is placed at the very end so that even if writing to disk
  # fails, the function can still return the result object to the caller.
  # The parent directory is created automatically; dir.create is silent if
  # the directory already exists.
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
  }
  
  return(out)
}


discrete_reassignment_sensitivity_regions <- function(
    finer_sf,
    coarser_sf,
    field_rules,
    k_regions,
    n_iterations,
    summary_function,
    summary_field,
    random_seed = NULL,
    parallel    = FALSE,
    n_workers   = default_parallel_workers(),
    save_path = NULL          # If non-NULL, results are saved as a .rds file
) {
  # Run discrete reassignment simulations constrained by the number of fine
  # regions moved and summarize the statistic after each simulation.
  
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(finer_sf,   "sf")) stop("finer_sf must be an sf object")
  if (!inherits(coarser_sf, "sf")) stop("coarser_sf must be an sf object")
  summary_field <- validate_summary_fields(coarser_sf, summary_field, "coarser_sf")
  if (!is.function(summary_function))
    stop("summary_function must be a function")
  if (!is.numeric(k_regions) || k_regions < 1)
    stop("k_regions must be a positive integer")
  if (!is.numeric(n_iterations) || n_iterations < 1)
    stop("n_iterations must be a positive integer")
  if (!is.null(save_path) && !is.character(save_path))
    stop("save_path must be a character string path or NULL")
  
  k_regions    <- as.integer(k_regions)
  n_iterations <- as.integer(n_iterations)
  # ── Baseline statistic & precomputation ─────────────────────────────────
  origin_value  <- call_summary_function(coarser_sf, summary_field, summary_function)
  neighbor_list <- st_touches(coarser_sf, sparse = TRUE)
  field_rules_local      <- field_rules
  summary_field_local    <- summary_field
  summary_function_local <- summary_function
  neighbor_list_local    <- neighbor_list
  reassignment_op        <- get("discrete_reassignment_operation", mode = "function", inherits = TRUE)
  call_summary_fn        <- get("call_summary_function", mode = "function", inherits = TRUE)
  
  run_one_discrete_region <- function(iter) {
    # Run one fixed-count discrete reassignment simulation.

    current_fine   <- finer_sf
    current_coarse <- coarser_sf
    
    sampled_ids <- sample(
      current_fine$ID,
      size    = min(k_regions, nrow(current_fine)),
      replace = FALSE
    )
    
    for (fine_id in sampled_ids) {
      
      fine_row  <- which(current_fine$ID == fine_id)
      coarse_id <- current_fine$coarse_ID[fine_row]
      
      result <- tryCatch(
        reassignment_op(
          sf_fine       = current_fine,
          sf_coarse     = current_coarse,
          i             = fine_row,
          neighbor_list = neighbor_list_local,
          fine_id       = fine_id,
          coarse_id     = coarse_id,
          field_rules   = field_rules_local
        ),
        error = function(e) NULL
      )
      
      if (is.null(result)) next
      
      current_fine   <- result$sf_fine
      current_coarse <- result$sf_coarse
    }

    call_summary_fn(current_coarse, summary_field_local, summary_function_local)
  }
  
  distribution <- unlist(run_monte_carlo(
    X              = seq_len(n_iterations),
    FUN            = run_one_discrete_region,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = summary_function_export_names(summary_function),
    export_env     = environment(summary_function),
    progress_label = "Discrete reassignment sensitivity by number of regions"
  ))
  
  out <- list(
    origin_value = origin_value,
    distribution = data.frame(summary_value = distribution)
  )
  
  # ── Persist to disk (optional) ──────────────────────────────────────────
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
  }
  
  return(out)
}

