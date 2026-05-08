# On the sensitivities to the modifiable areal unit problem
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# yexiang@nnu.edu.cn
# 2026-04-30

library(geojsonio)
library(readxl)
library(tidyverse)
library(sf)
library(lwgeom)
library(ggplot2)
library(spdep)

merging_operation <- function(sf, i, j, field_rules) {
  
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

merging_sensitivity <- function(sf,
                                field_rules,
                                summary_field,
                                summary_function,
                                k_order              = 1,
                                exhaustive_threshold = 1e3,
                                n_iterations         = 100,
                                random_seed          = NULL,
                                n_simulations        = 200,
                                save_path            = NULL) {
  
  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.function(summary_function))         stop("summary_function must be a function")
  if (!summary_field %in% colnames(sf))       stop(paste0("Field '", summary_field, "' does not exist in the sf object"))
  if (!is.numeric(sf[[summary_field]]))       stop(paste0("Field '", summary_field, "' must be numeric"))
  if (k_order < 1)                            stop("k_order must be >= 1")
  
  if (!is.null(random_seed)) set.seed(random_seed)
  
  # ── Compute baseline statistic ────────────────────────────────────────────
  origin_value <- tryCatch({
    summary_function(sf[[summary_field]])
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
    
    message(sprintf("Starting exhaustive enumeration of all %d-order merge paths...", k_order))
    
    all_results <- data.frame(
      path_id = integer(),
      k       = integer(),
      i       = character(),
      j       = character(),
      new     = integer(),
      value   = numeric(),
      stringsAsFactors = FALSE
    )
    
    path_counter <- 0L
    
    explore_merge_paths <- function(current_sf, current_k, path_history) {
      
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
            value   = step$value,
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
            summary_function(temp_sf[[summary_field]])
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
              value = current_value
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
    
    message(sprintf("Starting random sampling of %d-order merging, %d iterations...", k_order, n_iterations))
    
    results <- data.frame(
      path_id = integer(),
      k       = integer(),
      i       = character(),
      j       = character(),
      new     = integer(),
      value   = numeric(),
      stringsAsFactors = FALSE
    )
    
    skipped <- 0L
    
    for (iter in seq_len(n_iterations)) {
      temp_sf    <- sf
      iter_rows  <- vector("list", k_order)
      iter_valid <- TRUE
      
      for (merge_level in seq_len(k_order)) {
        
        G          <- poly2nb(temp_sf, queen = FALSE)
        candidates <- which(sapply(G, length) > 0)
        
        if (length(candidates) == 0) {
          warning(sprintf("Iteration %d has no mergeable regions at step %d; discarding", iter, merge_level))
          iter_valid <- FALSE
          break
        }
        
        # Randomly pick a pair of adjacent units
        i_index <- sample(candidates, 1)
        j_index <- G[[i_index]][sample(length(G[[i_index]]), 1)]
        id_i    <- temp_sf$ID[i_index]
        id_j    <- temp_sf$ID[j_index]
        
        temp_sf <- tryCatch(
          merging_operation(temp_sf, id_i, id_j, field_rules),
          error = function(e) {
            warning(sprintf("Iteration %d step %d merge failed; discarding this iteration: %s",
                            iter, merge_level, e$message))
            NULL
          }
        )
        
        if (is.null(temp_sf)) {
          iter_valid <- FALSE
          break
        }
        
        new_id <- max(temp_sf$ID, na.rm = TRUE)
        
        current_value <- tryCatch({
          summary_function(temp_sf[[summary_field]])
        }, error = function(e) {
          warning(paste0("summary_function failed (iteration ", iter,
                         ", k=", merge_level, "): ", e$message))
          NA
        })
        
        iter_rows[[merge_level]] <- data.frame(
          path_id = iter,
          k       = merge_level,
          i       = as.character(id_i),
          j       = as.character(id_j),
          new     = new_id,
          value   = current_value,
          stringsAsFactors = FALSE
        )
      }
      
      if (iter_valid) {
        results <- rbind(results, do.call(rbind, iter_rows))
      } else {
        skipped <- skipped + 1L
      }
      
      if (iter %% 10 == 0) {
        message(sprintf("Completed %d/%d iterations (%d discarded)", iter, n_iterations, skipped))
      }
    }
    
    message(sprintf("Random sampling finished: %d successful, %d discarded", n_iterations - skipped, skipped))
    return(results)
  }
  
  # ── Run the analysis ──────────────────────────────────────────────────────
  if (use_exhaustive) {
    message(sprintf("Combination count (%.2e) is below the threshold (%.2e); using exhaustive method",
                    estimated_combinations, exhaustive_threshold))
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

splitting_operation <- function(
    sf,
    target_id,
    field_rules = NULL,
    max_iter    = 100
) {
  if (!inherits(sf, "sf"))       stop("Input must be an sf object")
  if (!target_id %in% sf$ID)     stop(paste0("Target ID does not exist: ", target_id))
  
  target_index <- which(sf$ID == target_id)
  target       <- sf[target_index, ] |> st_make_valid()
  
  attempts <- 0
  parts    <- NULL
  
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
      st_split(st_geometry(target), split_line) |>
        st_collection_extract("POLYGON") |>
        st_make_valid()
    }, error = function(e) NULL)
  }
  
  # Compute area ratios, used by the area_weighted rule
  areas      <- st_area(parts)
  ratios     <- as.numeric(areas / sum(areas))
  base_attrs <- sf[target_index, , drop = FALSE] |> st_drop_geometry()
  
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
                                  save_path    = NULL) {
  
  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.function(summary_function))   stop("summary_function must be a function")
  if (!summary_field %in% colnames(sf)) stop(paste0("Field '", summary_field, "' does not exist in the sf object"))
  if (!is.numeric(sf[[summary_field]])) stop(paste0("Field '", summary_field, "' must be numeric"))
  if (k_order < 1)                      stop("k_order must be >= 1")
  
  if (!is.null(random_seed)) set.seed(random_seed)
  
  # ── Compute baseline statistic ────────────────────────────────────────────
  # Record the statistic on the original data as a reference baseline before
  # any split operation is performed
  origin_value <- tryCatch({
    summary_function(sf[[summary_field]])
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
    value   = numeric(),   # statistic after the split
    stringsAsFactors = FALSE
  )
  
  successful   <- 0L
  attempts     <- 0L
  # Cap the total number of attempts to prevent infinite loops when split
  # conditions are very restrictive
  max_attempts <- max(n_iterations * 10L, n_iterations + 1000L)
  
  while (successful < n_iterations && attempts < max_attempts) {
    attempts   <- attempts + 1L
    temp_sf    <- sf          # restart from the original data each path
    iter_rows  <- vector("list", k_order)
    iter_valid <- TRUE
    
    for (split_level in seq_len(k_order)) {
      
      # Safety check: if sf becomes unexpectedly empty, abort this path
      if (nrow(temp_sf) == 0) {
        iter_valid <- FALSE
        break
      }
      
      # Randomly pick a region to split
      i_index    <- sample(seq_len(nrow(temp_sf)), 1)
      id_i       <- temp_sf$ID[i_index]
      ids_before <- temp_sf$ID  # snapshot of IDs to detect newly created ones
      
      # Perform the split; abort this path if it fails (e.g. the region is
      # not splittable)
      temp_sf <- tryCatch(
        splitting_operation(temp_sf, id_i, field_rules),
        error = function(e) NULL
      )
      
      if (is.null(temp_sf)) {
        iter_valid <- FALSE
        break
      }
      
      # Identify the two newly created sub-region IDs by diffing the ID set
      new_ids <- sort(setdiff(temp_sf$ID, ids_before))
      new1    <- if (length(new_ids) >= 1) new_ids[1] else NA_integer_
      new2    <- if (length(new_ids) >= 2) new_ids[2] else NA_integer_
      
      # Compute the statistic after this step; on failure, record NA and
      # continue rather than aborting the whole path
      current_value <- tryCatch({
        summary_function(temp_sf[[summary_field]])
      }, error = function(e) {
        warning(paste0("summary_function failed (attempt ", attempts,
                       ", k=", split_level, "): ", e$message))
        NA
      })
      
      iter_rows[[split_level]] <- data.frame(
        path_id = successful + 1L,
        k       = split_level,
        i       = as.character(id_i),
        new1    = new1,
        new2    = new2,
        value   = current_value,
        stringsAsFactors = FALSE
      )
    }
    
    # Only paths that completed all k_order steps are kept
    if (iter_valid) {
      successful <- successful + 1L
      results    <- rbind(results, do.call(rbind, iter_rows))
      
      if (successful %% 10 == 0) {
        message(sprintf("Successfully completed %d/%d (total attempts: %d)",
                        successful, n_iterations, attempts))
      }
    }
  }
  
  # If the attempt cap is reached without meeting the target, raise a clear
  # warning rather than silently returning a deficient result
  if (successful < n_iterations) {
    warning(sprintf(
      "Reached the maximum attempt cap (%d); only %d/%d paths succeeded. Consider relaxing the split conditions or lowering k_order",
      max_attempts, successful, n_iterations
    ))
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

generate_smooth_guided_path <- function(start, end) {
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

continuous_reassignment_operation <- function(
    sf,
    alpha,
    field_rules,
    summary_field,
    summary_function,
    max_iter  = 1000,
    tol_ratio = 0.1
) {
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(sf, "sf"))              stop("Input must be an sf object")
  if (!summary_field %in% names(sf))    stop(sprintf("Field '%s' does not exist", summary_field))
  if (!is.numeric(sf[[summary_field]])) stop("summary_field must be numeric")
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
    summary_function(sf[[summary_field]]),
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
    
    # -- Pick a source polygon from the unblocked candidate pool --------
    candidates <- setdiff(current_sf$ID, source_blocked)
    if (length(candidates) == 0L) {
      termination_reason <- "candidate set is empty"
      break
    }
    
    source_id  <- sample(candidates, 1)
    source_row <- which(current_sf$ID == source_id)
    
    # -- Identify topological neighbors of the source -------------------
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
    
    # -- Split the source polygon into sub-parts ------------------------
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
    
    # -- Try to merge one sub-part into a neighboring polygon -----------
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
        summary_function(merged_result[[summary_field]]),
        error = function(e) NA_real_
      )
      
      # -- Log this successful step -----------------------------------
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
      
      # -- Update state for the next iteration ------------------------
      current_stat    <- new_stat
      remaining_alpha <- remaining_alpha - part_area
      # Block the target so it cannot be picked as a source in the next
      # iteration. Without this guard the algorithm could immediately
      # split off the same area it just merged in, producing a no-op
      # cycle and wasting iterations.
      source_blocked  <- union(source_blocked, target_id)
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
    summary_function(current_sf[[summary_field]]),
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


continuous_sensitivity <- function(
    sf,
    field_rules,
    summary_field,
    summary_function,
    alpha,
    n_iterations = 100,
    max_iter     = 1000,
    tol_ratio    = 0.1,
    random_seed  = NULL,
    save_path    = NULL
) {
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(sf, "sf"))
    stop("Input must be an sf object")
  if (!is.function(summary_function))
    stop("summary_function must be a function")
  if (!summary_field %in% colnames(sf))
    stop(paste0("Field '", summary_field, "' does not exist, please check column names"))
  if (!is.numeric(sf[[summary_field]]))
    stop(paste0("Field '", summary_field, "' must be numeric"))
  
  if (!is.null(random_seed)) set.seed(random_seed)
  
  # ── Compute baseline statistic ──────────────────────────────────────────
  # Recorded once on the untouched data so every iteration's result can be
  # compared against a single, stable reference value.
  origin_value <- tryCatch(
    summary_function(sf[[summary_field]]),
    error = function(e) stop(paste0("summary_function failed: ", e$message))
  )
  message(sprintf("Baseline statistic: %.6f", origin_value))
  
  # ── Initialize result container ─────────────────────────────────────────
  results <- data.frame(
    iteration       = integer(0),
    value           = numeric(0),
    total_failure   = integer(0),
    remaining_alpha = numeric(0),
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Starting continuous sensitivity analysis: %d iterations", n_iterations))
  
  # ── Main loop ───────────────────────────────────────────────────────────
  # Each iteration starts from the original `sf` (not the previous result),
  # so the n_iterations runs are mutually independent Monte-Carlo samples.
  for (run in seq_len(n_iterations)) {
    
    res <- tryCatch(
      continuous_reassignment_operation(
        sf               = sf,
        alpha            = alpha,
        field_rules      = field_rules,
        summary_field    = summary_field,
        summary_function = summary_function,
        max_iter         = max_iter,
        tol_ratio        = tol_ratio
      ),
      error = function(e) {
        warning(sprintf("Run %d failed: %s", run, e$message))
        NULL
      }
    )
    
    # Failed runs (NULL) are skipped silently; only successful runs
    # contribute a row to the distribution.
    if (!is.null(res)) {
      results <- rbind(results, data.frame(
        iteration       = run,
        value           = res$stat_after,
        total_failure   = res$total_failure,
        remaining_alpha = res$remaining_alpha,
        stringsAsFactors = FALSE
      ))
    }
    
    # Coarse-grained progress reporting: avoids flooding the console on
    # large n_iterations while still giving the user a sense of progress.
    if (run %% 10 == 0 || run == n_iterations) {
      message(sprintf("Progress: %d/%d completed (%d successful)",
                      run, n_iterations, nrow(results)))
    }
  }
  
  message(sprintf("Done: %d/%d successful runs collected",
                  nrow(results), n_iterations))
  
  out <- list(
    origin_value = origin_value,
    distribution = results,
    field_used   = summary_field
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

discrete_area_sensitivity <- function(
    finer_sf,
    coarser_sf,
    field_rules,
    n_iterations,
    alpha,
    summary_field,
    summary_function,
    tol_ratio  = 0.1,
    keep_maps  = FALSE,
    save_path  = NULL        # If non-NULL, results are saved as a .rds file
) {
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(finer_sf,   "sf")) stop("finer_sf must be an sf object")
  if (!inherits(coarser_sf, "sf")) stop("coarser_sf must be an sf object")
  if (!summary_field %in% names(coarser_sf))
    stop(sprintf("Field '%s' does not exist in coarser_sf", summary_field))
  if (!is.numeric(coarser_sf[[summary_field]]))
    stop("summary_field is not numeric")
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
  origin_value  <- summary_function(coarser_sf[[summary_field]])
  neighbor_list <- st_touches(coarser_sf, sparse = TRUE)
  fine_areas    <- as.numeric(st_area(finer_sf))
  
  results       <- vector("list", n_iterations)
  success_count <- 0L
  
  # ── Main Monte-Carlo loop ───────────────────────────────────────────────
  for (iter in seq_len(n_iterations)) {
    
    budget <- runif(1, min = alpha - delta_alpha,
                    max = alpha + delta_alpha)
    
    current_fine   <- finer_sf
    current_coarse <- coarser_sf
    remaining      <- budget
    available_ids  <- current_fine$ID
    any_reassigned <- FALSE
    
    while (remaining > 0 && length(available_ids) > 0L) {
      
      fine_id  <- sample(available_ids, 1)
      fine_row <- which(current_fine$ID == fine_id)
      area_i   <- fine_areas[fine_row]
      
      if (area_i > remaining) {
        available_ids <- setdiff(available_ids, fine_id)
        next
      }
      
      coarse_id <- current_fine$coarse_ID[fine_row]
      
      result <- tryCatch(
        discrete_reassignment_operation(
          sf_fine       = current_fine,
          sf_coarse     = current_coarse,
          i             = fine_row,
          neighbor_list = neighbor_list,
          fine_id       = fine_id,
          coarse_id     = coarse_id,
          field_rules   = field_rules
        ),
        error = function(e) NULL
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
    
    if (any_reassigned) {
      success_count   <- success_count + 1L
      results[[iter]] <- list(
        stat_after = summary_function(current_coarse[[summary_field]]),
        updated_sf = if (keep_maps) current_coarse else NULL
      )
    }
  }
  
  # ── Assemble distribution & deviation summary ───────────────────────────
  valid        <- Filter(Negate(is.null), results)
  distribution <- data.frame(value = sapply(valid, `[[`, "stat_after"))
  devs         <- distribution$value / origin_value
  
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


discrete_region_sensitivity <- function(
    finer_sf,
    coarser_sf,
    field_rules,
    k_regions,
    n_iterations,
    summary_function,
    summary_field,
    save_path = NULL          # If non-NULL, results are saved as a .rds file
) {
  # ── Input validation ────────────────────────────────────────────────────
  if (!inherits(finer_sf,   "sf")) stop("finer_sf must be an sf object")
  if (!inherits(coarser_sf, "sf")) stop("coarser_sf must be an sf object")
  if (!summary_field %in% names(coarser_sf))
    stop(sprintf("Field '%s' does not exist in coarser_sf", summary_field))
  if (!is.numeric(coarser_sf[[summary_field]]))
    stop("summary_field is not numeric")
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
  origin_value  <- summary_function(coarser_sf[[summary_field]])
  neighbor_list <- st_touches(coarser_sf, sparse = TRUE)
  
  distribution <- numeric(n_iterations)
  
  # ── Main Monte-Carlo loop ───────────────────────────────────────────────
  for (iter in seq_len(n_iterations)) {
    
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
        discrete_reassignment_operation(
          sf_fine       = current_fine,
          sf_coarse     = current_coarse,
          i             = fine_row,
          neighbor_list = neighbor_list,
          fine_id       = fine_id,
          coarse_id     = coarse_id,
          field_rules   = field_rules
        ),
        error = function(e) NULL
      )
      
      if (is.null(result)) next
      
      current_fine   <- result$sf_fine
      current_coarse <- result$sf_coarse
    }
    
    distribution[iter] <- summary_function(current_coarse[[summary_field]])
  }
  
  out <- list(
    origin_value = origin_value,
    distribution = data.frame(value = distribution)
  )
  
  # ── Persist to disk (optional) ──────────────────────────────────────────
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file = save_path)
    message(sprintf("Results saved to: %s", normalizePath(save_path, mustWork = FALSE)))
  }
  
  return(out)
}

plot_and_save <- function(
    result_obj      = NULL,    # in-memory result list, OR
    rds_path        = NULL,    # path to a .rds file produced by *_sensitivity()
    filename_prefix = NULL,    # File name prefix; _k{k} suffix is automatically appended
    x_range         = NULL,
    bins            = 50,
    density_adjust  = 2,
    base_family     = "Times New Roman",
    width           = 7,
    height          = 5.25,
    dpi             = 300,
    show_stats      = TRUE,
    vline_color     = "red",
    vline_at        = 1,
    theme_style     = "minimal",
    device          = "png"
) {
  # ── 0. Parse input source (in-memory object OR .rds file, choose one) ─────
  if (is.null(result_obj) && is.null(rds_path))
    stop("Either result_obj or rds_path must be provided.")
  if (!is.null(result_obj) && !is.null(rds_path))
    stop("Only one of result_obj or rds_path can be provided.")
  
  if (!is.null(rds_path)) {
    if (!is.character(rds_path) || length(rds_path) != 1L)
      stop("rds_path must be a single string path.")
    if (!file.exists(rds_path))
      stop(sprintf("File does not exist: %s", rds_path))
    result_obj <- readRDS(rds_path)
    message(sprintf("✓ Read successfully: %s", normalizePath(rds_path, mustWork = FALSE)))
    
    # If filename_prefix is not specified, use the rds filename (without extension) 
    # as the default output prefix, so the plots will be automatically saved in the 
    # same directory as the data for easy management.
    if (is.null(filename_prefix)) {
      filename_prefix <- tools::file_path_sans_ext(rds_path)
    }
  }
  
  # ── 1. Data validation and extraction ─────────────────────────────────────
  if (!is.list(result_obj))
    stop("result_obj must be a list (or a list read from an .rds file).")
  if (!"origin_value"  %in% names(result_obj)) stop("result_obj is missing 'origin_value'.")
  if (!"distribution"  %in% names(result_obj)) stop("result_obj is missing 'distribution'.")
  if (!"value"         %in% colnames(result_obj$distribution)) stop("distribution is missing the 'value' column.")
  
  origin_value <- result_obj$origin_value
  dist_df      <- result_obj$distribution
  
  has_k       <- "k_order" %in% names(result_obj) && "k" %in% colnames(dist_df)
  k_order_max <- if (has_k) result_obj$k_order else 1L
  if (!has_k) dist_df$k <- 1L
  
  # ── 2. Calculate the ratio relative to the baseline ───────────────────────
  dist_df <- dist_df %>%
    mutate(ratio = value / origin_value)
  
  # ── 3. Unify coordinate range & fix binwidth ──────────────────────────────
  if (is.null(x_range)) {
    x_min   <- quantile(dist_df$ratio, 0.001, na.rm = TRUE)
    x_max   <- quantile(dist_df$ratio, 0.999, na.rm = TRUE)
    x_range <- c(x_min, x_max)
  }
  fixed_binwidth <- (x_range[2] - x_range[1]) / bins
  
  # ── 4. Theme setup ────────────────────────────────────────────────────────
  base_theme <- switch(
    theme_style,
    "minimal" = theme_minimal(base_family = base_family),
    "bw"      = theme_bw(base_family = base_family),
    "classic" = theme_classic(base_family = base_family),
    theme_minimal(base_family = base_family)
  )
  
  # ── 5. Single plot factory function ───────────────────────────────────────
  make_panel <- function(df, k_level) {
    
    if (nrow(df) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste0("k = ", k_level, "\n(No data)"),
                   size = 5, family = base_family) +
          theme_void(base_family = base_family)
      )
    }
    
    std_val <- sd(df$ratio, na.rm = TRUE)
    
    hist_data <- ggplot_build(
      ggplot(df, aes(x = ratio)) + geom_histogram(binwidth = fixed_binwidth)
    )$data
    
    dens_data <- ggplot_build(
      ggplot(df, aes(x = ratio)) + geom_density(adjust = density_adjust)
    )$data
    
    scale_factor <- max(hist_data$count, na.rm = TRUE) /
      max(dens_data$density, na.rm = TRUE)
    
    stats_label <- if (show_stats) paste0("SD = ", sprintf("%.4f", std_val)) else NULL
    
    ggplot(df, aes(x = ratio)) +
      
      geom_histogram(
        binwidth  = fixed_binwidth,
        fill      = "lightblue",
        alpha     = 0.7,
        color     = "white",
        linewidth = 0.3,
        aes(y = after_stat(count))
      ) +
      
      geom_density(
        aes(y = after_stat(density) * scale_factor),
        color     = "blue",
        linewidth = 1.2,
        adjust    = density_adjust
      ) +
      
      { if (!is.null(vline_at))
        geom_vline(xintercept = vline_at, color = vline_color,
                   linetype = "dashed", linewidth = 1)
      } +
      
      coord_cartesian(xlim = x_range) +
      
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      
      scale_y_continuous(
        name     = "Frequency",
        labels   = function(x) as.integer(x),
        sec.axis = sec_axis(
          trans = ~ . / scale_factor,
          name  = "Density"
        )
      ) +
      
      labs(
        x        = "Summary function (standardized)",
        subtitle = stats_label
      ) +
      
      base_theme +
      theme(
        text             = element_text(family = base_family, size = 14),
        plot.subtitle    = element_text(hjust = 0, size = 9, face = "bold"),
        axis.title       = element_text(size = 10),
        axis.text        = element_text(size = 8),
        panel.grid.minor = element_blank()
      )
  }
  
  # ── 6. Plot by k and save independently ───────────────────────────────────
  plots <- vector("list", k_order_max)
  
  for (k in seq_len(k_order_max)) {
    df_k       <- dist_df[dist_df$k == k, , drop = FALSE]
    plots[[k]] <- make_panel(df_k, k_level = k)
    
    if (!is.null(filename_prefix)) {
      fname <- if (has_k) {
        paste0(tools::file_path_sans_ext(filename_prefix), "_k", k, ".", device)
      } else {
        paste0(tools::file_path_sans_ext(filename_prefix), ".", device)
      }
      
      # Automatically create parent directories: when filename_prefix is derived 
      # from rds_path (which might be in a nested directory), prevent ggsave 
      # from failing due to missing directories.
      dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
      
      ggsave(
        filename = fname,
        plot     = plots[[k]],
        width    = width,
        height   = height,
        units    = "cm",
        dpi      = dpi,
        device   = device,
        type     = "cairo"
      )
      message(sprintf("✓ Saved: %s", fname))
    }
  }
  
  invisible(if (k_order_max == 1) plots else plots)
}