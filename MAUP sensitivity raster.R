# On the sensitivities to the modifiable areal unit problem
# Raster sensitivity functions
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# 2026-07-01

# This file stores raster-based MAUP sensitivity functions.
# It is designed to be sourced from the raster walk-through script.

default_parallel_workers_raster <- function(reserve_cores = 1L) {
  # Choose a conservative number of parallel workers.
  n_cores <- parallel::detectCores(logical = TRUE)

  if (is.na(n_cores) || n_cores < 2L) {
    return(1L)
  }

  max(1L, n_cores - as.integer(reserve_cores))
}


run_monte_carlo_raster <- function(X,
                                   FUN,
                                   parallel = FALSE,
                                   n_workers = default_parallel_workers_raster(),
                                   random_seed = NULL,
                                   export_names = character(0),
                                   export_env = .GlobalEnv,
                                   progress_label = "Raster Monte Carlo") {
  # Run Monte Carlo simulations either serially or with a temporary PSOCK cluster.
  # The cluster is stopped automatically when the function exits.

  if (!isTRUE(parallel)) {
    if (!is.null(random_seed)) {
      set.seed(random_seed)
    }
    return(lapply(X, FUN))
  }

  n_workers <- min(as.integer(n_workers), length(X))
  n_workers <- max(1L, n_workers)

  message(sprintf(
    "%s: using %d parallel workers for %d iterations",
    progress_label,
    n_workers,
    length(X)
  ))

  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(random_seed)) {
    parallel::clusterSetRNGStream(cl, random_seed)
  }

  export_names <- export_names[nzchar(export_names)]
  if (length(export_names) > 0L) {
    export_values <- stats::setNames(
      lapply(export_names, function(name) get(name, envir = export_env, inherits = TRUE)),
      export_names
    )
    parallel::clusterCall(cl, function(values) {
      list2env(values, envir = .GlobalEnv)
      NULL
    }, export_values)
  }

  batch_size <- max(n_workers, ceiling(length(X) / 10))
  batch_ids <- split(seq_along(X), ceiling(seq_along(X) / batch_size))
  results <- vector("list", length(X))

  for (batch in batch_ids) {
    batch_results <- parallel::parLapply(cl, X[batch], FUN)
    results[batch] <- batch_results

    message(sprintf(
      "%s: completed %d / %d iterations",
      progress_label,
      max(batch),
      length(X)
    ))
  }

  results
}


summary_90th_percentile_raster <- function(x) {
  # Compute the 90th percentile of a numeric vector while ignoring missing values.
  as.numeric(stats::quantile(x, probs = 0.9, na.rm = TRUE))
}


call_raster_summary_function <- function(summary_function, zone_values) {
  # Apply a user-supplied summary function to county-level raster values.
  value <- summary_function(zone_values)

  if (!is.numeric(value) || length(value) != 1L) {
    stop("summary_function must return one numeric value.")
  }

  as.numeric(value)
}


build_raster_neighbor_index <- function(r_zones, valid_idx, directions = 4L) {
  # Build a compact neighbor list for valid raster cells.
  # Cell IDs are mapped from full raster cell numbers to dense valid-cell indices.

  adj_matrix <- terra::adjacent(
    r_zones,
    cells      = valid_idx,
    directions = directions,
    pairs      = TRUE
  )

  valid_adj <- adj_matrix[adj_matrix[, 2L] %in% valid_idx, , drop = FALSE]

  map_idx <- integer(terra::ncell(r_zones))
  map_idx[valid_idx] <- seq_along(valid_idx)

  adj_pairs <- cbind(
    map_idx[valid_adj[, 1L]],
    map_idx[valid_adj[, 2L]]
  )

  adj_pairs <- adj_pairs[adj_pairs[, 1L] > 0L & adj_pairs[, 2L] > 0L, , drop = FALSE]
  adj_pairs <- unique(rbind(adj_pairs, adj_pairs[, c(2L, 1L), drop = FALSE]))

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
  # Check whether removing one pixel preserves 4-neighbor connectivity
  # within the pixel's original zone.

  nbrs_same <- neighbor_list[[pixel_idx]]
  nbrs_same <- nbrs_same[zone_vals[nbrs_same] == zone_id]
  n_target_nbrs <- length(nbrs_same)

  if (n_target_nbrs <= 1L) {
    return(TRUE)
  }

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

  bfs_env$visited[touched[1L:n_touched]] <- FALSE

  is_connected
}


reassign_one_raster_pixel <- function(zone_vals,
                                      pixel_idx,
                                      zone_counts,
                                      neighbor_list,
                                      bfs_env) {
  # Try to move one boundary pixel from its current zone to an adjacent zone.

  own_zone <- zone_vals[pixel_idx]
  own_key <- as.character(own_zone)
  nbrs <- neighbor_list[[pixel_idx]]

  target_cands <- unique(zone_vals[nbrs[zone_vals[nbrs] != own_zone]])

  if (length(target_cands) == 0L) {
    return(list(result_status = "NoTarget"))
  }

  if (zone_counts[own_key] <= 1L) {
    return(list(result_status = "ConnectivityFail"))
  }

  if (!check_zone_connected_after_removal_raster(
    zone_vals      = zone_vals,
    pixel_idx      = pixel_idx,
    zone_id        = own_zone,
    neighbor_list  = neighbor_list,
    bfs_env        = bfs_env
  )) {
    return(list(result_status = "ConnectivityFail"))
  }

  target_zone <- sample(target_cands, 1L)

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

  n_cells <- length(zone_vals)
  remaining <- as.integer(k_pixels)
  consecutive_failure <- 0L
  total_failure <- 0L
  n_reassigned <- 0L

  is_candidate <- logical(n_cells)
  is_reassigned <- logical(n_cells)

  diff_mask <- zone_vals[adj_pairs[, 1L]] != zone_vals[adj_pairs[, 2L]]
  boundary_cells <- unique(adj_pairs[diff_mask, 1L])
  is_candidate[boundary_cells] <- TRUE

  bfs_env <- new.env(parent = emptyenv())
  bfs_env$visited <- logical(n_cells)
  bfs_env$queue <- integer(n_cells)

  while (consecutive_failure < max_failure && remaining > 0L) {
    avail <- which(is_candidate & !is_reassigned)

    if (length(avail) == 0L) {
      break
    }

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
      own_key <- as.character(zone_vals[sel])
      target_key <- as.character(res$target_zone)

      zone_vals[sel] <- res$target_zone
      zone_counts[own_key] <- zone_counts[own_key] - 1L
      zone_counts[target_key] <- zone_counts[target_key] + 1L

      remaining <- remaining - 1L
      n_reassigned <- n_reassigned + 1L
      is_reassigned[sel] <- TRUE
      consecutive_failure <- 0L

      for (nid in neighbor_list[[sel]]) {
        nn <- neighbor_list[[nid]]
        if (length(nn) > 0L && any(zone_vals[nn] != zone_vals[nid])) {
          is_candidate[nid] <- TRUE
        }
      }
    } else {
      consecutive_failure <- consecutive_failure + 1L
      total_failure <- total_failure + 1L
    }
  }

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

  message("Preparing raster values and adjacency index...")

  data_values <- terra::values(r_data, mat = FALSE)
  zone_values <- terra::values(r_zones, mat = FALSE)

  valid_mask <- !is.na(data_values) & !is.na(zone_values)
  valid_idx <- which(valid_mask)

  if (length(valid_idx) == 0L) {
    stop("No valid overlapping cells were found in r_data and r_zones.")
  }

  zone_vals_init <- zone_values[valid_idx]
  data_vals_init <- data_values[valid_idx]
  zone_counts_init <- table(zone_vals_init)

  neighbor_index <- build_raster_neighbor_index(
    r_zones     = r_zones,
    valid_idx   = valid_idx,
    directions  = directions
  )

  if (nrow(neighbor_index$adj_pairs) == 0L) {
    stop("No adjacent valid raster cells were found.")
  }

  zone_means_base <- as.numeric(tapply(
    data_vals_init,
    zone_vals_init,
    mean,
    na.rm = TRUE
  ))

  variate_before <- call_raster_summary_function(
    summary_function = summary_function,
    zone_values      = zone_means_base
  )

  message(sprintf("Baseline summary value: %.6f", variate_before))
  message(sprintf("Running raster pixel reassignment engine with k_pixels = %d", k_pixels))

  adj_pairs_local <- neighbor_index$adj_pairs
  neighbor_list_local <- neighbor_index$neighbor_list
  zone_vals_init_local <- zone_vals_init
  data_vals_init_local <- data_vals_init
  zone_counts_init_local <- zone_counts_init
  summary_function_local <- summary_function
  max_failure_local <- max_failure
  k_pixels_local <- k_pixels

  run_one_iteration <- function(iter) {
    zone_vals <- zone_vals_init_local
    zone_counts <- zone_counts_init_local

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
        message("Iteration error: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(res_run) || res_run$result_status != "Success") {
      return(data.frame(
        iteration           = iter,
        summary_value       = NA_real_,
        n_reassigned        = if (is.null(res_run)) NA_integer_ else res_run$n_reassigned,
        total_failure       = if (is.null(res_run)) NA_integer_ else res_run$total_failure,
        consecutive_failure = if (is.null(res_run)) NA_integer_ else res_run$consecutive_failure,
        status              = if (is.null(res_run)) "Error" else "Failed",
        stringsAsFactors    = FALSE
      ))
    }

    zone_means_new <- as.numeric(tapply(
      data_vals_init_local,
      res_run$zone_vals,
      mean,
      na.rm = TRUE
    ))

    summary_value <- call_raster_summary_function(
      summary_function = summary_function_local,
      zone_values      = zone_means_new
    )

    data.frame(
      iteration           = iter,
      summary_value       = summary_value,
      n_reassigned        = res_run$n_reassigned,
      total_failure       = res_run$total_failure,
      consecutive_failure = res_run$consecutive_failure,
      status              = "Success",
      stringsAsFactors    = FALSE
    )
  }

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
      "call_raster_summary_function"
    ),
    export_env     = environment(),
    progress_label = "Raster pixel reassignment"
  )

  distribution <- do.call(rbind, result_list)
  rownames(distribution) <- NULL

  out <- list(
    distribution      = distribution,
    results           = distribution,
    variate_before    = variate_before,
    baseline_value    = variate_before,
    k_pixels          = k_pixels,
    n_iterations      = n_iterations,
    directions        = directions,
    n_valid_cells     = length(valid_idx),
    summary_function  = deparse(substitute(summary_function)),
    sensitivity_type  = "pixel_reassignment_exact_count_raster",
    random_seed       = random_seed,
    parallel          = parallel,
    n_workers         = if (isTRUE(parallel)) n_workers else 1L,
    call              = match.call()
  )

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


discrete_reassignment_sensitivity_raster <- function(...) {
  # Backward-compatible alias for old scripts. New walk-through code should
  # use pixel_reassignment_sensitivity_raster(), which implements alpha * s.
  warning(
    paste(
      "discrete_reassignment_sensitivity_raster() is retained only as a",
      "backward-compatible exact-count helper. For pixel reassignment",
      "sensitivity with alpha = 0.2s, 0.5s, or s, use",
      "pixel_reassignment_sensitivity_raster()."
    ),
    call. = FALSE
  )
  exact_pixel_reassignment_sensitivity_raster(...)
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
  # Run raster pixel reassignment sensitivity.
  # alpha_multiplier is multiplied by s, where s is the number of pixels
  # covering the average county size.

  if (!inherits(r_zones, "SpatRaster")) {
    stop("r_zones must be a terra SpatRaster.")
  }

  zone_values <- terra::values(r_zones, mat = FALSE)
  zone_values <- zone_values[!is.na(zone_values)]

  if (length(zone_values) == 0L) {
    stop("No valid zone pixels were found in r_zones.")
  }

  county_pixel_counts <- table(zone_values)
  s_pixels <- mean(as.numeric(county_pixel_counts), na.rm = TRUE)
  alpha_pixels <- max(1L, as.integer(round(alpha_multiplier * s_pixels)))

  message(sprintf(
    "Pixel reassignment alpha = %.2f s | s = %.2f pixels | target pixels = %d",
    alpha_multiplier,
    s_pixels,
    alpha_pixels
  ))

  out <- exact_pixel_reassignment_sensitivity_raster(
    r_data           = r_data,
    r_zones          = r_zones,
    k_pixels         = alpha_pixels,
    summary_function = summary_function,
    n_iterations     = n_iterations,
    random_seed      = random_seed,
    max_failure      = max_failure,
    directions       = directions,
    parallel         = parallel,
    n_workers        = n_workers,
    save_path        = save_path
  )

  out$s_pixels <- s_pixels
  out$alpha_multiplier <- alpha_multiplier
  out$alpha_pixels <- alpha_pixels
  out$sensitivity_type <- "pixel_reassignment_raster"

  if (!is.null(save_path)) {
    saveRDS(out, file = save_path)
  }

  out
}


plot_pixel_reassignment_exact_count_raster <- function(
    result_obj,
    filename = NULL,
    bins = 50,
    density_adjust = 2,
    base_family = "Times New Roman",
    width = 8.5,
    height = 5.25,
    dpi = 300,
    device = "png"
) {
  # Plot the standardized distribution of one raster pixel reassignment result.

  if (!all(c("distribution", "variate_before") %in% names(result_obj))) {
    stop("result_obj must contain 'distribution' and 'variate_before'.")
  }

  plot_data <- result_obj$distribution
  plot_data <- plot_data[plot_data$status == "Success" & is.finite(plot_data$summary_value), , drop = FALSE]

  if (nrow(plot_data) < 2L) {
    stop("At least two successful iterations are required for plotting.")
  }

  plot_data$ratio <- plot_data$summary_value / result_obj$variate_before

  std_val <- stats::sd(plot_data$ratio, na.rm = TRUE)
  x_range <- range(plot_data$ratio, na.rm = TRUE)
  binwidth <- diff(x_range) / bins

  if (!is.finite(binwidth) || binwidth <= 0) {
    binwidth <- 1
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ratio)) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      fill = "lightblue",
      alpha = 0.7,
      color = "white",
      linewidth = 0.3,
      ggplot2::aes(y = ggplot2::after_stat(count))
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count * binwidth)),
      color = "blue",
      linewidth = 0.8,
      adjust = density_adjust
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::labs(
      x = "Summary function (standardized)",
      y = "Frequency",
      subtitle = paste0("SD = ", sprintf("%.4f", std_val))
    ) +
    ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = 9),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 8, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 9, hjust = 0.5, margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(size = 9, margin = ggplot2::margin(r = 6)),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 6, r = 12, b = 12, l = 10)
    )

  if (!is.null(filename)) {
    if (device == "png" && !grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(tools::file_path_sans_ext(filename), ".png")
    }

    if (device == "pdf" && !grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(tools::file_path_sans_ext(filename), ".pdf")
    }

    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      units = "cm",
      dpi = dpi,
      device = device
    )

    message(sprintf("Plot saved to: %s", normalizePath(filename, mustWork = FALSE)))
  }

  invisible(p)
}


plot_pixel_reassignment_grid_raster <- function(
    result_list,
    filename = NULL,
    panel_labels = NULL,
    x_range = NULL,
    bins = 50,
    density_adjust = 2,
    base_family = "Times New Roman",
    width = 8.5,
    height = 15.75,
    dpi = 300,
    device = "png"
) {
  # Plot alpha-based pixel reassignment sensitivity as a three-panel figure.

  if (!is.list(result_list) || length(result_list) == 0L) {
    stop("result_list must be a non-empty list of pixel reassignment result objects.")
  }

  plot_parts <- vector("list", length(result_list))
  sd_parts <- vector("list", length(result_list))

  if (is.null(panel_labels)) {
    panel_labels <- vapply(result_list, function(result_obj) {
      if (!is.null(result_obj$alpha_multiplier)) {
        paste0("alpha == ", format(result_obj$alpha_multiplier, nsmall = 1), "*s")
      } else if (!is.null(result_obj$alpha_pixels)) {
        paste0("alpha[pixels] == ", result_obj$alpha_pixels)
      } else {
        "alpha"
      }
    }, character(1L))
  }

  if (length(panel_labels) != length(result_list)) {
    stop("panel_labels must have the same length as result_list.")
  }

  for (i in seq_along(result_list)) {
    result_obj <- result_list[[i]]

    if (!all(c("distribution", "variate_before") %in% names(result_obj))) {
      stop("Each result object must contain 'distribution' and 'variate_before'.")
    }

    one_panel <- result_obj$distribution
    one_panel <- one_panel[
      one_panel$status == "Success" & is.finite(one_panel$summary_value),
      ,
      drop = FALSE
    ]

    if (nrow(one_panel) < 2L) {
      stop("Each panel requires at least two successful iterations for plotting.")
    }

    one_panel$ratio <- one_panel$summary_value / result_obj$variate_before
    one_panel$panel <- panel_labels[[i]]
    plot_parts[[i]] <- one_panel

    sd_parts[[i]] <- data.frame(
      panel = panel_labels[[i]],
      sd_label = paste0("SD = ", sprintf("%.4f", stats::sd(one_panel$ratio, na.rm = TRUE))),
      stringsAsFactors = FALSE
    )
  }

  plot_data <- do.call(rbind, plot_parts)
  plot_data$panel <- factor(plot_data$panel, levels = panel_labels)

  sd_data <- do.call(rbind, sd_parts)
  sd_data$panel <- factor(sd_data$panel, levels = panel_labels)

  if (is.null(x_range)) {
    x_range <- as.numeric(stats::quantile(plot_data$ratio, c(0.001, 0.999), na.rm = TRUE))
    x_pad <- diff(x_range) * 0.08
    x_range <- x_range + c(-x_pad, x_pad)
  }

  if (!all(is.finite(x_range)) || x_range[1L] == x_range[2L]) {
    x_range <- range(plot_data$ratio, na.rm = TRUE)
  }

  binwidth <- diff(x_range) / bins
  if (!is.finite(binwidth) || binwidth <= 0) {
    binwidth <- 1
  }

  n_ref <- max(table(plot_data$panel))
  density_axis_scale <- as.numeric(n_ref) * binwidth

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ratio)) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      fill = "lightblue",
      alpha = 0.7,
      color = "white",
      linewidth = 0.25,
      ggplot2::aes(y = ggplot2::after_stat(count))
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count * binwidth)),
      color = "blue",
      linewidth = 0.9,
      adjust = density_adjust
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = sd_data,
      ggplot2::aes(x = x_range[1L], y = Inf, label = sd_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.25,
      family = base_family,
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::coord_cartesian(xlim = x_range, clip = "off") +
    ggplot2::facet_wrap(
      ~ panel,
      ncol = 1,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_y_continuous(
      name = "Frequency",
      sec.axis = ggplot2::sec_axis(~ . / density_axis_scale, name = "Density"),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(x = "Summary function (standardized)") +
    ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "plain"),
      axis.title.x = ggplot2::element_text(size = 9, margin = ggplot2::margin(t = 7)),
      axis.title.y.left = ggplot2::element_text(size = 9, margin = ggplot2::margin(r = 6)),
      axis.title.y.right = ggplot2::element_text(size = 9, margin = ggplot2::margin(l = 6)),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(0.9, "lines"),
      plot.margin = ggplot2::margin(t = 6, r = 14, b = 12, l = 10)
    )

  if (!is.null(filename)) {
    if (device == "png" && !grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(tools::file_path_sans_ext(filename), ".png")
    }

    if (device == "pdf" && !grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(tools::file_path_sans_ext(filename), ".pdf")
    }

    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      units = "cm",
      dpi = dpi,
      device = device
    )

    message(sprintf("Plot saved to: %s", normalizePath(filename, mustWork = FALSE)))
  }

  invisible(p)
}


plot_pixel_reassignment_raster <- function(...) {
  # Alias for plotting pixel reassignment sensitivity results.
  plot_pixel_reassignment_exact_count_raster(...)
}


plot_discrete_reassignment_raster <- function(...) {
  # Backward-compatible alias for old scripts.
  warning(
    "Use plot_pixel_reassignment_raster() for the alpha-based pixel reassignment workflow.",
    call. = FALSE
  )
  plot_pixel_reassignment_raster(...)
}


calculate_raster_zone_values <- function(data_values,
                                         zone_values,
                                         area_values = NULL,
                                         merge_method = c("mean", "area_weighted")) {
  # Calculate one value for each raster zone from vectorized raster cells.
  # For PRISM precipitation, area_weighted returns area-weighted mean
  # precipitation for each zone.

  merge_method <- match.arg(merge_method)

  if (merge_method == "mean") {
    zone_means <- tapply(data_values, zone_values, mean, na.rm = TRUE)
    return(as.numeric(zone_means))
  }

  if (is.null(area_values)) {
    stop("area_values must be supplied when merge_method = 'area_weighted'.")
  }

  weighted_sum <- tapply(data_values * area_values, zone_values, sum, na.rm = TRUE)
  area_sum <- tapply(area_values, zone_values, sum, na.rm = TRUE)

  as.numeric(weighted_sum / area_sum)
}


build_raster_zone_adjacency <- function(zone_values,
                                        adj_pairs) {
  # Build a zone-level adjacency list from cell-level adjacency pairs.

  zone_pair_matrix <- cbind(
    zone_values[adj_pairs[, 1L]],
    zone_values[adj_pairs[, 2L]]
  )

  zone_pair_matrix <- zone_pair_matrix[
    !is.na(zone_pair_matrix[, 1L]) &
      !is.na(zone_pair_matrix[, 2L]) &
      zone_pair_matrix[, 1L] != zone_pair_matrix[, 2L],
    ,
    drop = FALSE
  ]

  zone_ids <- sort(unique(zone_values[!is.na(zone_values)]))
  adjacency_list <- vector("list", length(zone_ids))
  names(adjacency_list) <- as.character(zone_ids)

  if (nrow(zone_pair_matrix) == 0L) {
    return(adjacency_list)
  }

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


update_raster_zone_adjacency_after_merge <- function(adj_list,
                                                     target_id,
                                                     source_id) {
  # Update a zone-level adjacency list after source_id is merged into target_id.

  source_key <- as.character(source_id)
  target_key <- as.character(target_id)

  source_neighbors <- adj_list[[source_key]]
  target_neighbors <- adj_list[[target_key]]

  new_neighbors <- unique(c(target_neighbors, source_neighbors))
  new_neighbors <- new_neighbors[new_neighbors != target_id & new_neighbors != source_id]

  adj_list[[target_key]] <- new_neighbors
  adj_list[[source_key]] <- NULL

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
    directions       = 8L,
    parallel         = FALSE,
    n_workers        = default_parallel_workers_raster(),
    save_path        = NULL
) {
  # Run raster-based first- through k-order merging sensitivity.
  # Each Monte Carlo iteration performs sequential adjacent zone merges and
  # records the summary statistic after each merge order.

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

  message("Preparing raster merging sensitivity inputs...")

  data_values_all <- terra::values(r_data, mat = FALSE)
  zone_values_all <- terra::values(r_zones, mat = FALSE)

  valid_mask <- !is.na(data_values_all) & !is.na(zone_values_all)
  valid_idx <- which(valid_mask)

  if (length(valid_idx) == 0L) {
    stop("No valid overlapping cells were found in r_data and r_zones.")
  }

  data_values <- data_values_all[valid_idx]
  zone_values_init <- zone_values_all[valid_idx]

  area_values <- NULL
  if (merge_method == "area_weighted") {
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

  initial_adj_list <- build_raster_zone_adjacency(
    zone_values = zone_values_init,
    adj_pairs   = neighbor_index$adj_pairs
  )

  baseline_values <- calculate_raster_zone_values(
    data_values  = data_values,
    zone_values  = zone_values_init,
    area_values  = area_values,
    merge_method = merge_method
  )

  variate_before <- call_raster_summary_function(
    summary_function = summary_function,
    zone_values      = baseline_values
  )

  message(sprintf("Baseline summary value: %.6f", variate_before))
  message(sprintf("Running raster merging sensitivity for k = 1 to %d", k_order))

  data_values_local <- data_values
  zone_values_init_local <- zone_values_init
  area_values_local <- area_values
  initial_adj_list_local <- initial_adj_list
  summary_function_local <- summary_function
  merge_method_local <- merge_method
  k_order_local <- k_order

  run_one_merge_path <- function(iter) {
    current_zone_values <- zone_values_init_local
    current_adj_list <- initial_adj_list_local

    rows <- vector("list", k_order_local)

    for (merge_level in seq_len(k_order_local)) {
      zone_ids <- as.numeric(names(current_adj_list))
      valid_zones <- zone_ids[vapply(current_adj_list, length, integer(1L)) > 0L]

      if (length(valid_zones) == 0L) {
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

      target_id <- sample(valid_zones, 1L)
      source_id <- sample(current_adj_list[[as.character(target_id)]], 1L)

      current_zone_values[current_zone_values == source_id] <- target_id

      current_adj_list <- update_raster_zone_adjacency_after_merge(
        adj_list   = current_adj_list,
        target_id  = target_id,
        source_id  = source_id
      )

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

  result_list <- run_monte_carlo_raster(
    X              = seq_len(n_iterations),
    FUN            = run_one_merge_path,
    parallel       = parallel,
    n_workers      = n_workers,
    random_seed    = random_seed,
    export_names   = c(
      "calculate_raster_zone_values",
      "update_raster_zone_adjacency_after_merge",
      "call_raster_summary_function"
    ),
    export_env     = environment(),
    progress_label = "Raster merging sensitivity"
  )

  distribution <- do.call(rbind, result_list)
  rownames(distribution) <- NULL

  out <- list(
    distribution      = distribution,
    results           = distribution,
    variate_before    = variate_before,
    baseline_value    = variate_before,
    merge_method      = merge_method,
    merge_method_used = merge_method,
    k_order           = k_order,
    n_iterations      = n_iterations,
    directions        = directions,
    n_valid_cells     = length(valid_idx),
    summary_function  = deparse(substitute(summary_function)),
    sensitivity_type  = "merging_raster",
    random_seed       = random_seed,
    parallel          = parallel,
    n_workers         = if (isTRUE(parallel)) n_workers else 1L,
    call              = match.call()
  )

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


plot_merging_sensitivity_raster <- function(
    result_obj,
    filename = NULL,
    x_range = NULL,
    bins = 50,
    density_adjust = 2,
    base_family = "Times New Roman",
    width = 8.5,
    height = 15.75,
    dpi = 300,
    device = "png"
) {
  # Plot first- through k-order raster merging sensitivity distributions.

  if (!all(c("distribution", "variate_before", "k_order") %in% names(result_obj))) {
    stop("result_obj must contain 'distribution', 'variate_before', and 'k_order'.")
  }

  plot_data <- result_obj$distribution
  plot_data <- plot_data[
    plot_data$status == "Success" & is.finite(plot_data$summary_value),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) < 2L) {
    stop("At least two successful records are required for plotting.")
  }

  plot_data$ratio <- plot_data$summary_value / result_obj$variate_before
  plot_data$k <- factor(plot_data$k, levels = sort(unique(plot_data$k)))

  sd_data <- stats::aggregate(
    ratio ~ k,
    data = plot_data,
    FUN = stats::sd,
    na.rm = TRUE
  )
  sd_data$sd_label <- paste0("SD = ", sprintf("%.4f", sd_data$ratio))

  if (is.null(x_range)) {
    x_range <- as.numeric(stats::quantile(plot_data$ratio, c(0.001, 0.999), na.rm = TRUE))
    x_pad <- diff(x_range) * 0.08
    x_range <- x_range + c(-x_pad, x_pad)
  }

  if (!all(is.finite(x_range)) || x_range[1L] == x_range[2L]) {
    x_range <- range(plot_data$ratio, na.rm = TRUE)
  }

  binwidth <- diff(x_range) / bins
  if (!is.finite(binwidth) || binwidth <= 0) {
    binwidth <- 1
  }

  n_ref <- max(table(plot_data$k))
  density_axis_scale <- as.numeric(n_ref) * binwidth

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ratio)) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      fill = "lightblue",
      alpha = 0.7,
      color = "white",
      linewidth = 0.3,
      ggplot2::aes(y = ggplot2::after_stat(count))
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count * binwidth)),
      color = "blue",
      linewidth = 0.9,
      adjust = density_adjust
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = sd_data,
      ggplot2::aes(x = x_range[1L], y = Inf, label = sd_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.25,
      family = base_family,
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::coord_cartesian(xlim = x_range, clip = "off") +
    ggplot2::facet_wrap(
      ~ k,
      ncol = 1,
      labeller = ggplot2::labeller(k = function(x) paste0("k = ", x))
    ) +
    ggplot2::scale_y_continuous(
      name = "Frequency",
      sec.axis = ggplot2::sec_axis(~ . / density_axis_scale, name = "Density"),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      x = "Summary function (standardized)"
    ) +
    ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "plain"),
      axis.title.x = ggplot2::element_text(size = 9, hjust = 0.5, margin = ggplot2::margin(t = 7)),
      axis.title.y.left = ggplot2::element_text(size = 9, margin = ggplot2::margin(r = 6)),
      axis.title.y.right = ggplot2::element_text(size = 9, margin = ggplot2::margin(l = 6)),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(0.9, "lines"),
      plot.margin = ggplot2::margin(t = 6, r = 18, b = 14, l = 10)
    )

  if (!is.null(filename)) {
    if (device == "png" && !grepl("\\.png$", filename, ignore.case = TRUE)) {
      filename <- paste0(tools::file_path_sans_ext(filename), ".png")
    }

    if (device == "pdf" && !grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(tools::file_path_sans_ext(filename), ".pdf")
    }

    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      units = "cm",
      dpi = dpi,
      device = device
    )

    message(sprintf("Plot saved to: %s", normalizePath(filename, mustWork = FALSE)))
  }

  invisible(p)
}


split_raster_zone_values_once <- function(zone_values,
                                          target_id,
                                          neighbor_list,
                                          next_zone_id,
                                          n_attempts = 10L) {
  # Split one raster zone into two connected parts by two-seed region growing.
  # The function works on compact valid-cell indices rather than editing a
  # SpatRaster object, which keeps Monte Carlo iterations lightweight.

  zone_cells <- which(zone_values == target_id)

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

  for (attempt in seq_len(as.integer(n_attempts))) {
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

    unassigned <- which(assigned == 0L)
    for (node in unassigned) {
      nbr_labels <- assigned[adj_list[[node]]]
      nbr_labels <- nbr_labels[nbr_labels != 0L]

      assigned[node] <- if (length(nbr_labels) > 0L) {
        as.integer(names(which.max(table(nbr_labels))))
      } else {
        1L
      }
    }

    region1_cells <- zone_cells[assigned == 1L]
    region2_cells <- zone_cells[assigned == 2L]

    if (length(region1_cells) == 0L || length(region2_cells) == 0L) {
      next
    }

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
    directions       = 8L,
    parallel         = FALSE,
    n_workers        = default_parallel_workers_raster(),
    save_path        = NULL
) {
  # Run raster-based first- through k-order splitting sensitivity.
  # Each Monte Carlo iteration sequentially splits eligible raster zones and
  # records the summary statistic after each split order.

  split_method <- match.arg(split_method)

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

  message("Preparing raster splitting sensitivity inputs...")

  data_values_all <- terra::values(r_data, mat = FALSE)
  zone_values_all <- terra::values(r_zones, mat = FALSE)

  valid_mask <- !is.na(data_values_all) & !is.na(zone_values_all)
  valid_idx <- which(valid_mask)

  if (length(valid_idx) == 0L) {
    stop("No valid overlapping cells were found in r_data and r_zones.")
  }

  data_values <- data_values_all[valid_idx]
  zone_values_init <- zone_values_all[valid_idx]

  area_values <- NULL
  if (split_method == "area_weighted") {
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

  baseline_values <- calculate_raster_zone_values(
    data_values  = data_values,
    zone_values  = zone_values_init,
    area_values  = area_values,
    merge_method = split_method
  )

  variate_before <- call_raster_summary_function(
    summary_function = summary_function,
    zone_values      = baseline_values
  )

  initial_cell_counts <- table(zone_values_init)
  next_zone_id_init <- max(zone_values_init, na.rm = TRUE) + 1L

  message(sprintf("Baseline summary value: %.6f", variate_before))
  message(sprintf("Running raster splitting sensitivity for k = 1 to %d", k_order))

  data_values_local <- data_values
  zone_values_init_local <- zone_values_init
  area_values_local <- area_values
  neighbor_list_local <- neighbor_index$neighbor_list
  summary_function_local <- summary_function
  split_method_local <- split_method
  k_order_local <- k_order
  n_attempts_local <- as.integer(n_attempts)
  initial_cell_counts_local <- initial_cell_counts
  next_zone_id_init_local <- next_zone_id_init

  run_one_split_path <- function(iter) {
    current_zone_values <- zone_values_init_local
    cell_counts <- initial_cell_counts_local
    next_zone_id <- next_zone_id_init_local

    rows <- vector("list", k_order_local)

    for (split_level in seq_len(k_order_local)) {
      eligible_ids <- as.numeric(names(cell_counts[cell_counts >= 2L]))

      if (length(eligible_ids) == 0L) {
        rows[[split_level]] <- data.frame(
          iteration        = iter,
          k                = split_level,
          split_id         = NA_real_,
          new_id           = NA_real_,
          summary_value    = NA_real_,
          status           = "Failed",
          stringsAsFactors = FALSE
        )
        next
      }

      target_id <- sample(eligible_ids, 1L)

      split_result <- split_raster_zone_values_once(
        zone_values   = current_zone_values,
        target_id     = target_id,
        neighbor_list = neighbor_list_local,
        next_zone_id  = next_zone_id,
        n_attempts    = n_attempts_local
      )

      if (is.null(split_result)) {
        rows[[split_level]] <- data.frame(
          iteration        = iter,
          k                = split_level,
          split_id         = target_id,
          new_id           = NA_real_,
          summary_value    = NA_real_,
          status           = "Failed",
          stringsAsFactors = FALSE
        )
        next
      }

      current_zone_values <- split_result$zone_values

      cell_counts[as.character(target_id)] <- length(split_result$region1_cells)
      cell_counts[as.character(split_result$new_zone_id)] <- length(split_result$region2_cells)
      next_zone_id <- next_zone_id + 1L

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
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, rows)
  }

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

  out <- list(
    distribution      = distribution,
    results           = distribution,
    variate_before    = variate_before,
    baseline_value    = variate_before,
    split_method      = split_method,
    split_method_used = split_method,
    k_order           = k_order,
    n_iterations      = n_iterations,
    n_attempts        = as.integer(n_attempts),
    directions        = directions,
    n_valid_cells     = length(valid_idx),
    summary_function  = deparse(substitute(summary_function)),
    sensitivity_type  = "splitting_raster",
    random_seed       = random_seed,
    parallel          = parallel,
    n_workers         = if (isTRUE(parallel)) n_workers else 1L,
    call              = match.call()
  )

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


plot_splitting_sensitivity_raster <- function(
    result_obj,
    filename = NULL,
    x_range = NULL,
    bins = 50,
    density_adjust = 2,
    base_family = "Times New Roman",
    width = 8.5,
    height = 15.75,
    dpi = 300,
    device = "png"
) {
  # Plot first- through k-order raster splitting sensitivity distributions.

  plot_merging_sensitivity_raster(
    result_obj     = result_obj,
    filename       = filename,
    x_range        = x_range,
    bins           = bins,
    density_adjust = density_adjust,
    base_family    = base_family,
    width          = width,
    height         = height,
    dpi            = dpi,
    device         = device
  )
}
