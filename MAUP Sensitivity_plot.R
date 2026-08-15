# On the sensitivities to the modifiable areal unit problem
# Shared plotting functions for vector and raster sensitivity results
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# yexiang@nnu.edu.cn
# 2026-08-16

# Please cite the following reference when part or all of the code in this file
# is reused under the license of CC-BY-4.0:
# Ye, X., & Chen, J. (2026). On the sensitivities to the modifiable areal unit
# problem. Big Earth Data, 1–36. https://doi.org/10.1080/20964471.2026.2692263

# This file converts vector- and raster-based sensitivity results to one
# standardized table before plotting. Results may be supplied directly or
# read from .rds files produced by the sensitivity calculation functions.


read_sensitivity_result <- function(result_obj = NULL, rds_path = NULL) {
  # Read one sensitivity result from memory or disk. Exactly one input source
  # must be supplied so the origin of the plotted data remains unambiguous.

  if (is.null(result_obj) == is.null(rds_path)) {
    stop("Supply exactly one of result_obj or rds_path.")
  }

  if (!is.null(result_obj)) {
    return(result_obj)
  }

  if (!is.character(rds_path) || length(rds_path) != 1L || !nzchar(rds_path)) {
    stop("rds_path must be one non-empty file path.")
  }

  if (!file.exists(rds_path)) {
    stop(sprintf("File does not exist: %s", rds_path))
  }

  readRDS(rds_path)
}


is_sensitivity_result <- function(x) {
  # Identify one result object by its distribution and pre-MAUP summary value.

  is.list(x) &&
    "distribution" %in% names(x) &&
    "origin_value" %in% names(x)
}


sensitivity_result_baseline <- function(result_obj) {
  # Extract the summary-function value calculated before the MAUP operation.

  if (!"origin_value" %in% names(result_obj)) {
    stop("Result object must contain 'origin_value'.")
  }

  baseline <- result_obj$origin_value

  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline)) {
    stop("The baseline statistic must be one finite numeric value.")
  }

  if (baseline == 0) {
    stop("The baseline statistic is zero and cannot be used for standardization.")
  }

  as.numeric(baseline)
}


sensitivity_result_value_column <- function(distribution) {
  # Locate the summary-function value calculated after the MAUP operation.

  if (!"summary_value" %in% names(distribution)) {
    stop("distribution must contain a 'summary_value' column.")
  }

  "summary_value"
}


default_sensitivity_panel_label <- function(result_obj, fallback = "Sensitivity") {
  # Derive an informative panel label from result metadata when the caller has
  # not supplied one explicitly.

  if (!is.null(result_obj$alpha_multiplier)) {
    return(paste0("alpha = ", format(result_obj$alpha_multiplier), "s"))
  }

  if (!is.null(result_obj$k_pixels)) {
    return(paste0("k = ", result_obj$k_pixels, " pixels"))
  }

  fallback
}


standardize_one_sensitivity_result <- function(result_obj,
                                                panel_label = NULL,
                                                terminal_k = FALSE) {
  # Convert one vector or raster result into the common plotting schema:
  # panel, k, summary_value, baseline, and standardized ratio.

  if (!is_sensitivity_result(result_obj)) {
    stop("Each result must contain a distribution and a supported baseline field.")
  }

  distribution <- result_obj$distribution
  if (!is.data.frame(distribution)) {
    distribution <- as.data.frame(distribution)
  }

  value_column <- sensitivity_result_value_column(distribution)
  baseline <- sensitivity_result_baseline(result_obj)

  # Failed raster paths remain useful diagnostics in saved results but should
  # not be treated as observations in a sensitivity distribution.
  if ("status" %in% names(distribution)) {
    distribution <- distribution[
      is.na(distribution$status) | distribution$status == "Success",
      ,
      drop = FALSE
    ]
  }

  values <- as.numeric(distribution[[value_column]])
  keep <- is.finite(values)
  distribution <- distribution[keep, , drop = FALSE]
  values <- values[keep]

  if (nrow(distribution) == 0L) {
    stop("No successful finite sensitivity values are available for plotting.")
  }

  k_values <- if ("k" %in% names(distribution)) {
    as.numeric(distribution$k)
  } else {
    rep(NA_real_, nrow(distribution))
  }

  # A collection of separate k-order analyses normally represents one terminal
  # distribution per item. Selecting its final order prevents lower-order steps
  # from being duplicated across the collection.
  if (isTRUE(terminal_k) && any(is.finite(k_values))) {
    terminal_value <- if (!is.null(result_obj$k_order)) {
      as.numeric(result_obj$k_order)[1L]
    } else {
      max(k_values, na.rm = TRUE)
    }
    keep_terminal <- is.finite(k_values) & k_values == terminal_value
    distribution <- distribution[keep_terminal, , drop = FALSE]
    values <- values[keep_terminal]
    k_values <- k_values[keep_terminal]
  }

  if (is.null(panel_label)) {
    panel_label <- default_sensitivity_panel_label(result_obj)
  }

  data.frame(
    panel    = rep(as.character(panel_label), length(values)),
    k        = k_values,
    summary_value = values,
    baseline      = rep(baseline, length(values)),
    ratio         = values / baseline,
    stringsAsFactors = FALSE
  )
}


standardize_sensitivity_result <- function(result_obj = NULL,
                                           rds_path = NULL,
                                           panel_labels = NULL,
                                           collection_mode = c("terminal", "all")) {
  # Convert a single result or a named collection of results to the common
  # plotting schema used by every plotting function in this file.

  collection_mode <- match.arg(collection_mode)
  result_obj <- read_sensitivity_result(result_obj = result_obj, rds_path = rds_path)

  if (is_sensitivity_result(result_obj)) {
    standardized <- standardize_one_sensitivity_result(
      result_obj = result_obj,
      panel_label = if (length(panel_labels) > 0L) panel_labels[1L] else NULL,
      terminal_k = FALSE
    )

    # A single high-order result contains one distribution for each k level;
    # use k itself as the panel label so all orders remain visible.
    if (any(is.finite(standardized$k))) {
      standardized$panel <- paste0("k = ", standardized$k)
    }

    standardized$panel <- factor(
      standardized$panel,
      levels = unique(standardized$panel)
    )
    return(standardized)
  }

  if (!is.list(result_obj) || length(result_obj) == 0L ||
      !all(vapply(result_obj, is_sensitivity_result, logical(1L)))) {
    stop("Input must be one sensitivity result or a non-empty collection of results.")
  }

  item_names <- names(result_obj)
  if (is.null(item_names)) {
    item_names <- rep("", length(result_obj))
  }

  if (is.null(panel_labels)) {
    panel_labels <- vapply(seq_along(result_obj), function(i) {
      fallback <- if (nzchar(item_names[i])) item_names[i] else paste0("Result ", i)
      default_sensitivity_panel_label(result_obj[[i]], fallback = fallback)
    }, character(1L))
  }

  if (length(panel_labels) != length(result_obj)) {
    stop("panel_labels must have the same length as the result collection.")
  }

  standardized_parts <- lapply(seq_along(result_obj), function(i) {
    standardize_one_sensitivity_result(
      result_obj = result_obj[[i]],
      panel_label = panel_labels[i],
      terminal_k = identical(collection_mode, "terminal")
    )
  })

  standardized <- do.call(rbind, standardized_parts)
  rownames(standardized) <- NULL
  standardized$panel <- factor(standardized$panel, levels = panel_labels)
  standardized
}


calculate_sensitivity_plot_range <- function(ratio, x_range = NULL) {
  # Use robust tail quantiles by default so isolated extreme simulations do not
  # compress the main distribution. Add padding for legible panel boundaries.

  if (!is.null(x_range)) {
    if (!is.numeric(x_range) || length(x_range) != 2L ||
        !all(is.finite(x_range)) || x_range[1L] >= x_range[2L]) {
      stop("x_range must contain two finite increasing numeric values.")
    }
    return(as.numeric(x_range))
  }

  plot_range <- as.numeric(stats::quantile(
    ratio,
    probs = c(0.001, 0.999),
    na.rm = TRUE,
    names = FALSE
  ))

  if (!all(is.finite(plot_range)) || plot_range[1L] == plot_range[2L]) {
    center <- mean(ratio, na.rm = TRUE)
    padding <- max(abs(center) * 0.01, 0.01)
    return(center + c(-padding, padding))
  }

  padding <- diff(plot_range) * 0.08
  plot_range + c(-padding, padding)
}


build_sensitivity_distribution_plot <- function(plot_data,
                                                x_range = NULL,
                                                bins = 50L,
                                                density_adjust = 2,
                                                base_family = "Times New Roman",
                                                ncol = 1L,
                                                show_stats = TRUE,
                                                vline_at = 1,
                                                vline_color = "red",
                                                x_title = "Summary function (standardized)") {
  # Build the shared histogram-density figure after all result formats have
  # been normalized. This is the only function that defines plot geometry and
  # styling, preventing vector and raster plotting code from diverging.

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for sensitivity plotting.")
  }

  if (!is.numeric(bins) || length(bins) != 1L || bins < 1L) {
    stop("bins must be one positive integer.")
  }

  x_range <- calculate_sensitivity_plot_range(plot_data$ratio, x_range)
  binwidth <- diff(x_range) / as.integer(bins)

  panel_sd <- stats::aggregate(
    ratio ~ panel,
    data = plot_data,
    FUN = stats::sd,
    na.rm = TRUE
  )
  panel_sd$label <- if (isTRUE(show_stats)) {
    paste0("SD = ", sprintf("%.4f", panel_sd$ratio))
  } else {
    ""
  }

  # Histogram counts and probability density use different units. Multiplying
  # density by the largest panel size and bin width places both on a comparable
  # visual scale while retaining density units on the secondary axis.
  panel_sizes <- table(plot_data$panel)
  density_scale <- max(as.numeric(panel_sizes)) * binwidth
  if (!is.finite(density_scale) || density_scale <= 0) {
    density_scale <- 1
  }

  plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ratio)) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      fill = "lightblue",
      alpha = 0.7,
      color = "white",
      linewidth = 0.3,
      ggplot2::aes(y = ggplot2::after_stat(count))
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density) * density_scale),
      color = "blue",
      linewidth = 0.9,
      adjust = density_adjust,
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = vline_at,
      color = vline_color,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = panel_sd,
      ggplot2::aes(x = x_range[1L], y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.25,
      family = base_family,
      fontface = "bold",
      size = 3
    ) +
    ggplot2::coord_cartesian(xlim = x_range, clip = "off") +
    ggplot2::facet_wrap(~ panel, ncol = as.integer(ncol), scales = "fixed") +
    ggplot2::scale_y_continuous(
      name = "Frequency",
      sec.axis = ggplot2::sec_axis(~ . / density_scale, name = "Density"),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(x = x_title) +
    ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "plain"),
      axis.title.x = ggplot2::element_text(size = 9, margin = ggplot2::margin(t = 7)),
      axis.title.y.left = ggplot2::element_text(size = 9, margin = ggplot2::margin(r = 6)),
      axis.title.y.right = ggplot2::element_text(size = 9, margin = ggplot2::margin(l = 6)),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.9, "lines"),
      plot.margin = ggplot2::margin(t = 6, r = 14, b = 12, l = 10)
    )

  plot_object
}


save_sensitivity_plot <- function(plot_object,
                                  filename,
                                  width,
                                  height,
                                  dpi = 300,
                                  device = NULL) {
  # Save a completed sensitivity figure and create its parent directory when
  # necessary. If device is omitted, infer it from the filename extension.

  if (is.null(filename)) {
    return(invisible(NULL))
  }

  if (!is.character(filename) || length(filename) != 1L || !nzchar(filename)) {
    stop("filename must be one non-empty file path.")
  }

  if (is.null(device)) {
    device <- tolower(tools::file_ext(filename))
    if (!nzchar(device)) {
      device <- "png"
      filename <- paste0(filename, ".png")
    }
  }

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "cm",
    dpi = dpi,
    device = device
  )
  message(sprintf("Plot saved to: %s", normalizePath(filename, mustWork = FALSE)))
  invisible(filename)
}


plot_sensitivity_distribution <- function(
    result_obj = NULL,
    rds_path = NULL,
    filename = NULL,
    panel_labels = NULL,
    collection_mode = c("terminal", "all"),
    x_range = NULL,
    bins = 50L,
    density_adjust = 2,
    base_family = "Times New Roman",
    ncol = 1L,
    width = 8.5,
    height = NULL,
    dpi = 300,
    device = NULL,
    show_stats = TRUE,
    vline_at = 1,
    vline_color = "red",
    x_title = "Summary function (standardized)"
) {
  # Plot any supported vector or raster sensitivity result. A single result,
  # result collection, or corresponding .rds file is accepted through the same
  # public interface.

  collection_mode <- match.arg(collection_mode)
  plot_data <- standardize_sensitivity_result(
    result_obj = result_obj,
    rds_path = rds_path,
    panel_labels = panel_labels,
    collection_mode = collection_mode
  )

  n_panels <- nlevels(plot_data$panel)
  ncol <- max(1L, min(as.integer(ncol), n_panels))
  if (is.null(height)) {
    height <- 5.25 * ceiling(n_panels / ncol)
  }

  plot_object <- build_sensitivity_distribution_plot(
    plot_data = plot_data,
    x_range = x_range,
    bins = bins,
    density_adjust = density_adjust,
    base_family = base_family,
    ncol = ncol,
    show_stats = show_stats,
    vline_at = vline_at,
    vline_color = vline_color,
    x_title = x_title
  )

  save_sensitivity_plot(
    plot_object = plot_object,
    filename = filename,
    width = width,
    height = height,
    dpi = dpi,
    device = device
  )

  invisible(plot_object)
}


plot_sensitivity_grid <- function(..., ncol = 1L) {
  # Plot a multi-result or multi-order sensitivity figure. This convenience
  # interface uses the same data normalization and plot engine as
  # plot_sensitivity_distribution().

  plot_sensitivity_distribution(..., ncol = ncol)
}
