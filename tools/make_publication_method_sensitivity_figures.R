#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages before running this script: ", paste(missing_packages, collapse = ", "))
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- list(
    results_dir = NULL,
    output_dir = NULL,
    spectrum_rank = 10L,
    dpi = 600L
  )
  
  positional <- character(0)
  for (arg in args) {
    if (startsWith(arg, "--")) {
      parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
      key <- gsub("-", "_", parts[[1L]], fixed = TRUE)
      value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
      if (!key %in% names(cfg)) {
        stop("Unknown argument: ", arg)
      }
      cfg[[key]] <- value
    } else {
      positional <- c(positional, arg)
    }
  }
  
  if (length(positional) >= 1L && is.null(cfg$results_dir)) cfg$results_dir <- positional[[1L]]
  if (length(positional) >= 2L && is.null(cfg$output_dir)) cfg$output_dir <- positional[[2L]]
  if (length(positional) > 2L) stop("Too many positional arguments.")
  
  cfg$spectrum_rank <- as.integer(cfg$spectrum_rank)
  cfg$dpi <- as.integer(cfg$dpi)
  cfg
}

infer_results_dir <- function() {
  env_dir <- Sys.getenv("METHOD_SENSITIVITY_OUT_DIR", unset = NA_character_)
  candidates <- unique(stats::na.omit(c(
    env_dir,
    file.path(getwd(), "out", "method_sensitivity_16pf"),
    file.path(getwd(), "..", "out", "method_sensitivity_16pf")
  )))
  
  ok <- candidates[file.exists(file.path(candidates, "method_sensitivity_summary.csv"))]
  if (length(ok)) return(normalizePath(ok[[1L]], mustWork = TRUE))
  
  stop(
    "Cannot infer the method-sensitivity output directory. ",
    "Run from the study root or pass --results-dir=/path/to/out/method_sensitivity_16pf."
  )
}

normalise_dir <- function(path, must_work = TRUE) {
  normalizePath(path.expand(path), mustWork = must_work)
}

read_result <- function(results_dir, filename) {
  path <- file.path(results_dir, filename)
  if (!file.exists(path)) stop("Missing required cached result: ", path)
  suppressMessages(readr::read_csv2(path, show_col_types = FALSE, progress = FALSE))
}

read_optional_result <- function(results_dir, filename) {
  path <- file.path(results_dir, filename)
  if (!file.exists(path)) return(NULL)
  suppressMessages(readr::read_csv2(path, show_col_types = FALSE, progress = FALSE))
}

write_csv2_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv2(x, path)
  invisible(path)
}

format_method_label <- function(method, multiline = FALSE) {
  labels <- c(
    pca_uniform = if (multiline) "All-item\nPCA" else "All-item PCA",
    robpca_uniform = if (multiline) "All-item\nrobust PCA" else "All-item robust PCA",
    pca_id_guided = if (multiline) "Selected-item\nPCA" else "Selected-item PCA",
    robpca_id_guided = if (multiline) "Selected-item\nrobust PCA" else "Selected-item robust PCA"
  )
  out <- unname(labels[as.character(method)])
  out[is.na(out)] <- gsub("_", if (multiline) "\n" else " ", as.character(method[is.na(out)]), fixed = TRUE)
  out
}

format_subset_label <- function(subset_mode) {
  dplyr::recode(
    subset_mode,
    unstratified = "Unstratified",
    stratified_by_scale = "Scale-stratified",
    .default = gsub("_", " ", subset_mode, fixed = TRUE)
  )
}

format_metric_label <- function(metric) {
  dplyr::recode(
    metric,
    axis_corr_mean = "Axis correlation",
    item_rmse = "Item RMSE",
    scale_regression_vector_cosine = "Scale regression-vector cosine",
    selection_jaccard = "Selection Jaccard",
    aligned_weight_correlation = "Aligned weight Pearson r",
    .default = gsub("_", " ", metric, fixed = TRUE)
  )
}

format_scale_label <- function(scale_name) {
  gsub("_", " ", as.character(scale_name), fixed = TRUE)
}

method_levels <- c("pca_uniform", "robpca_uniform", "pca_id_guided", "robpca_id_guided")
method_shapes <- c(
  pca_uniform = 16,
  robpca_uniform = 17,
  pca_id_guided = 15,
  robpca_id_guided = 18
)
method_linetypes <- c(
  pca_uniform = "solid",
  robpca_uniform = "dashed",
  pca_id_guided = "dotdash",
  robpca_id_guided = "dotted"
)
method_fills <- c(
  pca_uniform = "white",
  robpca_uniform = "grey82",
  pca_id_guided = "grey55",
  robpca_id_guided = "grey25"
)

theme_pub_bw <- function(base_size = 11,
                         legend_position = "bottom",
                         y_grid = TRUE,
                         x_grid = FALSE) {
  grid_theme <- ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(2.5, "pt"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", colour = "black"),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.margin = grid::unit(c(5.5, 6, 5.5, 6), "pt")
    )
  
  if (isTRUE(y_grid)) {
    grid_theme <- grid_theme +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_line(
          colour = "grey82",
          linetype = "dotted",
          linewidth = 0.25
        )
      )
  }
  if (isTRUE(x_grid)) {
    grid_theme <- grid_theme +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(
          colour = "grey88",
          linetype = "dotted",
          linewidth = 0.25
        )
      )
  }
  
  grid_theme + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

open_png <- function(path, width, height, dpi) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "in", res = dpi, background = "white")
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = dpi, bg = "white")
  }
}

open_pdf <- function(path, width, height) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(path, width = width, height = height, bg = "white")
  } else {
    grDevices::pdf(path, width = width, height = height, useDingbats = FALSE, bg = "white")
  }
}

save_plot_pair <- function(plot, name, output_dir, width, height, dpi) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(output_dir, paste0(name, ".png"))
  pdf_path <- file.path(output_dir, paste0(name, ".pdf"))
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi, device = ragg::agg_png, bg = "white")
  } else {
    ggplot2::ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, device = pdf_device, bg = "white")
  invisible(c(png = png_path, pdf = pdf_path))
}

draw_plot_grid <- function(plots, widths = rep(1, length(plots))) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(plots), widths = grid::unit(widths, "null"))))
  on.exit(grid::popViewport(), add = TRUE)
  for (i in seq_along(plots)) {
    print(
      plots[[i]],
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i)
    )
  }
}

save_plot_grid <- function(plots, name, output_dir, width, height, dpi, widths = rep(1, length(plots))) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(output_dir, paste0(name, ".png"))
  pdf_path <- file.path(output_dir, paste0(name, ".pdf"))
  
  open_png(png_path, width = width, height = height, dpi = dpi)
  draw_plot_grid(plots, widths = widths)
  grDevices::dev.off()
  
  open_pdf(pdf_path, width = width, height = height)
  draw_plot_grid(plots, widths = widths)
  grDevices::dev.off()
  
  invisible(c(png = png_path, pdf = pdf_path))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "", paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

write_markdown_table <- function(x, path) {
  escape_cell <- function(z) gsub("\\|", "\\\\|", as.character(z))
  header <- paste0("| ", paste(escape_cell(names(x)), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) paste0("| ", paste(escape_cell(row), collapse = " | "), " |"))
  writeLines(c(header, divider, rows), path, useBytes = TRUE)
  invisible(path)
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na <- function(x) if (sum(!is.na(x)) <= 1L) NA_real_ else stats::sd(x, na.rm = TRUE)

summarise_metric_by_method <- function(data, value_cols) {
  data |>
    dplyr::select(dplyr::any_of(c("method", value_cols))) |>
    tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
    dplyr::group_by(method, metric) |>
    dplyr::summarise(
      n = sum(!is.na(value)),
      mean = mean_or_na(value),
      sd = sd_or_na(value),
      median = if (all(is.na(value))) NA_real_ else stats::median(value, na.rm = TRUE),
      q05 = if (all(is.na(value))) NA_real_ else stats::quantile(value, 0.05, na.rm = TRUE),
      q95 = if (all(is.na(value))) NA_real_ else stats::quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      method_label = format_method_label(method),
      metric_label = format_metric_label(metric),
      .after = method
    )
}

preferred_id_guided_method <- function(methods) {
  methods <- as.character(methods)
  if ("pca_id_guided" %in% methods) return("pca_id_guided")
  id_methods <- methods[grepl("id_guided", methods, fixed = TRUE)]
  if (length(id_methods)) return(id_methods[[1L]])
  methods[[1L]]
}

build_retained_item_table <- function(selected_item_content_tbl, item_weights_tbl) {
  source_tbl <- if (!is.null(selected_item_content_tbl) && nrow(selected_item_content_tbl)) {
    selected_item_content_tbl
  } else {
    item_weights_tbl |> dplyr::filter(.data$selected)
  }
  
  source_method <- preferred_id_guided_method(unique(source_tbl$method))
  source_tbl |>
    dplyr::filter(.data$method == source_method) |>
    dplyr::arrange(dplyr::desc(.data$weight), .data$item_id) |>
    dplyr::mutate(Rank = dplyr::row_number(), .before = 1L) |>
    dplyr::transmute(
      Rank,
      `Item ID` = .data$item_id,
      `Item text` = .data$item_text,
      Scale = .data$scale_name,
      `Original 16PF factor` = .data$original_16pf_factor,
      Keying = .data$keyed_sign,
      `Final weight` = .data$weight,
      `Active-survivor count` = .data$selection_count,
      `Active-survivor proportion` = .data$selection_prop,
      `Mean final weight` = .data$weight_mean
    )
}

build_readout_delta_table <- function(readout_tbl) {
  readout_tbl |>
    dplyr::select(scale_name, method, cv_r2) |>
    tidyr::pivot_wider(names_from = method, values_from = cv_r2) |>
    dplyr::transmute(
      Scale = format_scale_label(.data$scale_name),
      `All-item PCA R2` = .data$pca_uniform,
      `Selected-item PCA R2` = .data$pca_id_guided,
      `PCA delta` = .data$pca_id_guided - .data$pca_uniform,
      `All-item robust PCA R2` = .data$robpca_uniform,
      `Selected-item robust PCA R2` = .data$robpca_id_guided,
      `Robust PCA delta` = .data$robpca_id_guided - .data$robpca_uniform
    ) |>
    dplyr::arrange(.data$Scale)
}

build_stability_summary_table <- function(full_refit_stability_tbl,
                                          fixed_selection_stability_tbl,
                                          selection_stability_tbl,
                                          stability_metrics,
                                          selection_metrics) {
  dplyr::bind_rows(
    summarise_metric_by_method(fixed_selection_stability_tbl, stability_metrics) |>
      dplyr::mutate(analysis = "Fixed selected set", .before = 1L),
    summarise_metric_by_method(full_refit_stability_tbl, stability_metrics) |>
      dplyr::mutate(analysis = "Split-half re-selection maps", .before = 1L),
    summarise_metric_by_method(selection_stability_tbl, selection_metrics) |>
      dplyr::mutate(analysis = "Split-half re-selection item sets", .before = 1L)
  ) |>
    dplyr::transmute(
      Analysis = .data$analysis,
      Method = .data$method_label,
      Metric = .data$metric_label,
      N = .data$n,
      Mean = .data$mean,
      SD = .data$sd,
      Median = .data$median,
      `5th percentile` = .data$q05,
      `95th percentile` = .data$q95
    )
}

build_distance_geometry_table <- function(distance_preservation_tbl) {
  distance_preservation_tbl |>
    dplyr::mutate(method = factor(.data$method, levels = observed_method_levels)) |>
    dplyr::arrange(.data$method) |>
    dplyr::transmute(
      Method = format_method_label(as.character(.data$method)),
      `Evaluation N` = .data$n_eval,
      `Neighbour count` = .data$k_neighbors,
      Items = .data$n_items_selected,
      `Gower TwoNN ID` = .data$gower_twonn_id,
      `Gower-map Pearson r` = .data$map_distance_pearson,
      `Gower-map Spearman r` = .data$map_distance_spearman,
      `Gower-to-map neighbour overlap` = .data$neighbourhood_overlap_gower_to_map,
      `Map-to-Gower neighbour overlap` = .data$neighbourhood_overlap_map_to_gower,
      `PCoA two-axis positive-eigen share` = .data$pcoa2_positive_eigen_share
    )
}

build_selected_map_plot <- function(pc_scores_tbl,
                                    scale_vectors_tbl,
                                    summary_tbl,
                                    method = "robpca_id_guided",
                                    base_size = 11,
                                    fill_high = "grey42") {
  if (!method %in% unique(pc_scores_tbl$method)) {
    id_methods <- unique(pc_scores_tbl$method[grepl("id_guided", pc_scores_tbl$method, fixed = TRUE)])
    method <- if (length(id_methods)) id_methods[[1L]] else unique(pc_scores_tbl$method)[[1L]]
  }
  method_id <- method
  map_scores <- pc_scores_tbl |>
    dplyr::filter(.data$method == method_id) |>
    dplyr::select(participant_id, u1, u2)
  map_summary <- summary_tbl |>
    dplyr::filter(as.character(.data$method) == method_id) |>
    dplyr::slice(1L)
  padded_limits <- function(x, probs = c(0.01, 0.99), pad = 0.15) {
    q <- as.numeric(stats::quantile(x, probs, na.rm = TRUE))
    span <- diff(q)
    if (!is.finite(span) || span <= 0) span <- diff(range(x, na.rm = TRUE))
    if (!is.finite(span) || span <= 0) span <- 1
    q + c(-pad, pad) * span
  }
  x_limits <- padded_limits(map_scores$u1)
  y_limits <- padded_limits(map_scores$u2)
  axis_prefix <- if (nrow(map_summary) && identical(map_summary$decomposition[[1]], "robpca")) "RC" else "PC"
  
  vector_base <- scale_vectors_tbl |>
    dplyr::filter(.data$method == method_id) |>
    dplyr::mutate(vec_len = sqrt(beta1^2 + beta2^2)) |>
    dplyr::filter(is.finite(vec_len), vec_len > 0)
  arrow_scale <- if (nrow(vector_base)) {
    0.42 * min(diff(x_limits), diff(y_limits)) / max(vector_base$vec_len, na.rm = TRUE)
  } else {
    1
  }
  vector_plot_data <- vector_base |>
    dplyr::mutate(
      panel = "B. Scale-score gradients",
      x0 = 0,
      y0 = 0,
      xend = beta1 * arrow_scale,
      yend = beta2 * arrow_scale,
      label = format_scale_label(scale_name)
    )
  map_panel_data <- dplyr::bind_rows(
    map_scores |> dplyr::mutate(panel = "A. Respondent density"),
    map_scores |> dplyr::mutate(panel = "B. Scale-score gradients")
  )
  label_layer <- if (nrow(vector_plot_data) && requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_label_repel(
      data = vector_plot_data,
      ggplot2::aes(x = xend, y = yend, label = label),
      inherit.aes = FALSE,
      size = 2.0,
      min.segment.length = 0,
      box.padding = 0.28,
      point.padding = 0.12,
      force = 4,
      force_pull = 0.15,
      max.overlaps = Inf,
      label.size = 0.12,
      label.padding = grid::unit(0.08, "lines"),
      fill = "white",
      seed = 42
    )
  } else {
    ggplot2::geom_text(
      data = vector_plot_data,
      ggplot2::aes(x = xend, y = yend, label = label),
      inherit.aes = FALSE,
      size = 2.0,
      check_overlap = TRUE
    )
  }
  
  ggplot2::ggplot(map_panel_data, ggplot2::aes(u1, u2)) +
    ggplot2::stat_density_2d(
      ggplot2::aes(fill = ggplot2::after_stat(level)),
      geom = "polygon",
      contour = TRUE,
      bins = 18,
      colour = NA,
      alpha = 0.95
    ) +
    ggplot2::stat_density_2d(
      colour = "grey40",
      linewidth = 0.18,
      bins = 8
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.25, colour = "grey45") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.25, colour = "grey45") +
    ggplot2::geom_segment(
      data = vector_plot_data,
      ggplot2::aes(x = x0, y = y0, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.38,
      arrow = grid::arrow(length = grid::unit(0.07, "in"))
    ) +
    label_layer +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_fill_gradient(low = "grey97", high = fill_high, guide = "none") +
    ggplot2::coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +
    ggplot2::labs(
      x = sprintf("%s1 (%.1f%%)", axis_prefix, 100 * map_summary$pc1_explained_variance_ratio[[1]]),
      y = sprintf("%s2 (%.1f%%)", axis_prefix, 100 * map_summary$pc2_explained_variance_ratio[[1]])
    ) +
    theme_pub_bw(base_size = base_size, y_grid = FALSE, x_grid = FALSE, legend_position = "none")
}

cfg <- parse_args()
results_dir <- normalise_dir(cfg$results_dir %||% infer_results_dir(), must_work = TRUE)
output_dir <- normalise_dir(cfg$output_dir %||% file.path(results_dir, "publication_bw"), must_work = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading cached evaluation outputs from: ", results_dir)
message("Writing publication figures and tables to: ", output_dir)

summary_tbl <- read_result(results_dir, "method_sensitivity_summary.csv")
eigenspectrum_tbl <- read_result(results_dir, "method_sensitivity_eigenspectrum.csv")
baseline_tbl <- read_result(results_dir, "method_sensitivity_matched_random_subset_baseline.csv")
scale_prediction_tbl <- read_result(results_dir, "method_sensitivity_scale_prediction.csv")
leave_one_scale_out_tbl <- read_result(results_dir, "method_sensitivity_leave_one_scale_out_prediction.csv")
full_refit_stability_tbl <- read_result(results_dir, "method_sensitivity_split_half_full_refit_stability.csv")
fixed_selection_stability_tbl <- read_result(results_dir, "method_sensitivity_split_half_fixed_selection_stability.csv")
selection_stability_tbl <- read_result(results_dir, "method_sensitivity_selection_stability.csv")
item_weights_tbl <- read_result(results_dir, "method_sensitivity_item_weights.csv")
item_component_correlations_tbl <- read_result(results_dir, "method_sensitivity_item_component_correlations.csv")
pc_scores_tbl <- read_optional_result(results_dir, "method_sensitivity_pc_scores_2d.csv")
scale_regression_vectors_tbl <- read_optional_result(results_dir, "method_sensitivity_scale_regression_vectors.csv")
distance_preservation_tbl <- read_optional_result(results_dir, "method_sensitivity_distance_preservation.csv")
selected_item_content_tbl <- read_optional_result(results_dir, "method_sensitivity_selected_item_content.csv")

observed_method_levels <- method_levels[method_levels %in% unique(summary_tbl$method)]
if (length(observed_method_levels) != nrow(summary_tbl)) {
  observed_method_levels <- unique(summary_tbl$method)
}
method_label_lookup <- stats::setNames(format_method_label(observed_method_levels, multiline = TRUE), observed_method_levels)
pc12_method_label_lookup <- stats::setNames(
  c(
    pca_uniform = "All-item\nPCA",
    robpca_uniform = "All-item\nrobust PCA",
    pca_id_guided = "Selected\nPCA",
    robpca_id_guided = "Selected\nrobust PCA"
  )[observed_method_levels],
  observed_method_levels
)

summary_tbl <- summary_tbl |>
  dplyr::mutate(
    method = factor(method, levels = observed_method_levels),
    method_label = format_method_label(as.character(method)),
    method_label_multiline = format_method_label(as.character(method), multiline = TRUE)
  )

table1_raw <- summary_tbl |>
  dplyr::transmute(
    method = as.character(method),
    method_label,
    weighting_mode,
    decomposition,
    n_items_selected,
    pc1_explained_variance_ratio,
    pc2_explained_variance_ratio,
    pc12_explained_variance_ratio,
    pc5_explained_variance_ratio,
    pc10_explained_variance_ratio,
    cv_r2_scale_prediction_mean,
    cv_r2_leave_one_scale_out_mean,
    fixed_selection_axis_corr_mean,
    fixed_selection_scale_regression_vector_cosine_mean
  )

table1_pub <- table1_raw |>
  dplyr::transmute(
    Method = method_label,
    Items = n_items_selected,
    `PC1+2` = fmt_pct(pc12_explained_variance_ratio),
    `PC1-5` = fmt_pct(pc5_explained_variance_ratio),
    `PC1-10` = fmt_pct(pc10_explained_variance_ratio),
    `Direct R2` = fmt_num(cv_r2_scale_prediction_mean),
    `LOSO R2` = fmt_num(cv_r2_leave_one_scale_out_mean),
    `Axis r` = fmt_num(fixed_selection_axis_corr_mean),
    `Scale cos.` = fmt_num(fixed_selection_scale_regression_vector_cosine_mean)
  )

write_csv2_safe(table1_raw, file.path(output_dir, "TABLE_1_analytic_solution_summary_raw.csv"))
write_csv2_safe(table1_pub, file.path(output_dir, "TABLE_1_analytic_solution_summary.csv"))
write_markdown_table(table1_pub, file.path(output_dir, "TABLE_1_analytic_solution_summary.md"))

if (!is.null(pc_scores_tbl) && !is.null(scale_regression_vectors_tbl)) {
  fig1_map <- build_selected_map_plot(
    pc_scores_tbl = pc_scores_tbl,
    scale_vectors_tbl = scale_regression_vectors_tbl,
    summary_tbl = summary_tbl,
    method = "robpca_id_guided",
    base_size = 11
  )
  save_plot_pair(
    fig1_map,
    "FIG_3_selected_item_respondent_map_bw",
    output_dir = output_dir,
    width = 7.4,
    height = 4.2,
    dpi = cfg$dpi
  )
}

scree_data <- eigenspectrum_tbl |>
  dplyr::filter(component <= cfg$spectrum_rank) |>
  dplyr::mutate(method = factor(method, levels = observed_method_levels))

fig1_scree <- ggplot2::ggplot(
  scree_data,
  ggplot2::aes(
    x = component,
    y = explained_variance_ratio,
    group = method,
    linetype = method,
    shape = method
  )
) +
  ggplot2::geom_line(colour = "black", linewidth = 0.55) +
  ggplot2::geom_point(colour = "black", fill = "white", size = 1.8, stroke = 0.35) +
  ggplot2::scale_linetype_manual(values = method_linetypes[observed_method_levels], labels = method_label_lookup) +
  ggplot2::scale_shape_manual(values = method_shapes[observed_method_levels], labels = method_label_lookup) +
  ggplot2::scale_x_continuous(breaks = seq_len(cfg$spectrum_rank)) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(tag = "A", x = "Component", y = "Explained variance") +
  theme_pub_bw(base_size = 11, y_grid = TRUE, legend_position = "bottom") +
  ggplot2::guides(
    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  )

pc12_data <- summary_tbl |>
  dplyr::select(method, method_label_multiline, pc1_explained_variance_ratio, pc2_explained_variance_ratio, pc12_explained_variance_ratio) |>
  dplyr::mutate(method = factor(as.character(method), levels = observed_method_levels)) |>
  tidyr::pivot_longer(
    cols = c(pc1_explained_variance_ratio, pc2_explained_variance_ratio),
    names_to = "component",
    values_to = "explained_variance_ratio"
  ) |>
  dplyr::mutate(
    component = factor(
      component,
      levels = c("pc1_explained_variance_ratio", "pc2_explained_variance_ratio"),
      labels = c("PC1", "PC2")
    )
  )

pc12_labels <- summary_tbl |>
  dplyr::mutate(method = factor(as.character(method), levels = observed_method_levels))

fig1_pc12 <- ggplot2::ggplot(
  pc12_data,
  ggplot2::aes(x = method, y = explained_variance_ratio, fill = component)
) +
  ggplot2::geom_col(width = 0.72, colour = "black", linewidth = 0.35) +
  ggplot2::geom_text(
    data = pc12_labels,
    ggplot2::aes(x = method, y = pc12_explained_variance_ratio, label = scales::percent(pc12_explained_variance_ratio, accuracy = 0.1)),
    inherit.aes = FALSE,
    vjust = -0.45,
    size = 2.5
  ) +
  ggplot2::scale_x_discrete(labels = pc12_method_label_lookup) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.14))) +
  ggplot2::scale_fill_manual(values = c(PC1 = "grey35", PC2 = "grey82"), name = NULL) +
  ggplot2::labs(tag = "B", x = NULL, y = "PC1 + PC2 explained variance") +
  theme_pub_bw(base_size = 11, y_grid = TRUE, legend_position = "bottom") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, lineheight = 0.9))

save_plot_grid(
  list(fig1_scree, fig1_pc12),
  "FIG_1_eigenspectrum_two_component_concentration_bw",
  output_dir = output_dir,
  width = 8.2,
  height = 4.2,
  dpi = cfg$dpi,
  widths = c(1.25, 1)
)

target_methods <- intersect(unique(baseline_tbl$target_method), as.character(summary_tbl$method))
baseline_plot_data <- baseline_tbl |>
  dplyr::filter(target_method %in% target_methods) |>
  dplyr::mutate(
    target_method = factor(target_method, levels = target_methods, labels = format_method_label(target_methods)),
    subset_mode = factor(
      subset_mode,
      levels = c("unstratified", "stratified_by_scale"),
      labels = c("Unstratified", "Scale-stratified")
    )
  )

id_guided_reference <- summary_tbl |>
  dplyr::filter(as.character(method) %in% target_methods) |>
  dplyr::transmute(
    target_method = factor(as.character(method), levels = target_methods, labels = format_method_label(target_methods)),
    pc12_explained_variance_ratio
  )

baseline_tail_labels <- baseline_plot_data |>
  dplyr::group_by(target_method, subset_mode) |>
  dplyr::summarise(
    n_reps = dplyr::n(),
    q95 = stats::quantile(pc12_explained_variance_ratio, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    y_label = q95 + 0.002,
    label = sprintf("n=%d\n95th=%s", n_reps, scales::percent(q95, accuracy = 0.1))
  )

fig2 <- ggplot2::ggplot(
  baseline_plot_data,
  ggplot2::aes(x = subset_mode, y = pc12_explained_variance_ratio, fill = subset_mode)
) +
  ggplot2::geom_violin(trim = FALSE, width = 0.82, colour = "black", linewidth = 0.35) +
  ggplot2::geom_boxplot(width = 0.18, fill = "white", colour = "black", outlier.shape = NA, linewidth = 0.35) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.055, height = 0, seed = 42),
    shape = 21,
    colour = "black",
    fill = "white",
    alpha = 0.22,
    size = 0.85,
    stroke = 0.2
  ) +
  ggplot2::geom_hline(
    data = id_guided_reference,
    ggplot2::aes(yintercept = pc12_explained_variance_ratio),
    linetype = "longdash",
    linewidth = 0.45,
    colour = "black"
  ) +
  ggplot2::geom_point(
    data = id_guided_reference,
    ggplot2::aes(x = 1.5, y = pc12_explained_variance_ratio),
    inherit.aes = FALSE,
    shape = 23,
    fill = "white",
    colour = "black",
    size = 2.4,
    stroke = 0.45
  ) +
  ggplot2::geom_text(
    data = baseline_tail_labels |> dplyr::filter(.data$subset_mode == "Unstratified"),
    ggplot2::aes(x = subset_mode, y = y_label, label = label),
    inherit.aes = FALSE,
    position = ggplot2::position_nudge(x = -0.32),
    hjust = 0.5,
    size = 2.2,
    lineheight = 0.9
  ) +
  ggplot2::geom_text(
    data = baseline_tail_labels |> dplyr::filter(.data$subset_mode == "Scale-stratified"),
    ggplot2::aes(x = subset_mode, y = y_label, label = label),
    inherit.aes = FALSE,
    position = ggplot2::position_nudge(x = 0.32),
    hjust = 0.5,
    size = 2.2,
    lineheight = 0.9
  ) +
  ggplot2::facet_wrap(~target_method, nrow = 1) +
  ggplot2::scale_fill_manual(values = c("Unstratified" = "white", "Scale-stratified" = "grey78"), guide = "none") +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.46)) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
  ggplot2::labs(x = "Matched random-subset sampling scheme", y = "PC1 + PC2 explained variance") +
  theme_pub_bw(base_size = 11, y_grid = TRUE, legend_position = "none") +
  ggplot2::coord_cartesian(clip = "off")

save_plot_pair(
  fig2,
  "FIG_2_matched_random_subset_pc12_distribution_bw",
  output_dir = output_dir,
  width = 7.4,
  height = 4.3,
  dpi = cfg$dpi
)

loso_diff_data <- leave_one_scale_out_tbl |>
  dplyr::mutate(
    decomposition = dplyr::case_when(
      grepl("^robpca", method) ~ "Robust PCA",
      TRUE ~ "PCA"
    ),
    weighting = dplyr::case_when(
      grepl("id_guided", method, fixed = TRUE) ~ "id_guided",
      TRUE ~ "uniform"
    )
  ) |>
  dplyr::select(decomposition, weighting, scale_name, cv_r2) |>
  tidyr::pivot_wider(names_from = weighting, values_from = cv_r2) |>
  dplyr::mutate(delta = id_guided - uniform)

scale_order <- loso_diff_data |>
  dplyr::group_by(scale_name) |>
  dplyr::summarise(mean_delta = mean_or_na(delta), .groups = "drop") |>
  dplyr::arrange(mean_delta) |>
  dplyr::pull(scale_name)

fig3_data <- loso_diff_data |>
  dplyr::mutate(
    decomposition = factor(decomposition, levels = c("PCA", "Robust PCA")),
    scale_name = factor(scale_name, levels = scale_order),
    gain = delta >= 0
  )

fig3 <- ggplot2::ggplot(
  fig3_data,
  ggplot2::aes(x = delta, y = scale_name, shape = decomposition, fill = gain)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.35, colour = "grey35") +
  ggplot2::geom_point(
    colour = "black",
    size = 2.1,
    stroke = 0.45
  ) +
  ggplot2::facet_wrap(~decomposition, nrow = 1) +
  ggplot2::scale_shape_manual(values = c("PCA" = 21, "Robust PCA" = 24), guide = "none") +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "grey45", `FALSE` = "white"), guide = "none") +
  ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
  ggplot2::scale_y_discrete(labels = format_scale_label) +
  ggplot2::labs(x = "Delta LOSO R2 (selected-item minus all-item)", y = NULL) +
  theme_pub_bw(base_size = 11, y_grid = TRUE, x_grid = TRUE, legend_position = "none")

save_plot_pair(
  fig3,
  "FIG_4_leave_one_scale_out_delta_by_scale_bw",
  output_dir = output_dir,
  width = 7.5,
  height = 5.8,
  dpi = cfg$dpi
)

direct_scale_order <- scale_prediction_tbl |>
  dplyr::group_by(scale_name) |>
  dplyr::summarise(mean_cv_r2 = mean_or_na(cv_r2), .groups = "drop") |>
  dplyr::arrange(mean_cv_r2) |>
  dplyr::pull(scale_name)

fig_direct <- scale_prediction_tbl |>
  dplyr::mutate(
    method = factor(method, levels = observed_method_levels),
    scale_name = factor(scale_name, levels = direct_scale_order)
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = cv_r2, y = scale_name, shape = method)) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.62),
    colour = "black",
    fill = "white",
    size = 2.0,
    stroke = 0.45
  ) +
  ggplot2::scale_shape_manual(values = method_shapes[observed_method_levels], labels = method_label_lookup) +
  ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
  ggplot2::scale_y_discrete(labels = format_scale_label) +
  ggplot2::labs(x = "Direct readout R2", y = NULL) +
  theme_pub_bw(base_size = 11, y_grid = TRUE, x_grid = TRUE, legend_position = "bottom") +
  ggplot2::guides(shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

save_plot_pair(
  fig_direct,
  "SUPP_FIG_direct_scale_score_prediction_by_scale_method_bw",
  output_dir = output_dir,
  width = 7.5,
  height = 5.8,
  dpi = cfg$dpi
)
save_plot_pair(
  fig_direct,
  "SUPP_FIG_S1_direct_scale_score_readout_by_scale_bw",
  output_dir = output_dir,
  width = 7.5,
  height = 5.8,
  dpi = cfg$dpi
)

stability_metrics <- c("axis_corr_mean", "item_rmse", "scale_regression_vector_cosine")
fig_full_refit <- full_refit_stability_tbl |>
  dplyr::select(method, dplyr::all_of(stability_metrics)) |>
  tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    method = factor(method, levels = observed_method_levels),
    metric = factor(metric, levels = stability_metrics, labels = format_metric_label(stability_metrics))
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = method, y = value, fill = method)) +
  ggplot2::geom_boxplot(width = 0.62, colour = "black", outlier.shape = 21, outlier.size = 0.8, linewidth = 0.35) +
  ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
  ggplot2::scale_x_discrete(labels = method_label_lookup) +
  ggplot2::scale_fill_manual(values = method_fills[observed_method_levels], guide = "none") +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_pub_bw(base_size = 11, y_grid = TRUE, legend_position = "none") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(lineheight = 0.9))

save_plot_pair(
  fig_full_refit,
  "SUPP_FIG_split_half_full_refit_stability_bw",
  output_dir = output_dir,
  width = 9.3,
  height = 3.6,
  dpi = cfg$dpi
)

selection_metrics <- c("selection_jaccard", "aligned_weight_correlation")
fig_selection <- selection_stability_tbl |>
  dplyr::select(method, dplyr::all_of(selection_metrics)) |>
  tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    method = factor(method, levels = observed_method_levels),
    metric = factor(metric, levels = selection_metrics, labels = format_metric_label(selection_metrics))
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = method, y = value, fill = method)) +
  ggplot2::geom_boxplot(width = 0.62, colour = "black", outlier.shape = 21, outlier.size = 0.8, linewidth = 0.35) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.045, height = 0, seed = 42),
    shape = 21,
    colour = "black",
    fill = "white",
    size = 0.9,
    alpha = 0.55,
    stroke = 0.25
  ) +
  ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
  ggplot2::scale_x_discrete(labels = method_label_lookup) +
  ggplot2::scale_fill_manual(values = method_fills[observed_method_levels], guide = "none") +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_pub_bw(base_size = 11, y_grid = TRUE, legend_position = "none") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(lineheight = 0.9))

save_plot_pair(
  fig_selection,
  "SUPP_FIG_selection_jaccard_and_weight_stability_bw",
  output_dir = output_dir,
  width = 6.8,
  height = 3.6,
  dpi = cfg$dpi
)
save_plot_pair(
  fig_selection,
  "SUPP_FIG_S2_selection_stability_bw",
  output_dir = output_dir,
  width = 6.8,
  height = 3.6,
  dpi = cfg$dpi
)

id_weight_data <- item_weights_tbl |>
  dplyr::filter(selected, grepl("id_guided", method, fixed = TRUE)) |>
  dplyr::group_by(method) |>
  dplyr::arrange(dplyr::desc(weight), .by_group = TRUE) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(method = factor(method, levels = observed_method_levels))

fig_weights <- ggplot2::ggplot(id_weight_data, ggplot2::aes(x = rank, y = weight)) +
  ggplot2::geom_segment(ggplot2::aes(xend = rank, y = 0, yend = weight), colour = "black", linewidth = 0.28) +
  ggplot2::geom_point(shape = 21, fill = "white", colour = "black", size = 0.75, stroke = 0.2) +
  ggplot2::facet_wrap(~method, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(method_label_lookup)) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
  ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(x = "Selected item rank by weight", y = "Item weight") +
  theme_pub_bw(base_size = 11, y_grid = TRUE, x_grid = FALSE, legend_position = "none")

save_plot_pair(
  fig_weights,
  "SUPP_FIG_id_guided_item_weights_ranked_bw",
  output_dir = output_dir,
  width = 7.2,
  height = 5.0,
  dpi = cfg$dpi
)
save_plot_pair(
  fig_weights,
  "SUPP_FIG_S4_ranked_selected_item_weights_bw",
  output_dir = output_dir,
  width = 7.2,
  height = 5.0,
  dpi = cfg$dpi
)

corr_selected <- item_component_correlations_tbl |>
  dplyr::left_join(
    item_weights_tbl |> dplyr::select(method, item_id, selected),
    by = c("method", "item_id")
  ) |>
  dplyr::filter(isTRUE(selected) | selected) |>
  dplyr::mutate(method = factor(method, levels = observed_method_levels))

fig_corr <- ggplot2::ggplot(corr_selected, ggplot2::aes(x = u1, y = u2)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey45") +
  ggplot2::geom_point(shape = 21, fill = "grey70", colour = "black", alpha = 0.55, size = 1.05, stroke = 0.22) +
  ggplot2::facet_wrap(~method, labeller = ggplot2::as_labeller(method_label_lookup)) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = "Item-component correlation with PC1", y = "Item-component correlation with PC2") +
  theme_pub_bw(base_size = 11, y_grid = TRUE, x_grid = TRUE, legend_position = "none")

save_plot_pair(
  fig_corr,
  "SUPP_FIG_item_component_correlations_selected_items_bw",
  output_dir = output_dir,
  width = 7.4,
  height = 5.6,
  dpi = cfg$dpi
)

if (!is.null(distance_preservation_tbl) && nrow(distance_preservation_tbl)) {
  distance_metrics <- c(
    "map_distance_spearman",
    "neighbourhood_overlap_gower_to_map",
    "neighbourhood_overlap_map_to_gower",
    "pcoa2_positive_eigen_share"
  )
  distance_plot <- distance_preservation_tbl |>
    dplyr::select(method, dplyr::any_of(distance_metrics)) |>
    tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
    dplyr::mutate(
      method = factor(method, levels = observed_method_levels),
      metric = factor(
        metric,
        levels = distance_metrics,
        labels = c(
          "Gower-map distance Spearman r",
          "Gower-to-map neighbour overlap",
          "Map-to-Gower neighbour overlap",
          "PCoA two-axis positive-eigen share"
        )
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = method, y = value, fill = method)) +
    ggplot2::geom_col(width = 0.62, colour = "black", linewidth = 0.3) +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 2) +
    ggplot2::scale_x_discrete(labels = method_label_lookup) +
    ggplot2::scale_fill_manual(values = method_fills[observed_method_levels], guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_pub_bw(base_size = 11, y_grid = TRUE, legend_position = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(lineheight = 0.9))
  
  save_plot_pair(
    distance_plot,
    "SUPP_FIG_distance_geometry_checks_bw",
    output_dir = output_dir,
    width = 8.4,
    height = 5.4,
    dpi = cfg$dpi
  )
  save_plot_pair(
    distance_plot,
    "SUPP_FIG_S3_distance_geometry_checks_bw",
    output_dir = output_dir,
    width = 8.4,
    height = 5.4,
    dpi = cfg$dpi
  )
}

retained_item_supp_tbl <- build_retained_item_table(selected_item_content_tbl, item_weights_tbl)
loso_readout_supp_tbl <- build_readout_delta_table(leave_one_scale_out_tbl)
direct_readout_supp_tbl <- build_readout_delta_table(scale_prediction_tbl)
stability_summary_supp_tbl <- build_stability_summary_table(
  full_refit_stability_tbl = full_refit_stability_tbl,
  fixed_selection_stability_tbl = fixed_selection_stability_tbl,
  selection_stability_tbl = selection_stability_tbl,
  stability_metrics = stability_metrics,
  selection_metrics = selection_metrics
)
distance_geometry_supp_tbl <- if (!is.null(distance_preservation_tbl) && nrow(distance_preservation_tbl)) {
  build_distance_geometry_table(distance_preservation_tbl)
} else {
  NULL
}

write_csv2_safe(
  retained_item_supp_tbl,
  file.path(output_dir, "SUPP_TABLE_S1_retained_item_set.csv")
)
write_csv2_safe(
  loso_readout_supp_tbl,
  file.path(output_dir, "SUPP_TABLE_S2_leave_one_scale_out_readout_by_scale.csv")
)
write_csv2_safe(
  direct_readout_supp_tbl,
  file.path(output_dir, "SUPP_TABLE_S3_direct_scale_score_readout_by_scale.csv")
)
write_csv2_safe(
  stability_summary_supp_tbl,
  file.path(output_dir, "SUPP_TABLE_S4_stability_summary.csv")
)
if (!is.null(distance_geometry_supp_tbl)) {
  write_csv2_safe(
    distance_geometry_supp_tbl,
    file.path(output_dir, "SUPP_TABLE_S5_distance_geometry_checks.csv")
  )
}

write_csv2_safe(
  scale_prediction_tbl |>
    dplyr::mutate(
      method_label = format_method_label(method),
      scale_label = format_scale_label(scale_name),
      .after = method
    ),
  file.path(output_dir, "SUPP_TABLE_direct_scale_score_prediction_by_scale.csv")
)
write_csv2_safe(
  leave_one_scale_out_tbl |>
    dplyr::mutate(
      method_label = format_method_label(method),
      scale_label = format_scale_label(scale_name),
      .after = method
    ),
  file.path(output_dir, "SUPP_TABLE_leave_one_scale_out_prediction_by_scale.csv")
)
write_csv2_safe(
  summarise_metric_by_method(full_refit_stability_tbl, stability_metrics),
  file.path(output_dir, "SUPP_TABLE_full_refit_stability_summary.csv")
)
write_csv2_safe(
  full_refit_stability_tbl |>
    dplyr::mutate(method_label = format_method_label(method), .after = method),
  file.path(output_dir, "SUPP_TABLE_full_refit_stability_raw.csv")
)
write_csv2_safe(
  summarise_metric_by_method(fixed_selection_stability_tbl, stability_metrics),
  file.path(output_dir, "SUPP_TABLE_fixed_selection_stability_summary.csv")
)
write_csv2_safe(
  fixed_selection_stability_tbl |>
    dplyr::mutate(method_label = format_method_label(method), .after = method),
  file.path(output_dir, "SUPP_TABLE_fixed_selection_stability_raw.csv")
)
write_csv2_safe(
  summarise_metric_by_method(selection_stability_tbl, selection_metrics),
  file.path(output_dir, "SUPP_TABLE_selection_stability_summary.csv")
)
write_csv2_safe(
  selection_stability_tbl |>
    dplyr::mutate(method_label = format_method_label(method), .after = method),
  file.path(output_dir, "SUPP_TABLE_selection_stability_raw.csv")
)
write_csv2_safe(
  item_weights_tbl |>
    dplyr::mutate(method_label = format_method_label(method), .after = method),
  file.path(output_dir, "SUPP_TABLE_item_weights_all_methods.csv")
)
write_csv2_safe(
  item_weights_tbl |>
    dplyr::filter(selected) |>
    dplyr::mutate(method_label = format_method_label(method), .after = method),
  file.path(output_dir, "SUPP_TABLE_item_weights_selected_items.csv")
)
write_csv2_safe(
  retained_item_supp_tbl,
  file.path(output_dir, "SUPP_TABLE_selected_item_content.csv")
)
write_csv2_safe(
  item_component_correlations_tbl |>
    dplyr::mutate(method_label = format_method_label(method), .after = method),
  file.path(output_dir, "SUPP_TABLE_item_component_correlations_all_items.csv")
)
write_csv2_safe(
  corr_selected |>
    dplyr::mutate(method_label = format_method_label(as.character(method)), .after = method),
  file.path(output_dir, "SUPP_TABLE_item_component_correlations_selected_items.csv")
)

manifest_artifacts <- c(
  "TABLE_1_analytic_solution_summary.csv",
  "TABLE_1_analytic_solution_summary.md",
  "FIG_1_eigenspectrum_two_component_concentration_bw.pdf",
  "FIG_2_matched_random_subset_pc12_distribution_bw.pdf",
  if (!is.null(pc_scores_tbl) && !is.null(scale_regression_vectors_tbl)) "FIG_3_selected_item_respondent_map_bw.pdf",
  "FIG_4_leave_one_scale_out_delta_by_scale_bw.pdf",
  "SUPP_TABLE_S1_retained_item_set.csv",
  "SUPP_TABLE_S2_leave_one_scale_out_readout_by_scale.csv",
  "SUPP_TABLE_S3_direct_scale_score_readout_by_scale.csv",
  "SUPP_TABLE_S4_stability_summary.csv",
  if (!is.null(distance_geometry_supp_tbl)) "SUPP_TABLE_S5_distance_geometry_checks.csv",
  "SUPP_FIG_S1_direct_scale_score_readout_by_scale_bw.pdf",
  "SUPP_FIG_S2_selection_stability_bw.pdf",
  if (!is.null(distance_geometry_supp_tbl)) "SUPP_FIG_S3_distance_geometry_checks_bw.pdf",
  "SUPP_FIG_S4_ranked_selected_item_weights_bw.pdf",
  "SUPP_FIG_direct_scale_score_prediction_by_scale_method_bw.pdf",
  if (!is.null(distance_preservation_tbl) && nrow(distance_preservation_tbl)) "SUPP_FIG_distance_geometry_checks_bw.pdf",
  "SUPP_FIG_split_half_full_refit_stability_bw.pdf",
  "SUPP_FIG_selection_jaccard_and_weight_stability_bw.pdf",
  "SUPP_FIG_id_guided_item_weights_ranked_bw.pdf",
  "SUPP_FIG_item_component_correlations_selected_items_bw.pdf",
  "SUPP_TABLE_selected_item_content.csv"
)
manifest_descriptions <- c(
  "Formatted manuscript summary table.",
  "Markdown rendering of Table 1.",
  "Eigenspectrum and PC1+PC2 concentration across analytic solutions.",
  "Matched random-subset baseline for PC1+PC2 explained variance.",
  if (!is.null(pc_scores_tbl) && !is.null(scale_regression_vectors_tbl)) "Selected-item two-dimensional respondent map with scale-score gradients.",
  "Leave-one-scale-out readout differences by scale.",
  "Retained item set with scale membership, keying, final weights, and survivor frequencies.",
  "Leave-one-scale-out scale-score readout by scale.",
  "Direct scale-score readout by scale.",
  "Split-half fixed-selection, full-refit, and selection-stability summaries.",
  if (!is.null(distance_geometry_supp_tbl)) "Distance-geometry check values.",
  "Direct fixed-map scale-score readout by scale and method.",
  "Selection Jaccard and aligned-weight stability.",
  if (!is.null(distance_geometry_supp_tbl)) "Distance-geometry checks comparing Gower distances with two-dimensional map distances.",
  "Ranked retained-item weights.",
  "Direct fixed-map scale-score readout by scale and method.",
  if (!is.null(distance_preservation_tbl) && nrow(distance_preservation_tbl)) "Distance-geometry checks comparing Gower distances with two-dimensional map distances.",
  "Split-half full-refit stability.",
  "Selection Jaccard and aligned-weight stability.",
  "Ranked ID-guided item weights.",
  "Selected item-component correlations.",
  "Selected item content, scale membership, keying, weights, and selection frequency."
)
manifest <- tibble::tibble(
  artifact = manifest_artifacts,
  description = manifest_descriptions
)
write_csv2_safe(manifest, file.path(output_dir, "publication_bw_manifest.csv"))

message("Finished. Publication artifacts written to: ", output_dir)
