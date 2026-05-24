suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(rrcov)
  library(future)
})

default_method_sensitivity_cfg <- list(
  root_dir = NULL,
  data_path = "data/psychometric_matrix.csv",
  qbank_path = "data/QBank.csv",
  scale_spec_path = "data/16pf_scale_spec.csv",
  setup_script = "0_setup.R",
  out_dir = "out/method_sensitivity_16pf",
  bam_threads = "all_but_one",
  omp_threads = 1L,
  blas_threads = 1L,
  seed = 42L,
  min_scale_completion_prop = 0.70,
  cv_folds = 5L,
  optimisation_subsample_n = 250L,
  optimisation_w_min = 0.05,
  optimisation_multi_runs = 20L,
  selection_stability_optimisation_multi_runs = 5L,
  optimisation_multi_min_prop = 0.35,
  optimisation_step_grid = c(0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.05),
  optimisation_batch_k = 3L,
  optimisation_batch_factor = 0.75,
  optimisation_max_iter = 1000L,
  optimisation_eval_per_iter = 50L,
  stability_reps = 20L,
  stability_fixed_selection_reps = 50L,
  stability_half_n = NA_integer_,
  split_half_full_refit_parallel = TRUE,
  split_half_data_cache = TRUE,
  post_parallel_diagnostics = FALSE,
  spectrum_rank = 10L,
  baseline_reps = 500L,
  baseline_subset_modes = c("unstratified", "stratified_by_scale"),
  distance_eval_n = 1200L,
  distance_eval_k = 15L,
  ncores_par = 4L,
  future_scheduling = 1,
  save_plots = TRUE
)

`%||%` <- function(x, y) if (is.null(x)) y else x

messagef <- function(...) {
  message(sprintf(...))
}

source_gower_optimizer_core <- function(root_dir = NULL) {
  candidates <- unique(stats::na.omit(c(
    if (!is.null(root_dir) && nzchar(root_dir)) file.path(root_dir, "tools", "gower_optimizer_core.R") else NA_character_,
    file.path(getwd(), "tools", "gower_optimizer_core.R"),
    file.path(getwd(), "gower_optimizer_core.R"),
    file.path(dirname(getwd()), "tools", "gower_optimizer_core.R")
  )))
  core_path <- candidates[file.exists(candidates)][1]
  if (!length(core_path) || is.na(core_path)) {
    stop("Cannot find tools/gower_optimizer_core.R")
  }
  
  core_env <- new.env(parent = globalenv())
  source(core_path, local = core_env, chdir = FALSE)
  for (nm in ls(core_env, all.names = TRUE)) {
    assign(nm, get(nm, envir = core_env), envir = parent.frame())
  }
  invisible(core_path)
}

source_gower_optimizer_core()

limit_worker_threads <- function() {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(1)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    BLAS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  invisible(NULL)
}

future_worker_count <- function() {
  as.integer(future::nbrOfWorkers())
}

.drop_future_args <- function(dots) {
  if (is.null(names(dots))) return(dots)
  dots[!grepl("^future\\.", names(dots), perl = TRUE)]
}

load_setup_helpers <- function(cfg, root_dir) {
  setup_path <- normalise_path(cfg$setup_script, root_dir, must_work = TRUE)
  setup_env <- new.env(parent = globalenv())
  setup_env$SETUP_CFG <- list(
    BAM_THREADS = cfg$bam_threads,
    NCORES_PAR = cfg$ncores_par,
    OMP_THREADS = cfg$omp_threads,
    BLAS_THREADS = cfg$blas_threads,
    SET_ENV_THREADS = TRUE
  )
  source(setup_path, chdir = TRUE, local = setup_env)
  setup_env$OUTPUTS_DIR <- cfg$out_dir
  setup_env$SAVE_IMAGES <- isTRUE(cfg$save_plots)
  
  assign("FUTURE_LAPPLY", setup_env$FUTURE_LAPPLY, envir = environment(map_with_progress))
  invisible(list(
    nworkers = future_worker_count(),
    save_plot_gg = setup_env$save_plot_gg,
    theme_pub = setup_env$theme_pub,
    scale_prob_fill = setup_env$scale_prob_fill,
    cluster_colours = setup_env$cluster_colours
  ))
}

map_with_progress_worker <- function(x,
                                     worker_fun,
                                     dots = list(),
                                     p = function() NULL,
                                     progress_arg_name = NULL,
                                     progress_at_end = TRUE) {
  worker_dots <- dots
  if (!is.null(progress_arg_name)) {
    worker_dots[[progress_arg_name]] <- p
  }
  out <- do.call(worker_fun, c(list(x), worker_dots))
  if (isTRUE(progress_at_end)) p()
  out
}

map_with_progress <- function(X,
                              FUN,
                              cfg,
                              ...,
                              show_progress = TRUE,
                              parallel = TRUE,
                              future.seed = TRUE,
                              progress_callback = NULL,
                              progress_steps = NULL,
                              progress_arg_name = NULL,
                              progress_at_end = TRUE) {
  dots <- list(...)
  if (!is.function(progress_callback)) progress_callback <- NULL
  if (is.null(progress_steps)) progress_steps <- length(X)
  make_worker_dots <- function() {
    worker_dots <- .drop_future_args(dots)
    if (!is.null(progress_arg_name)) worker_dots[[progress_arg_name]] <- progress_callback %||% function() NULL
    worker_dots
  }
  
  run_sequential <- function() {
    if (show_progress) {
      pb <- utils::txtProgressBar(min = 0, max = length(X), style = 3)
      on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }

    out <- vector("list", length(X))
    for (i in seq_along(X)) {
      out[[i]] <- do.call(FUN, c(list(X[[i]]), make_worker_dots()))
      if (isTRUE(progress_at_end) && !is.null(progress_callback)) progress_callback()
      if (show_progress && isTRUE(progress_at_end)) utils::setTxtProgressBar(pb, i)
    }
    out
  }

  if (!parallel) {
    return(run_sequential())
  }

  nworkers <- future_worker_count()
  if (nworkers <= 1L) {
    return(run_sequential())
  }

  worker_dots <- .drop_future_args(dots)
  future_scheduling <- dots$future.scheduling %||% cfg$future_scheduling %||% 1L
  future_packages <- dots$future.packages %||% NULL
  future_globals <- dots$future.globals
  if (!is.null(progress_callback)) {
    future_packages <- unique(c(future_packages, "progressr"))
  }
  
  call_future_lapply <- function(progress_fun, packages) {
    call_args <- list(
      X = X,
      FUN = map_with_progress_worker,
      worker_fun = FUN,
      dots = worker_dots,
      p = progress_fun %||% function() NULL,
      progress_arg_name = progress_arg_name,
      progress_at_end = progress_at_end,
      future.packages = packages,
      future.seed = future.seed,
      future.scheduling = future_scheduling
    )
    if (!is.null(future_globals)) {
      call_args$future.globals <- future_globals
    }
    do.call(FUTURE_LAPPLY, call_args)
  }
  
  if (!show_progress) {
    return(call_future_lapply(progress_callback, future_packages))
  }
  
  progressr::with_progress({
    p <- progressr::progressor(steps = progress_steps)
    worker_p <- if (is.null(progress_callback)) {
      p
    } else {
      function() {
        p()
        progress_callback()
      }
    }
    call_future_lapply(worker_p, unique(c(future_packages, "progressr")))
  })
}

infer_root_dir <- function() {
  cwd <- normalizePath(".", mustWork = TRUE)
  if (file.exists(file.path(cwd, "0_setup.R")) &&
      file.exists(file.path(cwd, "1_dimension_psychometric.R"))) {
    return(cwd)
  }
  if (basename(cwd) == "tools" &&
      file.exists(file.path(cwd, "..", "0_setup.R")) &&
      file.exists(file.path(cwd, "..", "1_dimension_psychometric.R"))) {
    return(normalizePath(file.path(cwd, ".."), mustWork = TRUE))
  }
  stop("Run from the project root or pass root_dir explicitly.")
}

normalise_path <- function(path, root_dir, must_work = TRUE) {
  if (!grepl("^(/|~|[A-Za-z]:)", path)) {
    path <- file.path(root_dir, path)
  }
  normalizePath(path, mustWork = must_work)
}

write_csv <- function(x, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv2(x, file, ...)
}

parse_item_range <- function(prefix, spec) {
  spec <- trimws(spec)
  if (!nzchar(spec)) return(character(0))
  chunks <- strsplit(spec, ",", fixed = TRUE)[[1L]]
  items <- character(0)
  for (chunk in chunks) {
    chunk <- trimws(chunk)
    if (grepl(":", chunk, fixed = TRUE)) {
      bounds <- as.integer(strsplit(chunk, ":", fixed = TRUE)[[1L]])
      items <- c(items, paste0(prefix, seq.int(bounds[[1L]], bounds[[2L]])))
    } else {
      items <- c(items, paste0(prefix, as.integer(chunk)))
    }
  }
  unique(items)
}

build_scale_key <- function(scale_spec_path) {
  spec <- readr::read_csv(scale_spec_path, show_col_types = FALSE)
  rows <- vector("list", nrow(spec))
  for (i in seq_len(nrow(spec))) {
    prefix <- spec$dataset_prefix[[i]]
    pos_items <- parse_item_range(prefix, spec$positive_items[[i]])
    neg_items <- parse_item_range(prefix, spec$negative_items[[i]])
    rows[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        dataset_prefix = prefix,
        scale_name = spec$scale_name[[i]],
        original_16pf_factor = spec$original_16pf_factor[[i]],
        item_id = pos_items,
        keyed_sign = 1L
      ),
      tibble::tibble(
        dataset_prefix = prefix,
        scale_name = spec$scale_name[[i]],
        original_16pf_factor = spec$original_16pf_factor[[i]],
        item_id = neg_items,
        keyed_sign = -1L
      )
    )
  }
  dplyr::bind_rows(rows) |>
    dplyr::arrange(dataset_prefix, item_id)
}

read_qbank <- function(qbank_path) {
  qbank <- readr::read_csv2(
    qbank_path,
    col_names = c("item_id", "field_type", "item_text"),
    show_col_types = FALSE
  )
  qbank |>
    dplyr::filter(item_id != "Field", grepl("^[A-P][0-9]+$", item_id)) |>
    dplyr::mutate(
      item_text = gsub('^"+|"+$', "", item_text),
      item_text = gsub('""', '"', item_text)
    )
}

median_impute <- function(X) {
  medians <- vapply(X, function(v) stats::median(v, na.rm = TRUE), numeric(1))
  for (nm in names(X)) {
    idx <- is.na(X[[nm]])
    if (any(idx)) X[[nm]][idx] <- medians[[nm]]
  }
  list(data = X, medians = medians)
}

zscore_df <- function(X) {
  Xs <- scale(as.matrix(X), center = TRUE, scale = TRUE)
  Xs[, apply(Xs, 2, function(v) all(is.finite(v))), drop = FALSE]
}

prep_ord_gower <- function(X) {
  X1 <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(X1)) {
    x <- as.integer(round(X1[[nm]]))
    x[!is.finite(x)] <- NA_integer_
    X1[[nm]] <- factor(x, exclude = NULL)
  }
  list(X = X1, type = list())
}

gower_dist <- function(Xdf, type_list, weights) {
  cluster::daisy(Xdf, metric = "gower", type = type_list, weights = weights)
}

two_nn_from_dist <- function(D, eps = .Machine$double.eps) {
  M <- as.matrix(D)
  if (nrow(M) < 2L) return(list(d1 = numeric(0), d2 = numeric(0)))
  diag(M) <- Inf
  negM <- -M
  idx1 <- max.col(negM, ties.method = "first")
  d1 <- M[cbind(seq_len(nrow(M)), idx1)]
  negM[cbind(seq_len(nrow(M)), idx1)] <- -Inf
  idx2 <- max.col(negM, ties.method = "first")
  d2 <- M[cbind(seq_len(nrow(M)), idx2)]
  d1 <- pmax(d1, eps)
  d2 <- pmax(d2, d1 + eps)
  list(d1 = d1, d2 = d2)
}

twonn_id_from_dist <- function(D, trim = 0.02) {
  nn <- two_nn_from_dist(D)
  r <- nn$d2 / nn$d1
  r <- r[is.finite(r) & r > 1]
  if (!length(r)) return(NA_real_)
  lr <- sort(log(r))
  k <- floor(trim * length(lr))
  if (k > 0L && (2L * k) < length(lr)) {
    lr <- lr[(k + 1L):(length(lr) - k)]
  }
  1 / mean(lr)
}

make_ns_cache <- function(Xdf, type_list) {
  gower_opt_make_ns_cache(Xdf, type = type_list, gower_fun = gower_dist)
}

run_item_weight_search <- function(X,
                                   row_idx,
                                   n_rows_sub = 800L,
                                   w_min = 0.05,
                                   step_grid = c(0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.05),
                                   batch_k = 3L,
                                   batch_factor = 0.75,
                                   max_iter = 200L,
                                   eval_per_iter = 50L,
                                   seed = 42L) {
  set.seed(seed)
  row_idx <- as.integer(row_idx)
  row_idx <- row_idx[is.finite(row_idx)]

  n_sub <- suppressWarnings(as.integer(n_rows_sub))
  if (!is.finite(n_sub) || n_sub < 1L) n_sub <- length(row_idx)
  n_sub <- min(n_sub, length(row_idx))

  idx_sub <- if (length(row_idx) > n_sub) sample(row_idx, n_sub) else row_idx
  prep <- prep_ord_gower(X[idx_sub, , drop = FALSE])
  cache <- make_ns_cache(prep$X, prep$type)
  vars <- names(prep$X)
  opt <- gower_opt_stochastic_weights(
    cache = cache,
    vars = vars,
    init_weights = setNames(rep(1, length(vars)), vars),
    allow_update = rep(TRUE, length(vars)),
    id_fun = twonn_id_from_dist,
    w_min = w_min,
    step_grid = step_grid,
    batch_k = batch_k,
    batch_factor = batch_factor,
    max_iter = max_iter,
    eval_per_iter = eval_per_iter,
    verbose = FALSE,
    n_rows = cache$n
  )
  
  hist <- tibble::as_tibble(opt$history)
  hist <- dplyr::transmute(
    hist,
    iter = iter,
    twonn_id = ID,
    changed = dplyr::if_else(iter == 0L, "START", changed),
    weight = suppressWarnings(as.numeric(note))
  )
  
  list(weights = opt$weights, history = hist, twonn_id = opt$final_ID, idx_sub = idx_sub)
}

optimise_item_weights <- function(X,
                                  n_rows_sub = 800L,
                                  w_min = 0.05,
                                  step_grid = c(0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.05),
                                  batch_k = 3L,
                                  batch_factor = 0.75,
                                  max_iter = 200L,
                                  eval_per_iter = 50L,
                                  seed = 42L) {
  run_item_weight_search(
    X = X,
    row_idx = seq_len(nrow(X)),
    n_rows_sub = n_rows_sub,
    w_min = w_min,
    step_grid = step_grid,
    batch_k = batch_k,
    batch_factor = batch_factor,
    max_iter = max_iter,
    eval_per_iter = eval_per_iter,
    seed = seed
  )
}

knee_select_items <- function(weights,
                              w_min = 0.05,
                              kmin = NULL,
                              active_eps = 1e-8) {
  w_names <- names(weights)
  if (is.null(w_names)) w_names <- paste0("V", seq_along(weights))
  w <- pmax(as.numeric(weights), w_min)
  names(w) <- w_names
  w <- sort(w, decreasing = TRUE)
  n <- length(w)
  x <- seq_len(n)
  if (is.null(kmin)) kmin <- max(30L, ceiling(0.10 * n))
  active_eps <- max(active_eps, .Machine$double.eps)
  active_floor <- w_min + active_eps
  
  if (n < 3L) {
    survivors <- names(w)
    active_survivors <- survivors[w > active_floor]
    return(list(
      survivors = survivors,
      active_survivors = active_survivors,
      threshold = if (length(w)) min(w) else NA_real_,
      weights_sorted = w
    ))
  }
  
  w_range <- max(w, na.rm = TRUE) - min(w, na.rm = TRUE)
  if (!is.finite(w_range) || w_range <= .Machine$double.eps) {
    cut_idx <- if (all(w > active_floor)) n else min(n, kmin)
    survivors <- names(w)[seq_len(cut_idx)]
    active_survivors <- survivors[w[seq_len(cut_idx)] > active_floor]
    return(list(
      survivors = survivors,
      active_survivors = active_survivors,
      threshold = w[[cut_idx]],
      weights_sorted = w
    ))
  }
  
  y_norm <- (w - min(w)) / w_range
  x_norm <- (x - min(x)) / (max(x) - min(x))
  d <- abs(y_norm - (1 - x_norm))
  plateau_idx <- which(w >= 0.95)
  plateau_end <- if (length(plateau_idx)) max(plateau_idx) else NA_integer_
  if (!is.finite(plateau_end) || plateau_end == n) plateau_end <- 1L
  knee_local <- which.max(d[plateau_end:n]) + plateau_end - 1L
  keep_n <- min(n, max(knee_local, kmin))
  survivors <- names(w)[seq_len(keep_n)]
  active_survivors <- survivors[w[seq_len(keep_n)] > active_floor]
  list(
    survivors = survivors,
    active_survivors = active_survivors,
    threshold = w[[keep_n]],
    weights_sorted = w
  )
}

combine_multi_run_weights <- function(opt_list, cfg) {
  ref <- opt_list[[1L]]
  var_names <- names(ref$weights)
  n_runs <- length(opt_list)
  
  sel_info <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    sel_info[[r]] <- knee_select_items(opt_list[[r]]$weights, w_min = cfg$optimisation_w_min)
  }
  
  Wmat <- vapply(
    seq_len(n_runs),
    function(r) opt_list[[r]]$weights[var_names],
    numeric(length(var_names))
  )
  dimnames(Wmat) <- list(var_names, paste0("run", seq_len(n_runs)))
  
  sel_mat <- matrix(
    0L,
    nrow = length(var_names),
    ncol = n_runs,
    dimnames = list(var_names, paste0("run", seq_len(n_runs)))
  )
  for (r in seq_len(n_runs)) {
    active_r <- sel_info[[r]]$active_survivors %||% sel_info[[r]]$survivors
    sel_mat[active_r, r] <- 1L
  }
  
  stability_tbl <- tibble::tibble(
    item_id = var_names,
    selected_count = rowSums(sel_mat),
    selected_prop = rowSums(sel_mat) / n_runs,
    weight_ref = ref$weights[var_names],
    weight_mean = rowMeans(Wmat, na.rm = TRUE),
    weight_sd = apply(Wmat, 1L, stats::sd, na.rm = TRUE),
    weight_max = apply(Wmat, 1L, max, na.rm = TRUE)
  ) |>
    dplyr::arrange(dplyr::desc(selected_prop), dplyr::desc(weight_mean))
  
  active_floor <- cfg$optimisation_w_min + 1e-8
  robust_mask <- stability_tbl$selected_prop >= cfg$optimisation_multi_min_prop &
    stability_tbl$weight_max > active_floor
  
  if (!any(robust_mask)) {
    robust_set <- sel_info[[1L]]$survivors
    robust_mask <- stability_tbl$item_id %in% robust_set
  } else {
    robust_set <- stability_tbl$item_id[robust_mask]
  }
  
  weight_raw <- with(stability_tbl, weight_mean - cfg$optimisation_w_min * (1 - selected_prop))
  weight_raw[!is.finite(weight_raw)] <- 0
  weight_raw <- pmax(0, weight_raw)
  scale_fac <- max(weight_raw[robust_mask], na.rm = TRUE)
  if (!is.finite(scale_fac) || scale_fac <= 0) scale_fac <- 1
  
  weights <- stats::setNames(as.numeric(weight_raw / scale_fac), stability_tbl$item_id)
  
  list(
    weights = weights,
    selected_items = robust_set,
    stability = stability_tbl,
    per_run_weights = Wmat
  )
}

spectral_effective_rank <- function(values) {
  values <- values[is.finite(values) & values > 0]
  if (!length(values)) return(NA_real_)
  p <- values / sum(values)
  exp(-sum(p * log(p)))
}

fit_component_solution <- function(X_std_all,
                                   selected_items,
                                   item_weights,
                                   decomp_method = c("pca", "robpca"),
                                   spectrum_rank = 10L) {
  decomp_method <- match.arg(decomp_method)
  X_sel <- X_std_all[, selected_items, drop = FALSE]
  w_sel <- item_weights[selected_items]
  X_weighted <- sweep(X_sel, 2, sqrt(w_sel), "*")
  rank_eff <- min(as.integer(spectrum_rank), ncol(X_weighted), nrow(X_weighted) - 1L)
  if (rank_eff < 2L) stop("Need at least 2 components to fit the score map.")
  if (identical(decomp_method, "pca")) {
    fit <- prcomp(X_weighted, center = FALSE, scale. = FALSE, rank. = rank_eff)
    pc_scores_2d <- fit$x[, 1:2, drop = FALSE]
    spectrum <- fit$sdev[seq_len(rank_eff)]^2
  } else {
    fit <- rrcov::PcaHubert(
      x = X_weighted,
      k = rank_eff,
      kmax = rank_eff,
      scale = FALSE,
      signflip = TRUE
    )
    pc_scores_2d <- fit@scores[, 1:2, drop = FALSE]
    spectrum <- fit@eigenvalues[seq_len(rank_eff)]
  }
  colnames(pc_scores_2d) <- c("u1", "u2")
  total_var <- sum(apply(X_weighted, 2, stats::var))
  explained_variance_ratio <- spectrum / pmax(total_var, 1e-12)
  item_component_correlations <- cor(X_std_all, pc_scores_2d, use = "pairwise.complete.obs")
  item_component_correlations <- item_component_correlations[, 1:2, drop = FALSE]
  rownames(item_component_correlations) <- colnames(X_std_all)
  list(
    selected_items = selected_items,
    weights = item_weights,
    pc_scores_2d = pc_scores_2d,
    spectrum = spectrum,
    explained_variance_ratio = explained_variance_ratio,
    spectral_effective_rank = spectral_effective_rank(spectrum),
    item_component_correlations = item_component_correlations
  )
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  if (stats::sd(x) <= 1e-12 || stats::sd(y) <= 1e-12) return(NA_real_)
  suppressWarnings(stats::cor(x, y))
}

cumulative_explained_variance_ratio <- function(explained_variance_ratio, k) {
  if (!length(explained_variance_ratio)) return(NA_real_)
  sum(explained_variance_ratio[seq_len(min(k, length(explained_variance_ratio)))], na.rm = TRUE)
}

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

quantile_or_na <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE))
}

prepare_standardised_matrix <- function(X_raw) {
  X_imp <- median_impute(X_raw)$data
  X_std <- zscore_df(X_imp)
  list(std = X_std)
}

fit_fixed_selection_from_prepared <- function(X_std_all,
                                              selected_items,
                                              item_weights,
                                              decomp_method,
                                              spectrum_rank) {
  selected_available <- intersect(selected_items, colnames(X_std_all))
  if (length(selected_available) < 2L) {
    stop("Need at least two available selected items to fit a fixed-selection score map.")
  }
  weight_vec <- item_weights
  if (is.null(names(weight_vec))) {
    names(weight_vec) <- colnames(X_std_all)[seq_along(weight_vec)]
  }
  missing_weights <- setdiff(selected_available, names(weight_vec))
  if (length(missing_weights)) {
    extra <- setNames(rep(1, length(missing_weights)), missing_weights)
    weight_vec <- c(weight_vec, extra)
  }
  weight_vec <- weight_vec[selected_available]
  weight_vec[!is.finite(weight_vec) | weight_vec <= 0] <- 1
  fit_component_solution(
    X_std_all = X_std_all,
    selected_items = selected_available,
    item_weights = weight_vec,
    decomp_method = decomp_method,
    spectrum_rank = spectrum_rank
  )
}

selection_stability_metrics <- function(fit_a, fit_b) {
  sel_a <- unique(as.character(fit_a$selected_items))
  sel_b <- unique(as.character(fit_b$selected_items))
  inter_items <- intersect(sel_a, sel_b)
  union_items <- union(sel_a, sel_b)
  
  all_weight_items <- union(names(fit_a$weights), names(fit_b$weights))
  wa <- setNames(rep(0, length(all_weight_items)), all_weight_items)
  wb <- wa
  wa[names(fit_a$weights)] <- as.numeric(fit_a$weights)
  wb[names(fit_b$weights)] <- as.numeric(fit_b$weights)
  
  aligned_weight_correlation <- safe_cor(wa, wb)
  if (!is.finite(aligned_weight_correlation) &&
      identical(fit_a$weighting_mode, "uniform") &&
      identical(fit_b$weighting_mode, "uniform")) {
    aligned_weight_correlation <- 1
  }
  
  tibble::tibble(
    selection_jaccard = if (length(union_items)) length(inter_items) / length(union_items) else NA_real_,
    aligned_weight_correlation = aligned_weight_correlation,
    n_selected_a = length(sel_a),
    n_selected_b = length(sel_b),
    n_selected_intersection = length(inter_items)
  )
}

geometry_stability_metrics <- function(fit_a,
                                       fit_b,
                                       scale_scores_a,
                                       scale_scores_b,
                                       anchor_items) {
  anchor_items <- intersect(
    anchor_items,
    intersect(
      rownames(fit_a$item_component_correlations),
      rownames(fit_b$item_component_correlations)
    )
  )
  if (length(anchor_items) < 2L) {
    return(tibble::tibble(
      axis_corr_mean = NA_real_,
      item_rmse = NA_real_,
      scale_regression_vector_cosine = NA_real_,
      n_anchor_items = length(anchor_items)
    ))
  }
  
  item_a <- fit_a$item_component_correlations[anchor_items, , drop = FALSE]
  item_b <- fit_b$item_component_correlations[anchor_items, , drop = FALSE]
  Q <- orthogonal_procrustes(item_a, item_b)
  item_b_aligned <- item_b %*% Q
  
  axis_corrs <- c(
    safe_cor(item_a[, 1], item_b_aligned[, 1]),
    safe_cor(item_a[, 2], item_b_aligned[, 2])
  )
  item_rmse <- sqrt(mean((item_a - item_b_aligned)^2))
  
  vec_a <- scale_regression_vectors(scale_scores_a, fit_a$pc_scores_2d)
  vec_b <- scale_regression_vectors(scale_scores_b, fit_b$pc_scores_2d %*% Q)
  vec_join <- dplyr::inner_join(vec_a, vec_b, by = "scale_name", suffix = c("_a", "_b"))
  cosines <- apply(vec_join, 1, function(r) {
    a <- c(as.numeric(r[["beta1_a"]]), as.numeric(r[["beta2_a"]]))
    b <- c(as.numeric(r[["beta1_b"]]), as.numeric(r[["beta2_b"]]))
    if (any(!is.finite(c(a, b)))) return(NA_real_)
    den <- sqrt(sum(a^2) * sum(b^2))
    if (!is.finite(den) || den <= 1e-12) return(NA_real_)
    sum(a * b) / den
  })
  
  tibble::tibble(
    axis_corr_mean = mean_or_na(abs(axis_corrs)),
    item_rmse = item_rmse,
    scale_regression_vector_cosine = mean_or_na(cosines),
    n_anchor_items = length(anchor_items)
  )
}

make_cv_folds <- function(n, k, seed) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

cv_r2_for_scale <- function(y, pc_scores_2d, cv_folds, seed) {
  ok <- complete.cases(y, pc_scores_2d)
  if (sum(ok) < max(30L, cv_folds * 5L)) return(NA_real_)
  yy <- as.numeric(y[ok])
  UU <- as.data.frame(pc_scores_2d[ok, , drop = FALSE])
  dat <- data.frame(y = yy, UU)
  fold_id <- make_cv_folds(length(yy), cv_folds, seed)
  pred <- rep(NA_real_, length(yy))
  for (f in seq_len(cv_folds)) {
    tr <- fold_id != f
    te <- fold_id == f
    fit <- lm(y ~ u1 + u2, data = dat[tr, , drop = FALSE])
    pred[te] <- predict(fit, newdata = dat[te, c("u1", "u2"), drop = FALSE])
  }
  mse <- mean((yy - pred)^2)
  var_y <- stats::var(yy)
  if (!is.finite(var_y) || var_y <= 1e-12) return(NA_real_)
  max(0, 1 - mse / var_y)
}

scale_regression_vectors <- function(scale_scores, pc_scores_2d) {
  rows <- vector("list", ncol(scale_scores))
  for (i in seq_along(scale_scores)) {
    y <- scale_scores[[i]]
    ok <- complete.cases(y, pc_scores_2d)
    if (sum(ok) < 10L) {
      rows[[i]] <- tibble::tibble(
        scale_name = names(scale_scores)[[i]],
        beta1 = NA_real_,
        beta2 = NA_real_,
        r2 = NA_real_
      )
      next
    }
    yy_raw <- as.numeric(y[ok])
    if (!is.finite(stats::sd(yy_raw)) || stats::sd(yy_raw) <= 1e-12) {
      rows[[i]] <- tibble::tibble(
        scale_name = names(scale_scores)[[i]],
        beta1 = NA_real_,
        beta2 = NA_real_,
        r2 = NA_real_
      )
      next
    }
    yy <- as.numeric(scale(yy_raw))
    dat <- data.frame(y = yy, u1 = pc_scores_2d[ok, 1], u2 = pc_scores_2d[ok, 2])
    fit <- lm(y ~ u1 + u2, data = dat)
    coefs <- coef(fit)
    rows[[i]] <- tibble::tibble(
      scale_name = names(scale_scores)[[i]],
      beta1 = unname(coefs[["u1"]]),
      beta2 = unname(coefs[["u2"]]),
      r2 = summary(fit)$r.squared
    )
  }
  dplyr::bind_rows(rows)
}

orthogonal_procrustes <- function(target, source) {
  target <- as.matrix(target)
  source <- as.matrix(source)
  sv <- svd(t(source) %*% target)
  sv$u %*% t(sv$v)
}

resolve_split_half_n <- function(n, half_n) {
  n_max <- floor(as.integer(n) / 2L)
  if (!is.finite(n_max) || n_max < 1L) {
    stop("Need at least two rows for split-half stability.")
  }
  half_n <- suppressWarnings(as.integer(half_n))
  if (!length(half_n) || is.na(half_n) || !is.finite(half_n) || half_n <= 0L) {
    return(n_max)
  }
  min(half_n, n_max)
}

resolve_stability_reps <- function(cfg, kind = c("full_refit", "fixed_selection")) {
  kind <- match.arg(kind)
  if (identical(kind, "fixed_selection") && !is.null(cfg$stability_fixed_selection_reps)) {
    out <- suppressWarnings(as.integer(cfg$stability_fixed_selection_reps))
    if (length(out) && is.finite(out) && out > 0L) return(out)
  }
  out <- suppressWarnings(as.integer(cfg$stability_reps))
  if (!length(out) || is.na(out) || !is.finite(out) || out <= 0L) {
    stop("Stability repetitions must be a positive integer.")
  }
  out
}

resolve_selection_stability_optimisation_runs <- function(cfg) {
  val <- cfg$selection_stability_optimisation_multi_runs
  if (is.null(val) || !length(val)) {
    val <- cfg$optimisation_multi_runs
  }
  out <- suppressWarnings(as.integer(val[[1L]]))
  if (!length(out) || is.na(out) || !is.finite(out) || out <= 0L) {
    out <- suppressWarnings(as.integer(cfg$optimisation_multi_runs[[1L]]))
  }
  if (!length(out) || is.na(out) || !is.finite(out) || out <= 0L) out <- 1L
  out
}

align_fit_to_reference <- function(fit, reference_item_component_correlations) {
  common_items <- intersect(
    rownames(reference_item_component_correlations),
    rownames(fit$item_component_correlations)
  )
  if (length(common_items) < 2L) {
    stop("Need at least two common items to align fit ", fit$method)
  }
  Q <- orthogonal_procrustes(
    reference_item_component_correlations[common_items, , drop = FALSE],
    fit$item_component_correlations[common_items, , drop = FALSE]
  )
  fit$pc_scores_2d <- fit$pc_scores_2d %*% Q
  fit$item_component_correlations <- fit$item_component_correlations %*% Q
  colnames(fit$pc_scores_2d) <- c("u1", "u2")
  colnames(fit$item_component_correlations) <- c("u1", "u2")
  fit$rotation <- Q
  fit
}

score_scales <- function(X_raw, scale_key, min_completion_prop = 0.70) {
  split_key <- split(scale_key, scale_key$scale_name)
  out <- vector("list", length(split_key))
  names(out) <- names(split_key)
  for (nm in names(split_key)) {
    key_df <- split_key[[nm]]
    item_mat <- as.matrix(X_raw[, key_df$item_id, drop = FALSE])
    keyed_mat <- item_mat
    neg_idx <- key_df$keyed_sign < 0
    if (any(neg_idx)) keyed_mat[, neg_idx] <- 6 - keyed_mat[, neg_idx]
    n_obs <- rowSums(!is.na(keyed_mat))
    min_obs <- ceiling(min_completion_prop * ncol(keyed_mat))
    sc <- rowMeans(keyed_mat, na.rm = TRUE)
    sc[n_obs < min_obs] <- NA_real_
    out[[nm]] <- sc
  }
  as.data.frame(out, check.names = FALSE)
}

method_label <- function(weighting_mode, decomp_method) {
  paste(decomp_method, weighting_mode, sep = "_")
}

optimise_item_weights_run_worker <- function(r, X_raw, worker_cfg) {
  seed_r <- if (r == 1L) worker_cfg$seed else worker_cfg$seed + r
  optimise_item_weights(
    X = X_raw,
    n_rows_sub = worker_cfg$optimisation_subsample_n,
    w_min = worker_cfg$optimisation_w_min,
    step_grid = worker_cfg$optimisation_step_grid,
    batch_k = worker_cfg$optimisation_batch_k,
    batch_factor = worker_cfg$optimisation_batch_factor,
    max_iter = worker_cfg$optimisation_max_iter,
    eval_per_iter = worker_cfg$optimisation_eval_per_iter,
    seed = seed_r
  )
}

optimise_item_weights_indexed <- function(X_raw,
                                          row_idx,
                                          n_rows_sub = 800L,
                                          w_min = 0.05,
                                          step_grid = c(0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.05),
                                          batch_k = 3L,
                                          batch_factor = 0.75,
                                          max_iter = 200L,
                                          eval_per_iter = 50L,
                                          seed = 42L) {
  run_item_weight_search(
    X = X_raw,
    row_idx = row_idx,
    n_rows_sub = n_rows_sub,
    w_min = w_min,
    step_grid = step_grid,
    batch_k = batch_k,
    batch_factor = batch_factor,
    max_iter = max_iter,
    eval_per_iter = eval_per_iter,
    seed = seed
  )
}

optimise_split_half_side_opt_list <- function(X_raw,
                                              row_idx,
                                              split_rep,
                                              side,
                                              worker_cfg,
                                              progress_fun = NULL) {
  optimisation_runs <- as.integer(worker_cfg$optimisation_multi_runs)
  if (!is.finite(optimisation_runs) || optimisation_runs < 1L) optimisation_runs <- 1L
  side_offset <- if (identical(side, "b")) 2000L else 1000L
  out <- vector("list", optimisation_runs)
  
  for (r in seq_len(optimisation_runs)) {
    seed_r <- as.integer(worker_cfg$seed) +
      100000L * as.integer(split_rep) +
      side_offset +
      as.integer(r)
    opt <- optimise_item_weights_indexed(
      X_raw = X_raw,
      row_idx = row_idx,
      n_rows_sub = worker_cfg$optimisation_subsample_n,
      w_min = worker_cfg$optimisation_w_min,
      step_grid = worker_cfg$optimisation_step_grid,
      batch_k = worker_cfg$optimisation_batch_k,
      batch_factor = worker_cfg$optimisation_batch_factor,
      max_iter = worker_cfg$optimisation_max_iter,
      eval_per_iter = worker_cfg$optimisation_eval_per_iter,
      seed = seed_r
    )
    # Workers only need weights here; histories make parallel returns much larger.
    out[[r]] <- list(weights = opt$weights)
    if (is.function(progress_fun)) progress_fun()
  }
  
  out
}

fit_id_guided_from_compact_opt_list <- function(X_std_all,
                                                decomp_method,
                                                cfg,
                                                opt_list) {
  if (!length(opt_list)) stop("Empty id-guided optimisation list.")
  
  if (length(opt_list) > 1L) {
    combined <- combine_multi_run_weights(opt_list, cfg)
    selected_items <- intersect(combined$selected_items, colnames(X_std_all))
    weights <- setNames(rep(0, ncol(X_std_all)), colnames(X_std_all))
    nm <- intersect(names(combined$weights), colnames(X_std_all))
    weights[nm] <- combined$weights[nm]
    multi_run_summary <- NULL
  } else {
    opt <- opt_list[[1L]]
    selected <- knee_select_items(opt$weights, w_min = cfg$optimisation_w_min)
    selected_items <- intersect(selected$survivors, colnames(X_std_all))
    weights <- setNames(rep(cfg$optimisation_w_min, ncol(X_std_all)), colnames(X_std_all))
    nm <- intersect(names(opt$weights), colnames(X_std_all))
    weights[nm] <- opt$weights[nm]
    multi_run_summary <- NULL
  }
  
  fit <- fit_component_solution(
    X_std_all = X_std_all,
    selected_items = selected_items,
    item_weights = weights,
    decomp_method = decomp_method,
    spectrum_rank = cfg$spectrum_rank
  )
  fit$history <- tibble::tibble(
    iter = integer(),
    twonn_id = numeric(),
    changed = character(),
    weight = numeric()
  )
  fit$weighting_mode <- "id_guided"
  fit$decomp_method <- decomp_method
  fit$method <- method_label("id_guided", decomp_method)
  fit$multi_run_summary <- multi_run_summary
  fit
}

split_half_id_guided_full_refit_worker <- function(task,
                                                   decomp_method,
                                                   worker_cfg,
                                                   X_raw = NULL,
                                                   scale_scores = NULL,
                                                   data_cache_path = NULL,
                                                   progress_fun = NULL) {
  limit_worker_threads()
  
  if (is.null(X_raw) || is.null(scale_scores)) {
    dat <- read_split_half_data_cache(data_cache_path)
    X_raw <- dat$X_raw
    scale_scores <- dat$scale_scores
  }
  
  prep_a <- prepare_standardised_matrix(X_raw[task$idx_a, , drop = FALSE])
  prep_b <- prepare_standardised_matrix(X_raw[task$idx_b, , drop = FALSE])
  common_items <- intersect(colnames(prep_a$std), colnames(prep_b$std))
  
  Xa_std <- prep_a$std[, common_items, drop = FALSE]
  Xb_std <- prep_b$std[, common_items, drop = FALSE]
  
  opt_a <- tryCatch(
    optimise_split_half_side_opt_list(
      X_raw = X_raw,
      row_idx = task$idx_a,
      split_rep = task$rep,
      side = "a",
      worker_cfg = worker_cfg,
      progress_fun = progress_fun
    ),
    error = function(e) NULL
  )
  opt_b <- tryCatch(
    optimise_split_half_side_opt_list(
      X_raw = X_raw,
      row_idx = task$idx_b,
      split_rep = task$rep,
      side = "b",
      worker_cfg = worker_cfg,
      progress_fun = progress_fun
    ),
    error = function(e) NULL
  )
  
  fit_a <- if (is.null(opt_a)) NULL else tryCatch(
    fit_id_guided_from_compact_opt_list(
      X_std_all = Xa_std,
      decomp_method = decomp_method,
      cfg = worker_cfg,
      opt_list = opt_a
    ),
    error = function(e) NULL
  )
  fit_b <- if (is.null(opt_b)) NULL else tryCatch(
    fit_id_guided_from_compact_opt_list(
      X_std_all = Xb_std,
      decomp_method = decomp_method,
      cfg = worker_cfg,
      opt_list = opt_b
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_a) || is.null(fit_b)) {
    return(empty_stability_row(task$rep))
  }
  
  geom_metrics <- geometry_stability_metrics(
    fit_a = fit_a,
    fit_b = fit_b,
    scale_scores_a = scale_scores[task$idx_a, , drop = FALSE],
    scale_scores_b = scale_scores[task$idx_b, , drop = FALSE],
    anchor_items = common_items
  )
  sel_metrics <- selection_stability_metrics(fit_a, fit_b)
  
  dplyr::bind_cols(tibble::tibble(rep = task$rep), geom_metrics, sel_metrics)
}

fit_one_method <- function(X_raw,
                           X_std_all,
                           weighting_mode,
                           decomp_method,
                           cfg) {
  limit_worker_threads()
  if (identical(weighting_mode, "id_guided")) {
    run_ids <- seq_len(as.integer(cfg$optimisation_multi_runs))
    progress_callback <- cfg$progress_callback
    if (!is.function(progress_callback)) progress_callback <- NULL
    worker_cfg <- cfg
    worker_cfg$progress_callback <- NULL
    inner_parallel <- length(run_ids) > 1L &&
      !isTRUE(cfg$disable_inner_parallel) &&
      future_worker_count() > 1L
    opt_list <- if (inner_parallel) {
      map_with_progress(
        X = as.list(run_ids),
        FUN = optimise_item_weights_run_worker,
        cfg = cfg,
        X_raw = X_raw,
        worker_cfg = worker_cfg,
        show_progress = FALSE,
        parallel = TRUE,
        future.seed = TRUE,
        future.scheduling = cfg$future_scheduling %||% 1L,
        progress_callback = progress_callback
      )
    } else {
      out <- vector("list", length(run_ids))
      for (ii in seq_along(run_ids)) {
        out[[ii]] <- optimise_item_weights_run_worker(
          run_ids[[ii]],
          X_raw = X_raw,
          worker_cfg = worker_cfg
        )
        if (!is.null(progress_callback)) progress_callback()
      }
      out
    }
    if (length(opt_list) > 1L) {
      opt <- opt_list[[1L]]
      combined <- combine_multi_run_weights(opt_list, cfg)
      selected_items <- combined$selected_items
      weights <- setNames(rep(0, ncol(X_std_all)), colnames(X_std_all))
      weights[names(combined$weights)] <- combined$weights
      history <- opt$history
      multi_run_summary <- combined
    } else {
      opt <- opt_list[[1L]]
      selected <- knee_select_items(opt$weights, w_min = cfg$optimisation_w_min)
      selected_items <- selected$survivors
      weights <- setNames(rep(cfg$optimisation_w_min, ncol(X_std_all)), colnames(X_std_all))
      weights[names(opt$weights)] <- opt$weights
      history <- opt$history
      multi_run_summary <- NULL
    }
  } else {
    selected_items <- colnames(X_std_all)
    weights <- setNames(rep(1, ncol(X_std_all)), colnames(X_std_all))
    history <- tibble::tibble(iter = 0L, twonn_id = NA_real_, changed = "UNIFORM", weight = NA_real_)
    multi_run_summary <- NULL
  }
  fit <- fit_component_solution(
    X_std_all = X_std_all,
    selected_items = selected_items,
    item_weights = weights,
    decomp_method = decomp_method,
    spectrum_rank = cfg$spectrum_rank
  )
  fit$history <- history
  fit$weighting_mode <- weighting_mode
  fit$decomp_method <- decomp_method
  fit$method <- method_label(weighting_mode, decomp_method)
  fit$multi_run_summary <- multi_run_summary
  fit
}

scale_prediction_worker <- function(j, method, pc_scores_2d, scale_scores, cv_folds, seed) {
  tibble::tibble(
    method = method,
    scale_name = names(scale_scores)[[j]],
    cv_r2 = cv_r2_for_scale(scale_scores[[j]], pc_scores_2d, cv_folds, seed + j)
  )
}

leave_one_scale_out_prediction_worker <- function(j,
                                                  method,
                                                  weighting_mode,
                                                  decomp_method,
                                                  reference_selected_items,
                                                  reference_weights,
                                                  scale_scores,
                                                  scale_items_by_name,
                                                  X_std_all,
                                                  available_items,
                                                  uniform_weights,
                                                  spectrum_rank,
                                                  cv_folds,
                                                  seed) {
  scale_name <- names(scale_scores)[[j]]
  holdout_items <- unique(scale_items_by_name[[scale_name]])
  if (identical(weighting_mode, "uniform")) {
    selected_items <- setdiff(available_items, holdout_items)
    item_weights <- uniform_weights
  } else {
    selected_items <- setdiff(reference_selected_items, holdout_items)
    item_weights <- reference_weights
  }
  
  fit_loso <- tryCatch(
    fit_fixed_selection_from_prepared(
      X_std_all = X_std_all,
      selected_items = selected_items,
      item_weights = item_weights,
      decomp_method = decomp_method,
      spectrum_rank = spectrum_rank
    ),
    error = function(e) NULL
  )
  
  tibble::tibble(
    method = method,
    scale_name = scale_name,
    n_items_selected = if (is.null(fit_loso)) NA_integer_ else length(fit_loso$selected_items),
    cv_r2 = if (is.null(fit_loso)) NA_real_ else cv_r2_for_scale(
      scale_scores[[j]],
      fit_loso$pc_scores_2d,
      cv_folds,
      seed + j
    )
  )
}

empty_stability_row <- function(rep,
                                selection = TRUE,
                                n_fixed_items_used = NULL) {
  out <- tibble::tibble(
    rep = rep,
    axis_corr_mean = NA_real_,
    item_rmse = NA_real_,
    scale_regression_vector_cosine = NA_real_,
    n_anchor_items = NA_integer_
  )
  if (!is.null(n_fixed_items_used)) {
    out$n_fixed_items_used <- n_fixed_items_used
  }
  if (isTRUE(selection)) {
    out <- dplyr::bind_cols(
      out,
      tibble::tibble(
        selection_jaccard = NA_real_,
        aligned_weight_correlation = NA_real_,
        n_selected_a = NA_integer_,
        n_selected_b = NA_integer_,
        n_selected_intersection = NA_integer_
      )
    )
  }
  out
}

split_half_full_pipeline_refit_worker <- function(b,
                                                  weighting_mode,
                                                  decomp_method,
                                                  worker_cfg,
                                                  X_raw = NULL,
                                                  scale_scores = NULL,
                                                  data_cache_path = NULL) {
  limit_worker_threads()
  
  b <- materialise_split_half_task(
    b,
    X_raw = X_raw,
    scale_scores = scale_scores,
    data_cache_path = data_cache_path
  )
  
  Xa_raw <- b$Xa_raw
  Xb_raw <- b$Xb_raw
  
  prep_a <- prepare_standardised_matrix(Xa_raw)
  prep_b <- prepare_standardised_matrix(Xb_raw)
  
  common_items <- intersect(colnames(prep_a$std), colnames(prep_b$std))
  
  Xa_raw <- Xa_raw[, common_items, drop = FALSE]
  Xb_raw <- Xb_raw[, common_items, drop = FALSE]
  Xa_std <- prep_a$std[, common_items, drop = FALSE]
  Xb_std <- prep_b$std[, common_items, drop = FALSE]
  
  fit_a <- tryCatch(
    fit_one_method(Xa_raw, Xa_std, weighting_mode, decomp_method, worker_cfg),
    error = function(e) NULL
  )
  
  fit_b <- tryCatch(
    fit_one_method(Xb_raw, Xb_std, weighting_mode, decomp_method, worker_cfg),
    error = function(e) NULL
  )
  
  if (is.null(fit_a) || is.null(fit_b)) {
    return(empty_stability_row(b$rep))
  }
  
  geom_metrics <- geometry_stability_metrics(
    fit_a = fit_a,
    fit_b = fit_b,
    scale_scores_a = b$scale_scores_a,
    scale_scores_b = b$scale_scores_b,
    anchor_items = common_items
  )
  
  sel_metrics <- selection_stability_metrics(fit_a, fit_b)
  
  dplyr::bind_cols(tibble::tibble(rep = b$rep), geom_metrics, sel_metrics)
}

split_half_fixed_selection_worker <- function(b,
                                              fixed_items,
                                              fixed_weights,
                                              decomp_method,
                                              worker_cfg,
                                              X_raw = NULL,
                                              scale_scores = NULL,
                                              data_cache_path = NULL) {
  limit_worker_threads()
  
  b <- materialise_split_half_task(
    b,
    X_raw = X_raw,
    scale_scores = scale_scores,
    data_cache_path = data_cache_path
  )
  
  prep_a <- prepare_standardised_matrix(b$Xa_raw)
  prep_b <- prepare_standardised_matrix(b$Xb_raw)
  
  common_selected <- intersect(
    fixed_items,
    intersect(colnames(prep_a$std), colnames(prep_b$std))
  )
  
  fit_a <- tryCatch(
    fit_fixed_selection_from_prepared(
      X_std_all = prep_a$std,
      selected_items = common_selected,
      item_weights = fixed_weights,
      decomp_method = decomp_method,
      spectrum_rank = worker_cfg$spectrum_rank
    ),
    error = function(e) NULL
  )
  
  fit_b <- tryCatch(
    fit_fixed_selection_from_prepared(
      X_std_all = prep_b$std,
      selected_items = common_selected,
      item_weights = fixed_weights,
      decomp_method = decomp_method,
      spectrum_rank = worker_cfg$spectrum_rank
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_a) || is.null(fit_b)) {
    return(empty_stability_row(
      b$rep,
      selection = FALSE,
      n_fixed_items_used = length(common_selected)
    ))
  }
  
  geom_metrics <- geometry_stability_metrics(
    fit_a = fit_a,
    fit_b = fit_b,
    scale_scores_a = b$scale_scores_a,
    scale_scores_b = b$scale_scores_b,
    anchor_items = common_selected
  )
  
  dplyr::bind_cols(
    tibble::tibble(rep = b$rep, n_fixed_items_used = length(common_selected)),
    geom_metrics
  )
}

build_split_half_tasks <- function(X_raw, scale_scores, reps, half_n, seed) {
  n <- nrow(X_raw)
  
  if (!is.null(scale_scores) && nrow(scale_scores) != n) {
    stop("X_raw and scale_scores must have the same number of rows.")
  }
  
  half_n <- resolve_split_half_n(n, half_n)
  
  lapply(seq_len(as.integer(reps)), function(b) {
    set.seed(seed + b)
    idx <- sample(seq_len(n), size = 2L * half_n)
    
    list(
      rep = b,
      idx_a = idx[seq_len(half_n)],
      idx_b = idx[half_n + seq_len(half_n)]
    )
  })
}

.split_half_data_cache <- new.env(parent = emptyenv())

read_split_half_data_cache <- function(data_cache_path) {
  if (is.null(data_cache_path) || !nzchar(data_cache_path)) {
    stop("data_cache_path is required for indexed split-half tasks without in-memory data.")
  }
  key <- normalizePath(data_cache_path, mustWork = TRUE)
  if (!exists(key, envir = .split_half_data_cache, inherits = FALSE)) {
    assign(key, readRDS(key), envir = .split_half_data_cache)
  }
  get(key, envir = .split_half_data_cache, inherits = FALSE)
}

materialise_split_half_task <- function(task,
                                        X_raw = NULL,
                                        scale_scores = NULL,
                                        data_cache_path = NULL) {
  if (!is.null(task$Xa_raw)) return(task)
  if (is.null(task$idx_a) || is.null(task$idx_b)) {
    stop("Split-half task must contain either raw half matrices or row indices.")
  }
  if (is.null(X_raw) || is.null(scale_scores)) {
    dat <- read_split_half_data_cache(data_cache_path)
    X_raw <- dat$X_raw
    scale_scores <- dat$scale_scores
  }
  list(
    rep = task$rep,
    Xa_raw = X_raw[task$idx_a, , drop = FALSE],
    Xb_raw = X_raw[task$idx_b, , drop = FALSE],
    scale_scores_a = scale_scores[task$idx_a, , drop = FALSE],
    scale_scores_b = scale_scores[task$idx_b, , drop = FALSE]
  )
}

make_split_half_data_cache <- function(X_raw, scale_scores, prefix = "split_half_data_") {
  path <- tempfile(prefix, fileext = ".rds")
  saveRDS(
    list(X_raw = X_raw, scale_scores = scale_scores),
    file = path,
    compress = FALSE
  )
  path
}

matched_random_subset_baseline_worker <- function(task,
                                                  item_pool,
                                                  scale_key,
                                                  target_items,
                                                  target_method,
                                                  target_decomp_method,
                                                  scale_scores,
                                                  X_std_all,
                                                  uniform_weights,
                                                  worker_cfg) {
  mode <- task$subset_mode[[1]]
  b <- task$rep[[1]]
  iter <- ((match(mode, worker_cfg$baseline_subset_modes) - 1L) * as.integer(worker_cfg$baseline_reps)) + b
  set.seed(worker_cfg$seed + iter + match(mode, worker_cfg$baseline_subset_modes) * 100000L)
  subset_items <- sample_baseline_subset(
    item_pool = item_pool,
    scale_key = scale_key,
    target_items = target_items,
    mode = mode
  )
  fit_baseline <- tryCatch(
    fit_fixed_selection_from_prepared(
      X_std_all = X_std_all,
      selected_items = subset_items,
      item_weights = uniform_weights,
      decomp_method = target_decomp_method,
      spectrum_rank = worker_cfg$spectrum_rank
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_baseline)) {
    return(tibble::tibble(
      target_method = target_method,
      subset_mode = mode,
      rep = b,
      n_items_selected = length(target_items),
      pc1_explained_variance_ratio = NA_real_,
      pc2_explained_variance_ratio = NA_real_,
      pc12_explained_variance_ratio = NA_real_,
      spectral_effective_rank = NA_real_,
      cv_r2_scale_prediction_mean = NA_real_,
      cv_r2_scale_prediction_median = NA_real_
    ))
  }
  
  scale_prediction_tbl <- scale_prediction_for_method(
    method = target_method,
    pc_scores_2d = fit_baseline$pc_scores_2d,
    scale_scores = scale_scores,
    cfg = worker_cfg,
    show_progress = FALSE,
    parallel = FALSE
  )
  
  tibble::tibble(
    target_method = target_method,
    subset_mode = mode,
    rep = b,
    n_items_selected = length(fit_baseline$selected_items),
    pc1_explained_variance_ratio = fit_baseline$explained_variance_ratio[[1]],
    pc2_explained_variance_ratio = fit_baseline$explained_variance_ratio[[2]],
    pc12_explained_variance_ratio = cumulative_explained_variance_ratio(
      fit_baseline$explained_variance_ratio,
      2L
    ),
    spectral_effective_rank = fit_baseline$spectral_effective_rank,
    cv_r2_scale_prediction_mean = mean_or_na(scale_prediction_tbl$cv_r2),
    cv_r2_scale_prediction_median = median_or_na(scale_prediction_tbl$cv_r2)
  )
}

split_half_full_pipeline_refit_stability <- function(X_raw,
                                                     scale_scores,
                                                     weighting_mode,
                                                     decomp_method,
                                                     cfg,
                                                     show_progress = TRUE,
                                                     parallel = TRUE) {
  limit_worker_threads()
  
  n <- nrow(X_raw)
  half_n <- resolve_split_half_n(n, cfg$stability_half_n)
  reps <- resolve_stability_reps(cfg, "full_refit")
  
  tasks <- build_split_half_tasks(
    X_raw = X_raw,
    scale_scores = scale_scores,
    reps = reps,
    half_n = half_n,
    seed = cfg$seed
  )
  
  if (identical(weighting_mode, "id_guided")) {
    optimisation_runs <- resolve_selection_stability_optimisation_runs(cfg)
    
    worker_cfg <- cfg
    worker_cfg$optimisation_multi_runs <- optimisation_runs
    worker_cfg$progress_callback <- NULL
    worker_cfg$disable_inner_parallel <- TRUE
    
    use_cache <- isTRUE(parallel) &&
      future_worker_count() > 1L &&
      isTRUE(cfg$split_half_data_cache %||% TRUE)
    data_cache_path <- if (use_cache) make_split_half_data_cache(X_raw, scale_scores) else NULL
    if (!is.null(data_cache_path)) on.exit(unlink(data_cache_path), add = TRUE)
    
    if (isTRUE(show_progress)) {
      messagef(
        "  id-guided split-half full-refit: one FUTURE_LAPPLY split batch; workers=%d; splits=%d; optimisation_multi_runs_per_half=%d; optimiser_run_ticks=%d; data_cache=%s",
        future_worker_count(), length(tasks), optimisation_runs, 2L * optimisation_runs, !is.null(data_cache_path)
      )
    }
    
    rows <- map_with_progress(
      X = tasks,
      FUN = split_half_id_guided_full_refit_worker,
      cfg = cfg,
      decomp_method = decomp_method,
      worker_cfg = worker_cfg,
      X_raw = if (is.null(data_cache_path)) X_raw else NULL,
      scale_scores = if (is.null(data_cache_path)) scale_scores else NULL,
      data_cache_path = data_cache_path,
      show_progress = show_progress,
      parallel = isTRUE(parallel),
      future.seed = TRUE,
      future.scheduling = cfg$future_scheduling %||% 1L,
      progress_steps = length(tasks) * 2L * optimisation_runs,
      progress_arg_name = "progress_fun",
      progress_at_end = FALSE
    )
    
    return(dplyr::bind_rows(rows))
  }
  
  worker_cfg <- cfg
  worker_cfg$disable_inner_parallel <- TRUE
  worker_cfg$progress_callback <- NULL
  
  use_cache <- isTRUE(parallel) &&
    future_worker_count() > 1L &&
    isTRUE(cfg$split_half_data_cache %||% TRUE)
  data_cache_path <- if (use_cache) make_split_half_data_cache(X_raw, scale_scores) else NULL
  if (!is.null(data_cache_path)) on.exit(unlink(data_cache_path), add = TRUE)
  
  rows <- map_with_progress(
    X = tasks,
    FUN = split_half_full_pipeline_refit_worker,
    cfg = cfg,
    weighting_mode = weighting_mode,
    decomp_method = decomp_method,
    worker_cfg = worker_cfg,
    X_raw = if (is.null(data_cache_path)) X_raw else NULL,
    scale_scores = if (is.null(data_cache_path)) scale_scores else NULL,
    data_cache_path = data_cache_path,
    show_progress = show_progress,
    parallel = parallel,
    future.scheduling = cfg$future_scheduling %||% 1L
  )
  
  dplyr::bind_rows(rows)
}

split_half_fixed_selection_stability <- function(X_raw,
                                                 scale_scores,
                                                 fit_reference,
                                                 cfg,
                                                 show_progress = TRUE,
                                                 parallel = TRUE) {
  limit_worker_threads()
  n <- nrow(X_raw)
  half_n <- resolve_split_half_n(n, cfg$stability_half_n)
  reps <- resolve_stability_reps(cfg, "fixed_selection")
  fixed_weights <- fit_reference$weights
  fixed_items <- fit_reference$selected_items
  cfg_worker <- cfg
  tasks <- build_split_half_tasks(
    X_raw = X_raw,
    scale_scores = scale_scores,
    reps = reps,
    half_n = half_n,
    seed = cfg$seed
  )
  
  use_cache <- isTRUE(parallel) &&
    future_worker_count() > 1L &&
    isTRUE(cfg$split_half_data_cache %||% TRUE)
  data_cache_path <- if (use_cache) make_split_half_data_cache(X_raw, scale_scores) else NULL
  if (!is.null(data_cache_path)) on.exit(unlink(data_cache_path), add = TRUE)
  
  rows <- map_with_progress(
    X = tasks,
    FUN = split_half_fixed_selection_worker,
    cfg = cfg,
    fixed_items = fixed_items,
    fixed_weights = fixed_weights,
    decomp_method = fit_reference$decomp_method,
    worker_cfg = cfg_worker,
    X_raw = if (is.null(data_cache_path)) X_raw else NULL,
    scale_scores = if (is.null(data_cache_path)) scale_scores else NULL,
    data_cache_path = data_cache_path,
    show_progress = show_progress,
    parallel = parallel,
    future.scheduling = cfg$future_scheduling %||% 1L
  )
  dplyr::bind_rows(rows)
}

scale_prediction_for_method <- function(method,
                                        pc_scores_2d,
                                        scale_scores,
                                        cfg,
                                        show_progress = TRUE,
                                        parallel = TRUE) {
  rows <- map_with_progress(
    X = as.list(seq_along(scale_scores)),
    FUN = scale_prediction_worker,
    cfg = cfg,
    method = method,
    pc_scores_2d = pc_scores_2d,
    scale_scores = scale_scores,
    cv_folds = cfg$cv_folds,
    seed = cfg$seed,
    show_progress = show_progress,
    parallel = parallel
  )
  dplyr::bind_rows(rows)
}

leave_one_scale_out_prediction_for_method <- function(fit_reference,
                                                      X_raw,
                                                      scale_scores,
                                                      scale_key,
                                                      cfg,
                                                      prep = NULL,
                                                      show_progress = TRUE,
                                                      parallel = TRUE) {
  if (is.null(prep)) {
    prep <- prepare_standardised_matrix(X_raw)
  }
  available_items <- colnames(prep$std)
  uniform_weights <- setNames(rep(1, length(available_items)), available_items)
  scale_items_by_name <- split(scale_key$item_id, scale_key$scale_name)
  
  rows <- map_with_progress(
    X = as.list(seq_along(scale_scores)),
    FUN = leave_one_scale_out_prediction_worker,
    cfg = cfg,
    method = fit_reference$method,
    weighting_mode = fit_reference$weighting_mode,
    decomp_method = fit_reference$decomp_method,
    reference_selected_items = fit_reference$selected_items,
    reference_weights = fit_reference$weights,
    scale_scores = scale_scores,
    scale_items_by_name = scale_items_by_name,
    X_std_all = prep$std,
    available_items = available_items,
    uniform_weights = uniform_weights,
    spectrum_rank = cfg$spectrum_rank,
    cv_folds = cfg$cv_folds,
    seed = cfg$seed,
    show_progress = show_progress,
    parallel = parallel
  )
  
  dplyr::bind_rows(rows)
}

sample_baseline_subset <- function(item_pool,
                                   scale_key,
                                   target_items,
                                   mode) {
  target_n <- length(target_items)
  if (target_n < 2L) {
    stop("Baseline subset target size must be at least 2.")
  }
  
  if (identical(mode, "unstratified")) {
    return(sample(item_pool, target_n, replace = FALSE))
  }
  
  if (!identical(mode, "stratified_by_scale")) {
    stop("Unknown baseline subset mode: ", mode)
  }
  
  scale_map <- scale_key |>
    dplyr::distinct(item_id, scale_name) |>
    dplyr::filter(item_id %in% item_pool)
  target_counts <- scale_map |>
    dplyr::filter(item_id %in% target_items) |>
    dplyr::count(scale_name, name = "n_keep")
  
  sampled <- unlist(lapply(seq_len(nrow(target_counts)), function(i) {
    pool_i <- scale_map$item_id[scale_map$scale_name == target_counts$scale_name[[i]]]
    sample(pool_i, target_counts$n_keep[[i]], replace = FALSE)
  }), use.names = FALSE)
  
  sampled <- unique(sampled)
  if (length(sampled) < target_n) {
    remainder <- setdiff(item_pool, sampled)
    sampled <- c(sampled, sample(remainder, target_n - length(sampled), replace = FALSE))
  }
  sampled[seq_len(target_n)]
}

matched_random_subset_baseline <- function(X_raw,
                                           scale_scores,
                                           scale_key,
                                           target_fit,
                                           cfg,
                                           prep = NULL,
                                           show_progress = TRUE,
                                           parallel = TRUE) {
  limit_worker_threads()
  if (is.null(prep)) {
    prep <- prepare_standardised_matrix(X_raw)
  }
  cfg_worker <- cfg
  item_pool <- colnames(prep$std)
  uniform_weights <- setNames(rep(1, length(item_pool)), item_pool)
  target_items <- intersect(target_fit$selected_items, item_pool)
  target_n <- length(target_items)
  if (target_n < 2L) {
    return(tibble::tibble())
  }
  
  task_grid <- expand.grid(
    subset_mode = cfg$baseline_subset_modes,
    rep = seq_len(as.integer(cfg$baseline_reps)),
    stringsAsFactors = FALSE
  )
  
  rows <- map_with_progress(
    X = split(task_grid, seq_len(nrow(task_grid))),
    FUN = matched_random_subset_baseline_worker,
    cfg = cfg,
    item_pool = item_pool,
    scale_key = scale_key,
    target_items = target_items,
    target_method = target_fit$method,
    target_decomp_method = target_fit$decomp_method,
    scale_scores = scale_scores,
    X_std_all = prep$std,
    uniform_weights = uniform_weights,
    worker_cfg = cfg_worker,
    show_progress = show_progress,
    parallel = parallel
  )
  
  dplyr::bind_rows(rows)
}

knn_index_matrix <- function(D, k) {
  M <- as.matrix(D)
  n <- nrow(M)
  if (n < 2L) return(matrix(integer(0), nrow = n, ncol = 0L))
  diag(M) <- Inf
  k <- min(as.integer(k), n - 1L)
  t(apply(M, 1L, function(row) order(row, decreasing = FALSE)[seq_len(k)]))
}

mean_knn_overlap <- function(D_reference, D_map, k) {
  ref_nn <- knn_index_matrix(D_reference, k)
  map_nn <- knn_index_matrix(D_map, k)
  if (!ncol(ref_nn) || !ncol(map_nn)) return(NA_real_)
  mean(vapply(seq_len(nrow(ref_nn)), function(i) {
    length(intersect(ref_nn[i, ], map_nn[i, ])) / ncol(ref_nn)
  }, numeric(1)), na.rm = TRUE)
}

distance_preservation_for_method <- function(fit,
                                             X_raw,
                                             cfg,
                                             show_progress = TRUE) {
  n <- nrow(X_raw)
  n_eval <- suppressWarnings(as.integer(cfg$distance_eval_n))
  if (!length(n_eval) || is.na(n_eval) || !is.finite(n_eval) || n_eval <= 1L) {
    n_eval <- min(n, 1200L)
  }
  n_eval <- min(n_eval, n)
  set.seed(cfg$seed + match(fit$method, sort(unique(fit$method))) + 9000L)
  idx <- if (n > n_eval) sample(seq_len(n), n_eval) else seq_len(n)
  
  selected_items <- intersect(fit$selected_items, colnames(X_raw))
  if (length(selected_items) < 2L) {
    return(tibble::tibble(
      method = fit$method,
      n_eval = n_eval,
      k_neighbors = NA_integer_,
      n_items_selected = length(selected_items),
      gower_twonn_id = NA_real_,
      map_distance_pearson = NA_real_,
      map_distance_spearman = NA_real_,
      neighbourhood_overlap_gower_to_map = NA_real_,
      neighbourhood_overlap_map_to_gower = NA_real_,
      pcoa2_positive_eigen_share = NA_real_
    ))
  }
  
  weights <- fit$weights[selected_items]
  weights[!is.finite(weights) | weights <= 0] <- 1
  prep <- prep_ord_gower(X_raw[idx, selected_items, drop = FALSE])
  D_gower <- as.matrix(gower_dist(prep$X, prep$type, weights = weights))
  U <- fit$pc_scores_2d[idx, , drop = FALSE]
  D_map <- as.matrix(stats::dist(U))
  upper <- upper.tri(D_gower)
  g <- D_gower[upper]
  m <- D_map[upper]
  ok <- is.finite(g) & is.finite(m)
  k_neighbors <- suppressWarnings(as.integer(cfg$distance_eval_k))
  if (!length(k_neighbors) || is.na(k_neighbors) || !is.finite(k_neighbors) || k_neighbors < 1L) {
    k_neighbors <- 15L
  }
  k_neighbors <- min(k_neighbors, n_eval - 1L)
  
  eig <- tryCatch(
    stats::cmdscale(stats::as.dist(D_gower), k = 2L, eig = TRUE)$eig,
    error = function(e) numeric(0)
  )
  eig_pos <- eig[eig > 0 & is.finite(eig)]
  pcoa_share <- if (length(eig_pos)) sum(head(eig_pos, 2L)) / sum(eig_pos) else NA_real_
  
  tibble::tibble(
    method = fit$method,
    n_eval = n_eval,
    k_neighbors = k_neighbors,
    n_items_selected = length(selected_items),
    gower_twonn_id = twonn_id_from_dist(stats::as.dist(D_gower)),
    map_distance_pearson = if (sum(ok) > 2L) suppressWarnings(stats::cor(g[ok], m[ok], method = "pearson")) else NA_real_,
    map_distance_spearman = if (sum(ok) > 2L) suppressWarnings(stats::cor(g[ok], m[ok], method = "spearman")) else NA_real_,
    neighbourhood_overlap_gower_to_map = mean_knn_overlap(D_gower, D_map, k_neighbors),
    neighbourhood_overlap_map_to_gower = mean_knn_overlap(D_map, D_gower, k_neighbors),
    pcoa2_positive_eigen_share = pcoa_share
  )
}

format_method_label <- function(method) {
  dplyr::recode(
    method,
    pca_uniform = "All-item\nPCA",
    robpca_uniform = "All-item\nrobust PCA",
    pca_id_guided = "Selected-item\nPCA",
    robpca_id_guided = "Selected-item\nrobust PCA",
    .default = gsub("_", "\n", method, fixed = TRUE)
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

method_sensitivity_plot_theme <- function(theme_pub_fn,
                                          base_size = 11,
                                          legend_position = "right",
                                          y_grid = FALSE) {
  plot_theme <- theme_pub_fn(base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        size = base_size + 2
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0,
        size = base_size - 1,
        lineheight = 1.05
      ),
      plot.caption = ggplot2::element_text(size = base_size - 2),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(3, "pt"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      strip.text = ggplot2::element_text(face = "bold", colour = "black"),
      legend.position = legend_position,
      plot.margin = grid::unit(c(5.5, 5.5, 5.5, 5.5), "pt")
    )
  
  if (isTRUE(y_grid)) {
    plot_theme <- plot_theme +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_line(
          colour = "grey85",
          linetype = "dotted",
          linewidth = 0.25
        ),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  plot_theme
}

save_method_sensitivity_plot <- function(setup_helpers,
                                         name,
                                         plot,
                                         width = 8,
                                         height = 5,
                                         dpi = 300) {
  setup_helpers$save_plot_gg(
    name = name,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

add_pareto_labels <- function(data, y_col = "pc12_explained_variance_ratio") {
  ggrepel::geom_text_repel(
    data = data,
    ggplot2::aes(label = method_label, y = .data[[y_col]]),
    size = 3,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE
  )
}

load_fit_bundle <- function(path, weighting_mode, decomp_method) {
  bundle <- readRDS(path)
  if (is.null(bundle$item_component_correlations) || !nrow(bundle$item_component_correlations)) {
    stop("Bundle is missing item_component_correlations: ", path)
  }
  pc_scores_2d <- as.matrix(bundle$pc_scores_2d)
  if (ncol(pc_scores_2d) != 2L) stop("Expected two score columns in bundle: ", path)
  colnames(pc_scores_2d) <- c("u1", "u2")
  item_component_correlations <- as.matrix(bundle$item_component_correlations)
  colnames(item_component_correlations) <- c("u1", "u2")
  spectrum <- as.numeric(bundle$spectrum)
  explained_variance_ratio <- as.numeric(bundle$explained_variance_ratio)
  if (!length(explained_variance_ratio) && length(spectrum)) {
    explained_variance_ratio <- spectrum / pmax(sum(spectrum), 1e-12)
  }
  stability_path <- file.path(dirname(path), "gower_survivor_stability.csv")
  selection_stability <- NULL
  if (file.exists(stability_path)) {
    selection_stability <- suppressMessages(readr::read_csv2(
      stability_path,
      show_col_types = FALSE,
      progress = FALSE
    ))
    if ("var" %in% names(selection_stability)) {
      selection_stability <- dplyr::rename(selection_stability, item_id = var)
    }
  }
  list(
    participant_id = as.character(bundle$participant_id),
    selected_items = as.character(bundle$selected_items),
    weights = bundle$weights,
    pc_scores_2d = pc_scores_2d,
    spectrum = spectrum,
    explained_variance_ratio = explained_variance_ratio,
    spectral_effective_rank = spectral_effective_rank(spectrum),
    item_component_correlations = item_component_correlations,
    selection_stability = selection_stability,
    weighting_mode = weighting_mode,
    decomp_method = decomp_method,
    method = method_label(weighting_mode, decomp_method)
  )
}

refresh_spectral_summary <- function(fit_reference, X_raw, cfg, prep = NULL) {
  if (is.null(prep)) {
    prep <- prepare_standardised_matrix(X_raw)
  }
  spectral_fit <- tryCatch(
    fit_fixed_selection_from_prepared(
      X_std_all = prep$std,
      selected_items = fit_reference$selected_items,
      item_weights = fit_reference$weights,
      decomp_method = fit_reference$decomp_method,
      spectrum_rank = cfg$spectrum_rank
    ),
    error = function(e) NULL
  )
  
  if (is.null(spectral_fit)) return(fit_reference)
  
  fit_reference$spectrum <- spectral_fit$spectrum
  fit_reference$explained_variance_ratio <- spectral_fit$explained_variance_ratio
  fit_reference$spectral_effective_rank <- spectral_fit$spectral_effective_rank
  fit_reference
}

run_method_sensitivity_summary <- function(cfg = default_method_sensitivity_cfg) {
  cfg <- utils::modifyList(default_method_sensitivity_cfg, cfg)
  root_dir <- cfg$root_dir %||% infer_root_dir()
  cfg$data_path <- normalise_path(cfg$data_path, root_dir, must_work = TRUE)
  cfg$qbank_path <- normalise_path(cfg$qbank_path, root_dir, must_work = TRUE)
  cfg$scale_spec_path <- normalise_path(cfg$scale_spec_path, root_dir, must_work = TRUE)
  cfg$setup_script <- normalise_path(cfg$setup_script, root_dir, must_work = TRUE)
  cfg$out_dir <- normalise_path(cfg$out_dir, root_dir, must_work = FALSE)
  
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(cfg$seed)
  setup_helpers <- load_setup_helpers(cfg, root_dir)
  limit_worker_threads()
  messagef("[method_sensitivity_summary] workers=%d", setup_helpers$nworkers)
  
  messagef("Reading data from %s", cfg$data_path)
  raw_df <- readr::read_csv2(cfg$data_path, show_col_types = FALSE)
  if (!"participant_id" %in% names(raw_df)) {
    stop("Expected a participant_id column in ", cfg$data_path)
  }
  
  scale_key <- build_scale_key(cfg$scale_spec_path)
  qbank <- read_qbank(cfg$qbank_path)
  item_cols <- intersect(scale_key$item_id, names(raw_df))
  scale_key <- scale_key |>
    dplyr::filter(item_id %in% item_cols)
  
  method_specs <- tibble::tribble(
    ~weighting_mode, ~decomp_method, ~run_dir,
    "uniform", "pca", "pca_uniform",
    "uniform", "robpca", "robpca_uniform",
    "id_guided", "pca", "pca_id_guided",
    "id_guided", "robpca", "robpca_id_guided"
  )
  
  fits <- vector("list", nrow(method_specs))
  names(fits) <- vapply(seq_len(nrow(method_specs)), function(i) {
    method_label(method_specs$weighting_mode[[i]], method_specs$decomp_method[[i]])
  }, character(1))
  
  for (i in seq_len(nrow(method_specs))) {
    spec <- method_specs[i, , drop = FALSE]
    bundle_path <- file.path(cfg$out_dir, spec$run_dir[[1]], "method_sensitivity_fit_bundle.rds")
    if (!file.exists(bundle_path)) {
      stop("Missing fit bundle: ", bundle_path)
    }
    label <- method_label(spec$weighting_mode[[1]], spec$decomp_method[[1]])
    messagef("Loading %s", label)
    fits[[label]] <- load_fit_bundle(bundle_path, spec$weighting_mode[[1]], spec$decomp_method[[1]])
  }
  
  participant_id <- fits[[1L]]$participant_id
  for (nm in names(fits)) {
    if (!identical(fits[[nm]]$participant_id, participant_id)) {
      stop("Participant IDs differ across runs; cannot compare fits directly.")
    }
  }
  
  X_raw <- raw_df |>
    dplyr::select(tidyselect::all_of(item_cols)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ {
      x <- as.numeric(.x)
      x[x == 0] <- NA_real_
      x
    })) |>
    as.data.frame(check.names = FALSE)
  rownames(X_raw) <- make.unique(as.character(raw_df$participant_id))
  
  missing_ids <- setdiff(participant_id, rownames(X_raw))
  if (length(missing_ids)) {
    stop("Some fit participants were not found in raw data: ", paste(head(missing_ids, 5L), collapse = ", "))
  }
  X_raw <- X_raw[participant_id, item_cols, drop = FALSE]
  scale_scores <- score_scales(X_raw, scale_key, min_completion_prop = cfg$min_scale_completion_prop)
  prep_full <- prepare_standardised_matrix(X_raw)
  
  quality_df <- tibble::tibble(
    step = c("raw_rows", "dimension_psychometric_rows"),
    n = c(nrow(raw_df), length(participant_id))
  )
  write_csv(quality_df, file.path(cfg$out_dir, "qc_row_counts.csv"))
  
  for (nm in names(fits)) {
    fits[[nm]] <- refresh_spectral_summary(fits[[nm]], X_raw, cfg, prep = prep_full)
  }
  
  reference_method <- if ("robpca_id_guided" %in% names(fits)) "robpca_id_guided" else names(fits)[[1L]]
  reference_item_component_correlations <- fits[[reference_method]]$item_component_correlations
  for (nm in names(fits)) {
    fits[[nm]] <- align_fit_to_reference(fits[[nm]], reference_item_component_correlations)
  }
  
  eigenspectrum_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    tibble::tibble(
      method = fit$method,
      component = seq_along(fit$spectrum),
      eigenvalue = fit$spectrum,
      explained_variance_ratio = fit$explained_variance_ratio
    )
  }))
  write_csv(eigenspectrum_tbl, file.path(cfg$out_dir, "method_sensitivity_eigenspectrum.csv"))
  
  weights_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    base <- tibble::tibble(
      method = fit$method,
      item_id = names(fit$weights),
      weight = as.numeric(fit$weights),
      selected = names(fit$weights) %in% fit$selected_items
    )
    if (!is.null(fit$selection_stability) && nrow(fit$selection_stability)) {
      base <- dplyr::left_join(
        base,
        fit$selection_stability |>
          dplyr::select(
            item_id,
            selection_count = count,
            selection_prop = prop,
            cutpoint_count = count_selected,
            cutpoint_prop = prop_selected,
            selected_ref,
            active_ref = in_ref,
            weight_ref,
            weight_mean,
            weight_sd,
            weight_min,
            weight_max
          ),
        by = "item_id"
      )
    }
    base
  }))
  weights_tbl <- weights_tbl |>
    dplyr::left_join(
      scale_key |> dplyr::distinct(item_id, scale_name, original_16pf_factor, keyed_sign),
      by = "item_id"
    ) |>
    dplyr::left_join(qbank |> dplyr::select(item_id, item_text), by = "item_id")
  write_csv(weights_tbl, file.path(cfg$out_dir, "method_sensitivity_item_weights.csv"))
  
  selected_item_content_tbl <- weights_tbl |>
    dplyr::filter(selected, grepl("id_guided", method, fixed = TRUE)) |>
    dplyr::arrange(method, dplyr::desc(selection_prop), dplyr::desc(weight), item_id) |>
    dplyr::select(
      method,
      item_id,
      scale_name,
      original_16pf_factor,
      keyed_sign,
      weight,
      weight_mean,
      weight_sd,
      selection_prop,
      selection_count,
      cutpoint_prop,
      selected_ref,
      active_ref,
      item_text
    )
  write_csv(
    selected_item_content_tbl,
    file.path(cfg$out_dir, "method_sensitivity_selected_item_content.csv")
  )
  
  pc_scores_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    tibble::tibble(
      participant_id = participant_id,
      method = fit$method,
      u1 = fit$pc_scores_2d[, 1],
      u2 = fit$pc_scores_2d[, 2]
    )
  }))
  write_csv(pc_scores_tbl, file.path(cfg$out_dir, "method_sensitivity_pc_scores_2d.csv"))
  
  item_component_correlation_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    tibble::tibble(
      item_id = rownames(fit$item_component_correlations),
      method = fit$method,
      u1 = fit$item_component_correlations[, 1],
      u2 = fit$item_component_correlations[, 2]
    )
  }))
  item_component_correlation_tbl <- item_component_correlation_tbl |>
    dplyr::left_join(qbank |> dplyr::select(item_id, item_text), by = "item_id") |>
    dplyr::left_join(scale_key |> dplyr::distinct(item_id, scale_name, original_16pf_factor), by = "item_id")
  write_csv(
    item_component_correlation_tbl,
    file.path(cfg$out_dir, "method_sensitivity_item_component_correlations.csv")
  )
  
  scale_prediction_tbl <- dplyr::bind_rows(lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    messagef("Scale score prediction for %s", fit$method)
    scale_prediction_for_method(
      fit$method,
      fit$pc_scores_2d,
      scale_scores,
      cfg,
      parallel = isTRUE(cfg$post_parallel_diagnostics)
    )
  }))
  write_csv(scale_prediction_tbl, file.path(cfg$out_dir, "method_sensitivity_scale_prediction.csv"))
  write_csv(scale_prediction_tbl, file.path(cfg$out_dir, "method_sensitivity_fixed_map_scale_readout.csv"))
  
  leave_one_scale_out_prediction_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    messagef("Leave-one-scale-out prediction for %s", fit$method)
    leave_one_scale_out_prediction_for_method(
      fit_reference = fit,
      X_raw = X_raw,
      scale_scores = scale_scores,
      scale_key = scale_key,
      cfg = cfg,
      prep = prep_full,
      parallel = isTRUE(cfg$post_parallel_diagnostics)
    )
  }))
  write_csv(
    leave_one_scale_out_prediction_tbl,
    file.path(cfg$out_dir, "method_sensitivity_leave_one_scale_out_prediction.csv")
  )
  write_csv(
    leave_one_scale_out_prediction_tbl,
    file.path(cfg$out_dir, "method_sensitivity_leave_one_scale_out_readout.csv")
  )
  
  scale_regression_vectors_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    scale_regression_vectors(scale_scores, fit$pc_scores_2d) |>
      dplyr::mutate(method = fit$method, .before = 1)
  }))
  write_csv(
    scale_regression_vectors_tbl,
    file.path(cfg$out_dir, "method_sensitivity_scale_regression_vectors.csv")
  )
  
  distance_preservation_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    messagef("Gower-space distance preservation for %s", fit$method)
    distance_preservation_for_method(
      fit = fit,
      X_raw = X_raw,
      cfg = cfg
    )
  }))
  write_csv(
    distance_preservation_tbl,
    file.path(cfg$out_dir, "method_sensitivity_distance_preservation.csv")
  )
  
  full_pipeline_stability_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(method_specs)), function(i) {
    weighting_mode <- method_specs$weighting_mode[[i]]
    decomp_method <- method_specs$decomp_method[[i]]
    label <- method_label(weighting_mode, decomp_method)
    messagef("Split-half full-refit stability for %s", label)
    if (identical(weighting_mode, "id_guided")) {
      messagef("  full-refit parallel=%s; workers=%d; reps=%d; final_optimisation_multi_runs=%d; selection_stability_optimisation_multi_runs=%d", isTRUE(cfg$split_half_full_refit_parallel), future_worker_count(), resolve_stability_reps(cfg, "full_refit"), as.integer(cfg$optimisation_multi_runs), resolve_selection_stability_optimisation_runs(cfg))
    } else {
      messagef("  full-refit parallel=%s; workers=%d; reps=%d; optimisation_multi_runs=%d", isTRUE(cfg$split_half_full_refit_parallel), future_worker_count(), resolve_stability_reps(cfg, "full_refit"), as.integer(cfg$optimisation_multi_runs))
    }
    split_half_full_pipeline_refit_stability(
      X_raw = X_raw,
      scale_scores = scale_scores,
      weighting_mode = weighting_mode,
      decomp_method = decomp_method,
      cfg = cfg,
      parallel = isTRUE(cfg$split_half_full_refit_parallel)
    ) |>
      dplyr::mutate(method = label, .before = 1)
  }))
  write_csv(
    full_pipeline_stability_tbl,
    file.path(cfg$out_dir, "method_sensitivity_split_half_full_refit_stability.csv")
  )
  
  fixed_selection_stability_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    messagef("Split-half fixed-selection stability for %s", fit$method)
    split_half_fixed_selection_stability(
      X_raw = X_raw,
      scale_scores = scale_scores,
      fit_reference = fit,
      cfg = cfg,
      parallel = isTRUE(cfg$post_parallel_diagnostics)
    ) |>
      dplyr::mutate(method = fit$method, .before = 1)
  }))
  write_csv(
    fixed_selection_stability_tbl,
    file.path(cfg$out_dir, "method_sensitivity_split_half_fixed_selection_stability.csv")
  )
  
  id_guided_methods <- vapply(
    fits[vapply(fits, function(fit) identical(fit$weighting_mode, "id_guided"), logical(1))],
    function(fit) fit$method,
    character(1)
  )
  selection_stability_tbl <- full_pipeline_stability_tbl |>
    dplyr::filter(method %in% id_guided_methods) |>
    dplyr::select(
      method,
      rep,
      selection_jaccard,
      aligned_weight_correlation,
      n_selected_a,
      n_selected_b,
      n_selected_intersection
    )
  write_csv(selection_stability_tbl, file.path(cfg$out_dir, "method_sensitivity_selection_stability.csv"))
  
  baseline_tbl <- dplyr::bind_rows(lapply(id_guided_methods, function(method_name) {
    fit <- fits[[method_name]]
    messagef("Matched random-subset baseline for %s", fit$method)
    matched_random_subset_baseline(
      X_raw = X_raw,
      scale_scores = scale_scores,
      scale_key = scale_key,
      target_fit = fit,
      cfg = cfg,
      prep = prep_full,
      parallel = isTRUE(cfg$post_parallel_diagnostics)
    )
  }))
  write_csv(
    baseline_tbl,
    file.path(cfg$out_dir, "method_sensitivity_matched_random_subset_baseline.csv")
  )
  
  baseline_summary_tbl <- baseline_tbl |>
    dplyr::group_by(target_method, subset_mode) |>
    dplyr::summarise(
      n_reps = dplyr::n(),
      pc12_mean = mean_or_na(pc12_explained_variance_ratio),
      pc12_q05 = quantile_or_na(pc12_explained_variance_ratio, probs = 0.05),
      pc12_q95 = quantile_or_na(pc12_explained_variance_ratio, probs = 0.95),
      cv_r2_scale_prediction_mean = mean_or_na(cv_r2_scale_prediction_mean),
      cv_r2_scale_prediction_q05 = quantile_or_na(cv_r2_scale_prediction_mean, probs = 0.05),
      cv_r2_scale_prediction_q95 = quantile_or_na(cv_r2_scale_prediction_mean, probs = 0.95),
      .groups = "drop"
    )
  write_csv(
    baseline_summary_tbl,
    file.path(cfg$out_dir, "method_sensitivity_matched_random_subset_baseline_summary.csv")
  )
  
  summary_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    scale_prediction_rows <- scale_prediction_tbl |> dplyr::filter(method == fit$method)
    leave_one_scale_out_rows <- leave_one_scale_out_prediction_tbl |> dplyr::filter(method == fit$method)
    stability_full <- full_pipeline_stability_tbl |> dplyr::filter(method == fit$method)
    stability_fixed <- fixed_selection_stability_tbl |> dplyr::filter(method == fit$method)
    selection_rows <- selection_stability_tbl |> dplyr::filter(method == fit$method)
    distance_rows <- distance_preservation_tbl |> dplyr::filter(method == fit$method)
    
    tibble::tibble(
      method = fit$method,
      weighting_mode = fit$weighting_mode,
      decomposition = fit$decomp_method,
      n_items_selected = length(fit$selected_items),
      pc1_explained_variance_ratio = fit$explained_variance_ratio[[1]],
      pc2_explained_variance_ratio = fit$explained_variance_ratio[[2]],
      pc12_explained_variance_ratio = cumulative_explained_variance_ratio(fit$explained_variance_ratio, 2L),
      pc5_explained_variance_ratio = cumulative_explained_variance_ratio(fit$explained_variance_ratio, 5L),
      pc10_explained_variance_ratio = cumulative_explained_variance_ratio(fit$explained_variance_ratio, 10L),
      spectral_effective_rank = fit$spectral_effective_rank,
      cv_r2_scale_prediction_mean = mean_or_na(scale_prediction_rows$cv_r2),
      cv_r2_scale_prediction_median = median_or_na(scale_prediction_rows$cv_r2),
      cv_r2_leave_one_scale_out_mean = mean_or_na(leave_one_scale_out_rows$cv_r2),
      cv_r2_leave_one_scale_out_median = median_or_na(leave_one_scale_out_rows$cv_r2),
      gower_map_distance_spearman = mean_or_na(distance_rows$map_distance_spearman),
      gower_map_distance_pearson = mean_or_na(distance_rows$map_distance_pearson),
      gower_to_map_neighbourhood_overlap = mean_or_na(distance_rows$neighbourhood_overlap_gower_to_map),
      map_to_gower_neighbourhood_overlap = mean_or_na(distance_rows$neighbourhood_overlap_map_to_gower),
      pcoa2_positive_eigen_share = mean_or_na(distance_rows$pcoa2_positive_eigen_share),
      full_refit_axis_corr_mean = mean_or_na(stability_full$axis_corr_mean),
      full_refit_item_rmse_mean = mean_or_na(stability_full$item_rmse),
      full_refit_scale_regression_vector_cosine_mean = mean_or_na(stability_full$scale_regression_vector_cosine),
      fixed_selection_axis_corr_mean = mean_or_na(stability_fixed$axis_corr_mean),
      fixed_selection_item_rmse_mean = mean_or_na(stability_fixed$item_rmse),
      fixed_selection_scale_regression_vector_cosine_mean = mean_or_na(stability_fixed$scale_regression_vector_cosine),
      selection_jaccard_mean = mean_or_na(selection_rows$selection_jaccard),
      selection_aligned_weight_correlation_mean = mean_or_na(selection_rows$aligned_weight_correlation)
    )
  }))
  write_csv(summary_tbl, file.path(cfg$out_dir, "method_sensitivity_summary.csv"))
  
  if (isTRUE(cfg$save_plots)) {
    method_levels <- summary_tbl$method
    method_labels <- stats::setNames(format_method_label(method_levels), method_levels)
    method_palette <- stats::setNames(
      unname(setup_helpers$cluster_colours(method_levels)),
      method_levels
    )
    decomp_shape_values <- c("PCA" = 21, "Robust PCA" = 24)
    decomp_names <- c("pca" = "PCA", "robpca" = "Robust PCA")
    
    map_method <- if ("robpca_id_guided" %in% names(fits)) "robpca_id_guided" else {
      id_methods <- names(fits)[grepl("id_guided", names(fits), fixed = TRUE)]
      if (length(id_methods)) id_methods[[1L]] else names(fits)[[1L]]
    }
    map_fit <- fits[[map_method]]
    map_scores <- tibble::tibble(
      participant_id = participant_id,
      u1 = map_fit$pc_scores_2d[, 1],
      u2 = map_fit$pc_scores_2d[, 2]
    )
    padded_limits <- function(x, probs = c(0.001, 0.999), pad = 0.18) {
      q <- as.numeric(stats::quantile(x, probs, na.rm = TRUE))
      span <- diff(q)
      if (!is.finite(span) || span <= 0) span <- diff(range(x, na.rm = TRUE))
      if (!is.finite(span) || span <= 0) span <- 1
      q + c(-pad, pad) * span
    }
    x_limits <- padded_limits(map_scores$u1)
    y_limits <- padded_limits(map_scores$u2)
    axis_prefix <- if (identical(map_fit$decomp_method, "robpca")) "RC" else "PC"
    vector_base <- scale_regression_vectors_tbl |>
      dplyr::filter(method == map_fit$method) |>
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
        label = gsub("_", " ", scale_name, fixed = TRUE)
      )
    map_panel_data <- dplyr::bind_rows(
      map_scores |> dplyr::mutate(panel = "A. Respondent density"),
      map_scores |> dplyr::mutate(panel = "B. Scale-score gradients")
    )
    label_layer <- if (nrow(vector_plot_data) && requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_label_repel(
        data = vector_plot_data,
        aes(x = xend, y = yend, label = label),
        inherit.aes = FALSE,
        size = 2.15,
        min.segment.length = 0,
        box.padding = 0.28,
        point.padding = 0.12,
        force = 4,
        force_pull = 0.15,
        max.overlaps = Inf,
        label.size = 0.12,
        label.padding = grid::unit(0.08, "lines"),
        fill = "white",
        seed = cfg$seed
      )
    } else {
      geom_text(
        data = vector_plot_data,
        aes(x = xend, y = yend, label = label),
        inherit.aes = FALSE,
        size = 2.05,
        check_overlap = TRUE
      )
    }
    respondent_map_plot <- ggplot(map_panel_data, aes(u1, u2)) +
      stat_density_2d(
        aes(fill = after_stat(level)),
        geom = "polygon",
        contour = TRUE,
        bins = 18,
        colour = NA,
        alpha = 0.95
      ) +
      stat_density_2d(
        colour = "grey40",
        linewidth = 0.18,
        bins = 8
      ) +
      geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.25, colour = "grey45") +
      geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.25, colour = "grey45") +
      geom_segment(
        data = vector_plot_data,
        aes(x = x0, y = y0, xend = xend, yend = yend),
        inherit.aes = FALSE,
        linewidth = 0.42,
        arrow = grid::arrow(length = grid::unit(0.075, "in"))
      ) +
      label_layer +
      facet_wrap(~panel, nrow = 1) +
      scale_fill_gradient(low = "grey97", high = "grey42", guide = "none") +
      coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +
      labs(
        title = "Selected-item respondent map",
        x = sprintf("%s1 (%.1f%%)", axis_prefix, 100 * map_fit$explained_variance_ratio[[1]]),
        y = sprintf("%s2 (%.1f%%)", axis_prefix, 100 * map_fit$explained_variance_ratio[[2]])
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 10, legend_position = "none", y_grid = FALSE)
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_selected_item_respondent_map",
      respondent_map_plot,
      width = 9.4,
      height = 4.8
    )
    
    scree_plot <- eigenspectrum_tbl |>
      dplyr::filter(component <= cfg$spectrum_rank) |>
      dplyr::mutate(method = factor(method, levels = method_levels)) |>
      ggplot(aes(component, explained_variance_ratio, colour = method, group = method)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      scale_colour_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        name = NULL
      ) +
      scale_x_continuous(breaks = seq_len(cfg$spectrum_rank)) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(
        title = "Explained variance profile",
        x = "Component",
        y = "Explained variance"
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11, y_grid = TRUE)
    save_method_sensitivity_plot(setup_helpers, "FIG_method_sensitivity_scree", scree_plot, width = 8, height = 5)
    
    scale_prediction_plot <- scale_prediction_tbl |>
      dplyr::mutate(method = factor(method, levels = method_levels)) |>
      ggplot(aes(method, cv_r2, fill = method)) +
      geom_boxplot(width = 0.65, alpha = 0.85, outlier.size = 0.8) +
      scale_fill_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        guide = "none"
      ) +
      labs(
        title = "Fixed-map scale-score readout",
        x = NULL,
        y = "Five-fold readout R2"
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11, legend_position = "none", y_grid = TRUE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(lineheight = 0.95)
      ) +
      scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_scale_prediction",
      scale_prediction_plot,
      width = 8,
      height = 5
    )
    
    leave_one_scale_out_plot <- leave_one_scale_out_prediction_tbl |>
      dplyr::mutate(method = factor(method, levels = method_levels)) |>
      ggplot(aes(method, cv_r2, fill = method)) +
      geom_boxplot(width = 0.65, alpha = 0.85, outlier.size = 0.8) +
      scale_fill_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        guide = "none"
      ) +
      labs(
        title = "Leave-one-scale-out scale-score readout",
        x = NULL,
        y = "Five-fold readout R2"
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11, legend_position = "none", y_grid = TRUE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(lineheight = 0.95)
      ) +
      scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_leave_one_scale_out_prediction",
      leave_one_scale_out_plot,
      width = 8,
      height = 5
    )
    
    full_stability_long <- full_pipeline_stability_tbl |>
      dplyr::select(method, axis_corr_mean, item_rmse, scale_regression_vector_cosine) |>
      tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
      dplyr::mutate(
        method = factor(method, levels = method_levels),
        metric = factor(
          metric,
          levels = c("axis_corr_mean", "item_rmse", "scale_regression_vector_cosine"),
          labels = format_metric_label(c("axis_corr_mean", "item_rmse", "scale_regression_vector_cosine"))
        )
      )
    full_stability_plot <- full_stability_long |>
      ggplot(aes(method, value, fill = method)) +
      geom_boxplot(width = 0.65, alpha = 0.85, outlier.size = 0.8) +
      scale_fill_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        guide = "none"
      ) +
      facet_wrap(~metric, scales = "free_y") +
      labs(
        title = "Split-half full-refit stability",
        x = NULL,
        y = NULL
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11, legend_position = "none", y_grid = TRUE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(lineheight = 0.95)
      )
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_split_half_full_refit_stability",
      full_stability_plot,
      width = 10,
      height = 6
    )
    
    fixed_stability_long <- fixed_selection_stability_tbl |>
      dplyr::select(method, axis_corr_mean, item_rmse, scale_regression_vector_cosine) |>
      tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
      dplyr::mutate(
        method = factor(method, levels = method_levels),
        metric = factor(
          metric,
          levels = c("axis_corr_mean", "item_rmse", "scale_regression_vector_cosine"),
          labels = format_metric_label(c("axis_corr_mean", "item_rmse", "scale_regression_vector_cosine"))
        )
      )
    fixed_stability_plot <- fixed_stability_long |>
      ggplot(aes(method, value, fill = method)) +
      geom_boxplot(width = 0.65, alpha = 0.85, outlier.size = 0.8) +
      scale_fill_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        guide = "none"
      ) +
      facet_wrap(~metric, scales = "free_y") +
      labs(
        title = "Split-half fixed-selection stability",
        x = NULL,
        y = NULL
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11, legend_position = "none", y_grid = TRUE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(lineheight = 0.95)
      )
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_split_half_fixed_selection_stability",
      fixed_stability_plot,
      width = 10,
      height = 6
    )
    
    distance_plot <- distance_preservation_tbl |>
      dplyr::select(
        method,
        map_distance_spearman,
        neighbourhood_overlap_gower_to_map,
        neighbourhood_overlap_map_to_gower,
        pcoa2_positive_eigen_share
      ) |>
      tidyr::pivot_longer(-method, names_to = "metric", values_to = "value") |>
      dplyr::mutate(
        method = factor(method, levels = method_levels),
        metric = factor(
          metric,
          levels = c(
            "map_distance_spearman",
            "neighbourhood_overlap_gower_to_map",
            "neighbourhood_overlap_map_to_gower",
            "pcoa2_positive_eigen_share"
          ),
          labels = c(
            "Gower-map distance Spearman r",
            "Gower-to-map neighbour overlap",
            "Map-to-Gower neighbour overlap",
            "PCoA two-axis positive-eigen share"
          )
        )
      ) |>
      ggplot(aes(method, value, fill = method)) +
      geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
      facet_wrap(~metric, scales = "free_y", nrow = 2) +
      scale_fill_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        guide = "none"
      ) +
      labs(
        title = "Distance-geometry checks",
        x = NULL,
        y = NULL
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 10, legend_position = "none", y_grid = TRUE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(lineheight = 0.95))
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_distance_geometry",
      distance_plot,
      width = 10,
      height = 6
    )
    
    selection_plot <- selection_stability_tbl |>
      tidyr::pivot_longer(
        cols = c(selection_jaccard, aligned_weight_correlation),
        names_to = "metric",
        values_to = "value"
      ) |>
      dplyr::mutate(
        method = factor(method, levels = method_levels),
        metric = factor(
          metric,
          levels = c("selection_jaccard", "aligned_weight_correlation"),
          labels = format_metric_label(c("selection_jaccard", "aligned_weight_correlation"))
        )
      ) |>
      ggplot(aes(method, value, fill = method)) +
      geom_boxplot(width = 0.65, alpha = 0.85, outlier.size = 0.8) +
      scale_fill_manual(
        values = method_palette,
        breaks = method_levels,
        labels = method_labels,
        guide = "none"
      ) +
      facet_wrap(~metric, scales = "free_y") +
      labs(
        title = "Selection stability",
        x = NULL,
        y = NULL
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11, legend_position = "none", y_grid = TRUE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(lineheight = 0.95)
      )
    save_method_sensitivity_plot(
      setup_helpers,
      "FIG_method_sensitivity_selection_stability",
      selection_plot,
      width = 8,
      height = 5
    )
    
    pareto_data <- summary_tbl |>
      dplyr::mutate(
        method_label = format_method_label(method),
        decomposition = factor(
          decomposition,
          levels = names(decomp_names),
          labels = unname(decomp_names)
        )
      )
    pareto_plot <- pareto_data |>
      ggplot(aes(
        n_items_selected,
        pc12_explained_variance_ratio,
        fill = cv_r2_leave_one_scale_out_mean,
        shape = decomposition
      )) +
      geom_point(size = 3.6, colour = "black", stroke = 0.3) +
      add_pareto_labels(data = pareto_data) +
      setup_helpers$scale_prob_fill(
        limits = range(summary_tbl$cv_r2_leave_one_scale_out_mean, na.rm = TRUE),
        name = "LOSO readout R2"
      ) +
      scale_shape_manual(values = decomp_shape_values, name = "Decomposition") +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(
        title = "Item count and leave-one-scale-out readout",
        x = "Selected items",
        y = "PC1 + PC2 explained variance",
        fill = "LOSO readout R2"
      ) +
      method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11)
    save_method_sensitivity_plot(setup_helpers, "FIG_method_sensitivity_pareto", pareto_plot, width = 8, height = 5)
    
    if (nrow(baseline_tbl)) {
      selected_fit_data <- summary_tbl |>
        dplyr::filter(method %in% unique(baseline_tbl$target_method)) |>
        dplyr::transmute(
          target_method = factor(method, levels = method_levels, labels = unname(method_labels)),
          subset_mode = "id_guided_selection",
          pc12_explained_variance_ratio = pc12_explained_variance_ratio,
          cv_r2_scale_prediction_mean = cv_r2_scale_prediction_mean
        )
      
      null_palette <- c(
        "Unstratified" = method_palette[[method_levels[[1L]]]],
        "Stratified by scale" = method_palette[[method_levels[[length(method_levels)]]]]
      )
      
      baseline_plot <- baseline_tbl |>
        dplyr::mutate(
          subset_mode = factor(
            subset_mode,
            levels = c("unstratified", "stratified_by_scale"),
            labels = names(null_palette)
          ),
          target_method = factor(target_method, levels = method_levels, labels = unname(method_labels))
        ) |>
        ggplot(aes(cv_r2_scale_prediction_mean, pc12_explained_variance_ratio, colour = subset_mode)) +
        geom_point(alpha = 0.22, size = 1.3) +
        scale_colour_manual(values = null_palette, name = "Subset mode") +
        scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
        geom_point(
          data = selected_fit_data,
          aes(cv_r2_scale_prediction_mean, pc12_explained_variance_ratio),
          inherit.aes = FALSE,
          colour = "black",
          fill = "white",
          shape = 23,
          stroke = 0.5,
          size = 3.1
        ) +
        facet_wrap(~target_method, scales = "free") +
        labs(
          title = "Matched random-subset baseline",
          x = "Mean fixed-map readout R2",
          y = "PC1 + PC2 explained variance"
        ) +
        method_sensitivity_plot_theme(setup_helpers$theme_pub, base_size = 11)
      save_method_sensitivity_plot(
        setup_helpers,
        "FIG_method_sensitivity_matched_random_subset_baseline",
        baseline_plot,
        width = 9,
        height = 5
      )
    }
  }
  
  cfg_lines <- capture.output(str(cfg, max.level = 1))
  writeLines(cfg_lines, con = file.path(cfg$out_dir, "run_config.txt"))
  write_csv(scale_key, file.path(cfg$out_dir, "expanded_scale_key.csv"))
  
  session_path <- file.path(cfg$out_dir, "sessionInfo.txt")
  zz <- file(session_path, open = "wt")
  sink(zz)
  print(sessionInfo())
  sink()
  close(zz)
  
  messagef(
    "Finished. Outputs written to %s (workers=%d)",
    cfg$out_dir,
    future_worker_count()
  )
  invisible(cfg$out_dir)
}
