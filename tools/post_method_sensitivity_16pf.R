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
  optimisation_subsample_n = 1200L,
  optimisation_w_min = 0.05,
  optimisation_multi_runs = 10L,
  optimisation_multi_min_prop = 0.35,
  optimisation_step_grid = c(0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.05),
  optimisation_batch_k = 3L,
  optimisation_batch_factor = 0.75,
  optimisation_max_iter = 1000L,
  optimisation_eval_per_iter = 50L,
  stability_reps = 20L,
  stability_half_n = 400L,
  spectrum_rank = 10L,
  baseline_reps = 250L,
  baseline_subset_modes = c("unstratified", "stratified_by_scale"),
  ncores_par = "all_but_one",
  future_scheduling = 1,
  save_plots = TRUE
)

`%||%` <- function(x, y) if (is.null(x)) y else x

messagef <- function(...) {
  message(sprintf(...))
}

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

map_with_progress_worker <- function(x, worker_fun, dots = list(), p = function() NULL) {
  out <- do.call(worker_fun, c(list(x), dots))
  p()
  out
}

map_with_progress <- function(X,
                              FUN,
                              cfg,
                              ...,
                              show_progress = TRUE,
                              parallel = TRUE,
                              future.seed = TRUE) {
  dots <- list(...)

  if (!parallel) {
    if (!show_progress) {
      return(do.call(base::lapply, c(list(X, FUN), .drop_future_args(dots))))
    }
    pb <- utils::txtProgressBar(min = 0, max = length(X), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
      out[[i]] <- do.call(FUN, c(list(X[[i]]), .drop_future_args(dots)))
      utils::setTxtProgressBar(pb, i)
    }
    return(out)
  }

  nworkers <- future_worker_count()
  if (nworkers <= 1L) {
    if (!show_progress) {
      return(do.call(base::lapply, c(list(X, FUN), .drop_future_args(dots))))
    }
    pb <- utils::txtProgressBar(min = 0, max = length(X), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
      out[[i]] <- do.call(FUN, c(list(X[[i]]), .drop_future_args(dots)))
      utils::setTxtProgressBar(pb, i)
    }
    return(out)
  }

  worker_dots <- .drop_future_args(dots)

  if (!show_progress) {
    return(FUTURE_LAPPLY(
      X,
      map_with_progress_worker,
      worker_fun = FUN,
      dots = worker_dots,
      future.seed = future.seed,
      future.scheduling = cfg$future_scheduling %||% 1L
    ))
  }

  progressr::with_progress({
    p <- progressr::progressor(steps = length(X))
    FUTURE_LAPPLY(
      X,
      map_with_progress_worker,
      worker_fun = FUN,
      dots = worker_dots,
      p = p,
      future.packages = "progressr",
      future.seed = future.seed,
      future.scheduling = cfg$future_scheduling %||% 1L
    )
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
    levs <- sort(unique(stats::na.omit(x)))
    X1[[nm]] <- ordered(x, levels = levs)
  }
  list(X = X1, type = list(ordratio = names(X1)))
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
  p <- ncol(Xdf)
  N_list <- vector("list", p)
  S_list <- vector("list", p)
  for (j in seq_len(p)) {
    wj <- rep(0, p)
    wj[[j]] <- 1
    Dj <- gower_dist(Xdf, type_list, weights = wj)
    dvec <- as.numeric(Dj)
    finite_mask <- as.numeric(is.finite(dvec))
    dvec[!is.finite(dvec)] <- 0
    N_list[[j]] <- dvec
    S_list[[j]] <- finite_mask
  }
  list(N = N_list, S = S_list, n = nrow(Xdf))
}

calc_id_from_cache <- function(num, den, n_rows) {
  Dvec <- num / pmax(den, .Machine$double.eps)
  attr(Dvec, "Size") <- n_rows
  attr(Dvec, "Diag") <- FALSE
  attr(Dvec, "Upper") <- FALSE
  class(Dvec) <- "dist"
  twonn_id_from_dist(Dvec)
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
  set.seed(seed)
  row_pool <- seq_len(nrow(X))
  n_sub <- min(as.integer(n_rows_sub), length(row_pool))
  idx_sub <- if (length(row_pool) > n_sub) sample(row_pool, n_sub) else row_pool
  prep <- prep_ord_gower(X[idx_sub, , drop = FALSE])
  cache <- make_ns_cache(prep$X, prep$type)
  vars <- names(prep$X)
  w <- setNames(rep(1, length(vars)), vars)
  num_cur <- Reduce(`+`, Map(`*`, cache$N, as.list(w)))
  den_cur <- Reduce(`+`, Map(`*`, cache$S, as.list(w)))
  id_cur <- calc_id_from_cache(num_cur, den_cur, cache$n)
  hist <- tibble::tibble(iter = 0L, twonn_id = id_cur, changed = "START", weight = NA_real_)
  hot_vars <- integer(0)
  step_grid <- sort(unique(step_grid), decreasing = FALSE)
  for (it in seq_len(as.integer(max_iter))) {
    candidates_all <- which(w > (w_min + 1e-12))
    if (!length(candidates_all)) break
    n_random <- max(10L, eval_per_iter - length(hot_vars))
    draw <- sample(candidates_all, min(length(candidates_all), n_random))
    candidates <- unique(c(hot_vars, draw))
    best <- list(id = Inf)
    for (j in candidates) {
      w_base <- w[[j]]
      for (factor in step_grid) {
        w_new <- max(w_min, w_base * factor)
        if (w_new >= (w_base - 1e-12)) next
        delta <- w_new - w_base
        id_try <- calc_id_from_cache(
          num_cur + delta * cache$N[[j]],
          den_cur + delta * cache$S[[j]],
          cache$n
        )
        if (is.finite(id_try) && id_try < best$id) {
          best <- list(id = id_try, j = j, w_new = w_new)
        }
      }
    }
    if (!is.finite(best$id) || best$id >= (id_cur - 1e-6)) {
      if (length(hot_vars)) {
        hot_vars <- integer(0)
        next
      }
      break
    }
    delta <- best$w_new - w[[best$j]]
    num_cur <- num_cur + delta * cache$N[[best$j]]
    den_cur <- den_cur + delta * cache$S[[best$j]]
    w[[best$j]] <- best$w_new
    id_cur <- best$id
    hot_vars <- unique(c(best$j, hot_vars))
    if (length(hot_vars) > 10L) hot_vars <- head(hot_vars, 10L)
    hist <- dplyr::bind_rows(
      hist,
      tibble::tibble(iter = it, twonn_id = id_cur, changed = vars[[best$j]], weight = best$w_new)
    )
    if (batch_k > 1L && (it %% 5L == 0L)) {
      remain <- setdiff(candidates_all, best$j)
      if (length(remain)) {
        remain <- sample(remain, min(length(remain), eval_per_iter))
        scores <- vapply(remain, function(j) {
          delta_b <- max(w_min, w[[j]] * batch_factor) - w[[j]]
          calc_id_from_cache(
            num_cur + delta_b * cache$N[[j]],
            den_cur + delta_b * cache$S[[j]],
            cache$n
          )
        }, numeric(1))
        take <- head(remain[order(scores)], min(batch_k, length(scores)))
        if (length(take)) {
          num_b <- num_cur
          den_b <- den_cur
          w_b <- w
          for (j in take) {
            w_new <- max(w_min, w[[j]] * batch_factor)
            delta_b <- w_new - w[[j]]
            num_b <- num_b + delta_b * cache$N[[j]]
            den_b <- den_b + delta_b * cache$S[[j]]
            w_b[[j]] <- w_new
          }
          id_b <- calc_id_from_cache(num_b, den_b, cache$n)
          if (is.finite(id_b) && id_b < (id_cur - 1e-6)) {
            num_cur <- num_b
            den_cur <- den_b
            w <- w_b
            id_cur <- id_b
            hist <- dplyr::bind_rows(
              hist,
              tibble::tibble(iter = it, twonn_id = id_cur, changed = "BATCH", weight = NA_real_)
            )
          }
        }
      }
    }
  }
  list(weights = w, history = hist, twonn_id = id_cur, idx_sub = idx_sub)
}

knee_select_items <- function(weights, w_min = 0.05) {
  w <- sort(pmax(as.numeric(weights), w_min), decreasing = TRUE)
  names(w) <- names(sort(weights, decreasing = TRUE))
  n <- length(w)
  x <- seq_len(n)
  if (n < 3L) {
    survivors <- names(w)
    return(list(survivors = survivors, threshold = min(w), weights_sorted = w))
  }
  y_norm <- (w - min(w)) / (max(w) - min(w) + 1e-12)
  x_norm <- (x - min(x)) / (max(x) - min(x) + 1e-12)
  d <- abs(y_norm - (1 - x_norm))
  plateau_end <- max(which(w >= 0.95))
  if (!is.finite(plateau_end) || plateau_end == n) plateau_end <- 1L
  knee_local <- which.max(d[plateau_end:n]) + plateau_end - 1L
  keep_n <- min(n, max(knee_local, ceiling(0.10 * n), 5L))
  survivors <- names(w)[seq_len(keep_n)]
  list(
    survivors = survivors,
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
    sel_mat[sel_info[[r]]$survivors, r] <- 1L
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

fit_fixed_selection_from_raw <- function(X_raw,
                                         selected_items,
                                         item_weights,
                                         decomp_method,
                                         spectrum_rank) {
  prep <- prepare_standardised_matrix(X_raw)
  fit_fixed_selection_from_prepared(
    X_std_all = prep$std,
    selected_items = selected_items,
    item_weights = item_weights,
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

fit_one_method <- function(X_raw,
                           X_std_all,
                           weighting_mode,
                           decomp_method,
                           cfg) {
  limit_worker_threads()
  if (identical(weighting_mode, "id_guided")) {
    run_ids <- seq_len(as.integer(cfg$optimisation_multi_runs))
    opt_list <- lapply(run_ids, function(r) {
      seed_r <- if (r == 1L) cfg$seed else cfg$seed + r
      optimise_item_weights(
        X = X_raw,
        n_rows_sub = cfg$optimisation_subsample_n,
        w_min = cfg$optimisation_w_min,
        step_grid = cfg$optimisation_step_grid,
        batch_k = cfg$optimisation_batch_k,
        batch_factor = cfg$optimisation_batch_factor,
        max_iter = cfg$optimisation_max_iter,
        eval_per_iter = cfg$optimisation_eval_per_iter,
        seed = seed_r
      )
    })
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

split_half_full_pipeline_refit_worker <- function(b,
                                                  X_raw,
                                                  scale_scores,
                                                  weighting_mode,
                                                  decomp_method,
                                                  worker_cfg,
                                                  n,
                                                  half_n) {
  limit_worker_threads()
  set.seed(worker_cfg$seed + b)
  idx <- sample(seq_len(n), size = 2L * half_n)
  idx_a <- idx[seq_len(half_n)]
  idx_b <- idx[half_n + seq_len(half_n)]
  Xa_raw <- X_raw[idx_a, , drop = FALSE]
  Xb_raw <- X_raw[idx_b, , drop = FALSE]
  prep_a <- prepare_standardised_matrix(Xa_raw)
  prep_b <- prepare_standardised_matrix(Xb_raw)
  common_items <- intersect(colnames(prep_a$std), colnames(prep_b$std))
  Xa_raw <- Xa_raw[, common_items, drop = FALSE]
  Xb_raw <- Xb_raw[, common_items, drop = FALSE]
  Xa_std <- prep_a$std[, common_items, drop = FALSE]
  Xb_std <- prep_b$std[, common_items, drop = FALSE]

  cfg_refit <- worker_cfg
  cfg_refit$optimisation_multi_runs <- 1L
  fit_a <- tryCatch(
    fit_one_method(Xa_raw, Xa_std, weighting_mode, decomp_method, cfg_refit),
    error = function(e) NULL
  )
  fit_b <- tryCatch(
    fit_one_method(Xb_raw, Xb_std, weighting_mode, decomp_method, cfg_refit),
    error = function(e) NULL
  )

  if (is.null(fit_a) || is.null(fit_b)) {
    return(tibble::tibble(
      rep = b,
      axis_corr_mean = NA_real_,
      item_rmse = NA_real_,
      scale_regression_vector_cosine = NA_real_,
      n_anchor_items = NA_integer_,
      selection_jaccard = NA_real_,
      aligned_weight_correlation = NA_real_,
      n_selected_a = NA_integer_,
      n_selected_b = NA_integer_,
      n_selected_intersection = NA_integer_
    ))
  }

  geom_metrics <- geometry_stability_metrics(
    fit_a = fit_a,
    fit_b = fit_b,
    scale_scores_a = scale_scores[idx_a, , drop = FALSE],
    scale_scores_b = scale_scores[idx_b, , drop = FALSE],
    anchor_items = common_items
  )
  sel_metrics <- selection_stability_metrics(fit_a, fit_b)

  dplyr::bind_cols(tibble::tibble(rep = b), geom_metrics, sel_metrics)
}

split_half_fixed_selection_worker <- function(b,
                                              X_raw,
                                              scale_scores,
                                              fixed_items,
                                              fixed_weights,
                                              decomp_method,
                                              worker_cfg,
                                              n,
                                              half_n) {
  limit_worker_threads()
  set.seed(worker_cfg$seed + b)
  idx <- sample(seq_len(n), size = 2L * half_n)
  idx_a <- idx[seq_len(half_n)]
  idx_b <- idx[half_n + seq_len(half_n)]
  prep_a <- prepare_standardised_matrix(X_raw[idx_a, , drop = FALSE])
  prep_b <- prepare_standardised_matrix(X_raw[idx_b, , drop = FALSE])
  common_selected <- intersect(fixed_items, intersect(colnames(prep_a$std), colnames(prep_b$std)))

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
    return(tibble::tibble(
      rep = b,
      axis_corr_mean = NA_real_,
      item_rmse = NA_real_,
      scale_regression_vector_cosine = NA_real_,
      n_anchor_items = NA_integer_,
      n_fixed_items_used = length(common_selected)
    ))
  }

  geom_metrics <- geometry_stability_metrics(
    fit_a = fit_a,
    fit_b = fit_b,
    scale_scores_a = scale_scores[idx_a, , drop = FALSE],
    scale_scores_b = scale_scores[idx_b, , drop = FALSE],
    anchor_items = common_selected
  )

  dplyr::bind_cols(
    tibble::tibble(rep = b, n_fixed_items_used = length(common_selected)),
    geom_metrics
  )
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
  half_n <- min(as.integer(cfg$stability_half_n), floor(n / 2))
  cfg_worker <- cfg

  rows <- map_with_progress(
    X = as.list(seq_len(as.integer(cfg$stability_reps))),
    FUN = split_half_full_pipeline_refit_worker,
    cfg = cfg,
    X_raw = X_raw,
    scale_scores = scale_scores,
    weighting_mode = weighting_mode,
    decomp_method = decomp_method,
    worker_cfg = cfg_worker,
    n = n,
    half_n = half_n,
    show_progress = show_progress,
    parallel = parallel
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
  half_n <- min(as.integer(cfg$stability_half_n), floor(n / 2))
  fixed_weights <- fit_reference$weights
  fixed_items <- fit_reference$selected_items
  cfg_worker <- cfg

  rows <- map_with_progress(
    X = as.list(seq_len(as.integer(cfg$stability_reps))),
    FUN = split_half_fixed_selection_worker,
    cfg = cfg,
    X_raw = X_raw,
    scale_scores = scale_scores,
    fixed_items = fixed_items,
    fixed_weights = fixed_weights,
    decomp_method = fit_reference$decomp_method,
    worker_cfg = cfg_worker,
    n = n,
    half_n = half_n,
    show_progress = show_progress,
    parallel = parallel
  )
  dplyr::bind_rows(rows)
}

split_half_stability <- function(X_raw,
                                 scale_scores,
                                 weighting_mode,
                                 decomp_method,
                                 cfg) {
  split_half_full_pipeline_refit_stability(
    X_raw = X_raw,
    scale_scores = scale_scores,
    weighting_mode = weighting_mode,
    decomp_method = decomp_method,
    cfg = cfg
  )
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

format_method_label <- function(method) {
  dplyr::recode(
    method,
    pca_uniform = "PCA\nuniform",
    robpca_uniform = "Robust PCA\nuniform",
    pca_id_guided = "PCA\nID-guided",
    robpca_id_guided = "Robust PCA\nID-guided",
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

normalise_bundle_method <- function(decomp_method) {
  decomp_method <- tolower(decomp_method)
  if (identical(decomp_method, "robust_pca")) return("robpca")
  decomp_method
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
  list(
    participant_id = as.character(bundle$participant_id),
    selected_items = as.character(bundle$selected_items),
    weights = bundle$weights,
    pc_scores_2d = pc_scores_2d,
    spectrum = spectrum,
    explained_variance_ratio = explained_variance_ratio,
    spectral_effective_rank = spectral_effective_rank(spectrum),
    item_component_correlations = item_component_correlations,
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
    tibble::tibble(
      method = fit$method,
      item_id = names(fit$weights),
      weight = as.numeric(fit$weights),
      selected = names(fit$weights) %in% fit$selected_items
    )
  }))
  weights_tbl <- dplyr::left_join(weights_tbl, qbank |> dplyr::select(item_id, item_text), by = "item_id")
  write_csv(weights_tbl, file.path(cfg$out_dir, "method_sensitivity_item_weights.csv"))

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
    scale_prediction_for_method(fit$method, fit$pc_scores_2d, scale_scores, cfg)
  }))
  write_csv(scale_prediction_tbl, file.path(cfg$out_dir, "method_sensitivity_scale_prediction.csv"))

  leave_one_scale_out_prediction_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    messagef("Leave-one-scale-out prediction for %s", fit$method)
    leave_one_scale_out_prediction_for_method(
      fit_reference = fit,
      X_raw = X_raw,
      scale_scores = scale_scores,
      scale_key = scale_key,
      cfg = cfg,
      prep = prep_full
    )
  }))
  write_csv(
    leave_one_scale_out_prediction_tbl,
    file.path(cfg$out_dir, "method_sensitivity_leave_one_scale_out_prediction.csv")
  )

  scale_regression_vectors_tbl <- dplyr::bind_rows(lapply(fits, function(fit) {
    scale_regression_vectors(scale_scores, fit$pc_scores_2d) |>
      dplyr::mutate(method = fit$method, .before = 1)
  }))
  write_csv(
    scale_regression_vectors_tbl,
    file.path(cfg$out_dir, "method_sensitivity_scale_regression_vectors.csv")
  )

  full_pipeline_stability_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(method_specs)), function(i) {
    weighting_mode <- method_specs$weighting_mode[[i]]
    decomp_method <- method_specs$decomp_method[[i]]
    label <- method_label(weighting_mode, decomp_method)
    messagef("Split-half full-refit stability for %s", label)
    split_half_full_pipeline_refit_stability(
      X_raw = X_raw,
      scale_scores = scale_scores,
      weighting_mode = weighting_mode,
      decomp_method = decomp_method,
      cfg = cfg
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
      cfg = cfg
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
      prep = prep_full
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
        y = "Explained variance ratio"
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
        title = "Cross-validated scale score prediction",
        x = NULL,
        y = "Cross-validated R-squared"
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
        title = "Leave-one-scale-out prediction",
        x = NULL,
        y = "Cross-validated R-squared"
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
        name = "Leave-one-scale-out CV R-squared"
      ) +
      scale_shape_manual(values = decomp_shape_values, name = "Decomposition") +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(
        title = "Compression vs leave-one-scale-out prediction",
        x = "Selected items",
        y = "PC1 + PC2 explained variance ratio",
        fill = "Leave-one-scale-out CV R-squared"
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
          x = "Mean cross-validated R-squared",
          y = "PC1 + PC2 explained variance ratio"
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
