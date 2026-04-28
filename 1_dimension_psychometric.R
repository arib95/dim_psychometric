# ==============================================================================
#                               dimension_psychometric.R
# ==============================================================================

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

# --- Data Preparation for Gower Distance ---
prep_X_for_gower <- function(X, rare_prop = 0.01, do_jitter = TRUE, seed = NULL) {
  X1 <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
  
  # Ensure character cols are factors
  for (nm in names(X1)) {
    v <- X1[[nm]]
    if (is.character(v)) X1[[nm]] <- factor(v)
  }
  
  # Collapses rare factor levels into NA or a residual category
  drop_rare <- function(f, prop) {
    if (!is.factor(f) || is.ordered(f)) {
      return(f)
    }
    tb <- prop.table(table(f))
    keep <- names(tb)[tb >= prop]
    f <- factor(ifelse(f %in% keep, as.character(f), NA), exclude = NULL)
    droplevels(f)
  }
  
  X1 <- as.data.frame(lapply(X1, drop_rare, prop = rare_prop), stringsAsFactors = FALSE)

  if (isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
    X1 <- as.data.frame(lapply(X1, function(v) {
      if (is.ordered(v)) factor(v, levels = levels(v), ordered = FALSE, exclude = NULL) else v
    }), stringsAsFactors = FALSE)
  }
  
  # Add microscopic jitter to numerics to prevent ties in distance
  if (isTRUE(do_jitter)) {
    for (nm in names(X1)) {
      if (is.numeric(X1[[nm]])) {
        sdv <- stats::sd(X1[[nm]], na.rm = TRUE)
        if (is.finite(sdv) && sdv > 0) {
          eps <- 1e-6 * sdv
          if (!is.null(seed)) {
            s_nm <- .seed_from_key(seed, paste0("prep_jitter::", nm))
            X1[[nm]] <- X1[[nm]] + .with_seed(s_nm, stats::rnorm(length(X1[[nm]]), 0, eps))
          } else {
            X1[[nm]] <- X1[[nm]] + stats::rnorm(length(X1[[nm]]), 0, eps)
          }
        }
      }
    }
  }
  
  # Identify column types for daisy()
  ord_cols <- names(X1)[vapply(X1, is.ordered, logical(1))]
  fac_cols <- names(X1)[vapply(X1, function(z) is.factor(z) && !is.ordered(z), logical(1))]
  
  # Internal helper to detect binary
  .is_binary <- function(x) length(unique(na.omit(x))) == 2
  bin_cols <- if (isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
    character(0)
  } else {
    fac_cols[vapply(X1[fac_cols], .is_binary, logical(1))]
  }
  
  type_list <- list()
  if (length(bin_cols)) type_list$asymm <- bin_cols
  if (length(ord_cols)) type_list$ordratio <- ord_cols
  
  w <- setNames(rep(1, ncol(X1)), names(X1))
  
  list(X = X1, type = type_list, weights = w)
}

# --- Gower Distance Wrapper ---
gower_dist <- function(Xdf, type_list = NULL, weights = NULL) {
  if (!is.null(type_list)) {
    type_list <- lapply(type_list, function(cols) intersect(cols, names(Xdf)))
    type_list <- type_list[lengths(type_list) > 0]
    if (!length(type_list)) type_list <- NULL
  }
  
  if (is.null(weights)) {
    weights <- rep(1, ncol(Xdf))
  } else if (length(weights) == 1) {
    weights <- rep(weights, ncol(Xdf))
  } else if (!is.null(names(weights))) {
    weights <- weights[names(Xdf)]
  }
  
  stopifnot(length(weights) == ncol(Xdf))
  cluster::daisy(Xdf, metric = "gower", type = type_list, weights = weights)
}

dist_health_console <- function(D, name = "D", eps0 = 1e-12) {
  M <- as.matrix(D)
  diag(M) <- Inf
  n <- nrow(M)

  negM <- -M
  i1 <- max.col(negM, ties.method = "first")
  d1 <- M[cbind(seq_len(n), i1)]
  negM[cbind(seq_len(n), i1)] <- -Inf
  i2 <- max.col(negM, ties.method = "first")
  d2 <- M[cbind(seq_len(n), i2)]

  v <- as.numeric(D)
  v <- v[is.finite(v)]
  uq <- length(unique(round(v, 10)))
  prop0 <- mean(v <= eps0)

  qd1 <- stats::quantile(
    d1[is.finite(d1)],
    probs = c(0, 0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 1),
    na.rm = TRUE
  )
  eps_adapt <- max(1e-12, as.numeric(qd1[[2]]) * 0.1)

  d1c <- pmax(d1, eps_adapt)
  d2c <- pmax(d2, d1c + eps_adapt)
  r <- d2c / d1c
  lr <- log(r[is.finite(r) & r > 1])

  cat("\n================ DIST HEALTH:", name, "================\n")
  cat(sprintf("n=%d | unique_dist(rounded 1e-10)=%d | prop(dist<=1e-12)=%.4f\n", n, uq, prop0))
  cat("d1 quantiles:\n")
  print(qd1)
  cat(sprintf(
    "eps_adapt=%.3g | mean(log r)=%.3f | sd(log r)=%.3f | TwoNN=%.3f\n",
    eps_adapt,
    mean(lr, na.rm = TRUE),
    stats::sd(lr, na.rm = TRUE),
    1 / mean(lr, na.rm = TRUE)
  ))
  cat(sprintf("tie rate (d1 <= eps_adapt): %.4f\n", mean(d1 <= eps_adapt, na.rm = TRUE)))
}

knn_from_dist <- function(D, k = 15) {
  M <- as.matrix(D)
  diag(M) <- Inf
  t(apply(M, 1, function(r) order(r)[1:k]))
}

mean_jaccard_knn <- function(idxA, idxB) {
  stopifnot(nrow(idxA) == nrow(idxB), ncol(idxA) == ncol(idxB))
  n <- nrow(idxA)
  s <- 0
  for (i in seq_len(n)) {
    a <- idxA[i, ]
    b <- idxB[i, ]
    s <- s + length(intersect(a, b)) / length(union(a, b))
  }
  s / n
}

constant_profile_value <- function(Xdf) {
  n <- nrow(Xdf)
  if (!n || !ncol(Xdf)) return(rep(NA_character_, n))

  ref <- as.character(Xdf[[1L]])
  same <- rep(TRUE, n)
  if (ncol(Xdf) > 1L) {
    for (j in 2:ncol(Xdf)) {
      same <- same & (as.character(Xdf[[j]]) == ref)
    }
  }

  out <- rep(NA_character_, n)
  out[same] <- ref[same]
  out
}

warn_diag_subsample <- function(tag, n_used, n_total, cap, pool_label = "rows") {
  warning(sprintf(
    "[%s] Using %d/%d %s because full distance diagnostics materialise dist objects and dense matrices (O(n^2) memory/time); cap=%d.",
    tag, n_used, n_total, pool_label, cap
  ))
}

# --- Nearest Neighbor Helpers (Vectorized/Fast) ---

# Fast single pass: first and second NN using max.col
.two_nn_from_distvec <- function(D, eps = .Machine$double.eps) {
  # Convert compact dist vector to full matrix (Internal C function, very fast)
  # Memory note: For N=1000, this is only ~8MB.
  M <- as.matrix(D)
  n <- nrow(M)
  
  if (n < 2L) return(list(d1 = numeric(0), d2 = numeric(0)))
  
  # Set diagonal to Inf so we don't find the point itself
  diag(M) <- Inf
  
  # Find 1st NN (min distance is max negative distance)
  negM <- -M
  idx1 <- max.col(negM, ties.method = "first")
  d1 <- M[cbind(seq_len(n), idx1)]
  
  # Find 2nd NN by masking the first
  negM[cbind(seq_len(n), idx1)] <- -Inf
  idx2 <- max.col(negM, ties.method = "first")
  d2 <- M[cbind(seq_len(n), idx2)]
  
  # Safety clamp
  d1 <- pmax(d1, eps)
  d2 <- pmax(d2, d1 + eps)
  
  list(d1 = d1, d2 = d2)
}

# TwoNN intrinsic dimension without materialising the full matrix
twonn_id_from_dist <- function(D, eps = 1e-8, trim = 0.02) {
  n <- attr(D, "Size")
  if (n < 3L) return(NA_real_)
  
  nn <- .two_nn_from_distvec(D, eps)
  r  <- nn$d2 / nn$d1
  r  <- r[is.finite(r) & r > 1]
  
  if (!length(r)) return(NA_real_)
  
  lr <- sort(log(r))
  k  <- floor(trim * length(lr))
  if (k > 0 && 2 * k < length(lr)) lr <- lr[(k + 1):(length(lr) - k)]
  
  1 / mean(lr)
}

first_nn_d1 <- function(D) {
  n <- attr(D, "Size")
  if (n < 2L) return(rep(Inf, n))
  .two_nn_from_distvec(D)$d1
}

# --- Dedup & Coreset Helpers ---

collapse_curve <- function(D, eps_grid) {
  n0 <- attr(D, "Size")
  
  get_ngroups <- function(eps) {
    length(unique(hclust(D, method = "complete") |> cutree(h = eps)))
  }
  
  data.frame(
    eps = eps_grid,
    n_groups = sapply(eps_grid, get_ngroups),
    prop_retained = NA_real_
  ) |>
    transform(prop_retained = n_groups / n0)
}

complete_groups <- function(D, eps) {
  hclust(D, method = "complete") |> cutree(h = eps)
}

# Core-set finder without matrix expansion
twonn_core_by_slope <- function(D, min_frac = 0.30, w = 20, slope_tol = 0.08, rmse_tol = 0.10) {
  n <- attr(D, "Size")
  if (n < 8L) return(seq_len(n))
  
  nn <- .two_nn_from_distvec(D, eps = .Machine$double.eps)
  mu <- pmax(nn$d2 / pmax(nn$d1, .Machine$double.eps), 1 + 1e-12)
  
  ord <- order(mu)
  x <- log(mu[ord])
  m <- length(x)
  
  u <- (seq_len(m) - 0.5) / (m + 1)
  y <- log(1 - u)
  
  k0 <- max(20L, floor(min_frac * m))
  slope <- rmse <- rep(NA_real_, m)
  
  for (k in k0:(m - 2L)) {
    fit <- stats::lm(y[1:k] ~ x[1:k])
    slope[k] <- coef(fit)[2]
    rmse[k] <- sqrt(mean(residuals(fit)^2))
  }
  
  ok <- which(is.finite(slope) & is.finite(rmse))
  if (!length(ok)) return(ord[seq_len(max(3L, k0))])
  
  pick <- function(k) {
    L <- max(k0, k - w + 1L)
    s <- slope[L:k]
    (max(s, na.rm = TRUE) - min(s, na.rm = TRUE) <= slope_tol) && (rmse[k] <= rmse_tol)
  }
  
  ks <- ok[vapply(ok, pick, logical(1))]
  k_star <- if (length(ks)) max(ks) else floor(0.6 * m)
  
  ord[seq_len(max(3L, min(k_star, m - 2L)))]
}

core_band_idx <- function(D, k = 10, band = c(0.20, 0.70)) {
  M <- as.matrix(D)
  diag(M) <- Inf
  
  kth <- function(r, k) {
    rf <- r[is.finite(r)]
    if (!length(rf)) return(NA_real_)
    k_eff <- min(k, length(rf))
    sort(rf, partial = k_eff)[k_eff]
  }
  
  rk <- apply(M, 1, kth, k = k)
  ok <- is.finite(rk)
  if (!any(ok)) return(integer(0))
  
  q <- stats::quantile(rk[ok], band, na.rm = TRUE)
  which(ok & rk >= q[1] & rk <= q[2])
}

# TwoNN + LB estimators
lb_mle_id <- function(Dm, k_lo = 5, k_hi = 15) {
  Dm <- as.matrix(Dm)
  n <- nrow(Dm)
  diag(Dm) <- Inf
  if (n <= k_lo) return(NA_real_)
  
  k_hi <- max(k_lo, min(k_hi, n - 1))
  
  ids <- sapply(k_lo:k_hi, function(k) {
    nn <- t(apply(Dm, 1L, function(r) {
      rf <- r[is.finite(r)]
      m <- length(rf)
      if (m < k) return(rep(NA_real_, k))
      sort(rf, partial = k)[1:k]
    }))
    
    if (!nrow(nn)) return(NA_real_)
    
    l <- log(nn[, k, drop = TRUE] / nn[, 1:(k - 1), drop = FALSE])
    d <- 1 / rowMeans(l, na.rm = TRUE)
    mean(d[is.finite(d)], na.rm = TRUE)
  })
  
  mean(ids, na.rm = TRUE)
}

# --- Weight Optimization Helpers ---

make_NS_cache <- function(Xdf, type = NULL) {
  p <- ncol(Xdf)
  N_list <- vector("list", p)
  S_list <- vector("list", p)
  
  for (j in seq_len(p)) {
    w1 <- rep(0, p)
    w1[j] <- 1
    Dj <- cluster::daisy(Xdf, metric = "gower", type = type, weights = w1)
    vD <- as.numeric(Dj)
    ok <- as.numeric(is.finite(vD))
    vD[!is.finite(vD)] <- 0
    N_list[[j]] <- vD
    S_list[[j]] <- ok
  }
  list(N = N_list, S = S_list, n = nrow(Xdf))
}

optimise_gower_weights_constrained <- function(X, init_weights, allow_update,
                                               objective = "TwoNN_all",
                                               w_min = W_MIN,
                                               step_grid = W_STEP_GRID, 
                                               batch_k = W_BATCH_K,
                                               batch_factor = W_BATCH_FACTOR,
                                               max_iter = W_MAX_ITERS,
                                               n_rows_sub = N_ROWS_SUB,
                                               ncores = NULL, 
                                               seed_jitter = SEED_JITTER,
                                               reps_idx = NULL,
                                               core_idx_rep = NULL,
                                               base_seed = SEED_GLOBAL,
                                               verbose = TRUE,
                                               plot_progress = TRUE,
                                               progress_fun = NULL) {
  seed_sub  <- .seed_from_key(base_seed, "optim::subsample_rows")
  seed_iter <- .seed_from_key(base_seed, "optim::iter_sampling")
  seed_prep <- seed_jitter
  
  # Internal ID calc from numerator/denominator
  calc_id_fast <- function(num, den, n_rows) {
    Dvec <- num / pmax(den, .Machine$double.eps)
    attr(Dvec, "Size")  <- n_rows
    attr(Dvec, "Diag")  <- FALSE
    attr(Dvec, "Upper") <- FALSE
    class(Dvec) <- "dist"
    twonn_id_from_dist(Dvec) 
  }
  
  # Prepare subset
  px  <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE, seed = seed_prep)
  X0  <- px$X
  typ <- px$type
  
  row_pool <- if (!is.null(reps_idx)) reps_idx else seq_len(nrow(X0))
  if (is.null(n_rows_sub) || !length(n_rows_sub)) {
    eff_n_sub <- 10000L
  } else if (!is.finite(n_rows_sub) || as.integer(n_rows_sub) <= 0L) {
    eff_n_sub <- length(row_pool)
  } else {
    eff_n_sub <- min(as.integer(n_rows_sub), length(row_pool))
  }

  if (eff_n_sub > 1500L) {
    warning(sprintf(
      "[optim] Large n_sub=%d requested; Gower optimisation can be slow, especially with multi-run enabled.",
      eff_n_sub
    ))
  }
  
  ix_sub_from_reps <- if (length(row_pool) > eff_n_sub) {
    if (isTRUE(FIX_REP_SUBSET)) {
      head(row_pool, eff_n_sub)
    } else {
      .with_seed(seed_sub, sample(row_pool, eff_n_sub))
    }
  } else {
    row_pool
  }
  
  Xs <- X0[ix_sub_from_reps, , drop = FALSE]
  n_sub <- nrow(Xs)
  
  vars <- colnames(Xs)
  w <- init_weights[vars]
  allow_update <- allow_update[vars]
  w[!allow_update] <- pmax(w_min, w[!allow_update])
  
  # Build cache
  cache <- make_NS_cache(Xs, type = typ)
  
  # Initialize State
  num_cur <- Reduce(`+`, Map(`*`, cache$N, as.list(w)))
  den_cur <- Reduce(`+`, Map(`*`, cache$S, as.list(w)))
  
  id0 <- calc_id_fast(num_cur, den_cur, n_sub)
  hist <- data.frame(iter = 0L, ID = id0, changed = NA_character_, note = NA_character_)
  
  if (verbose) cat(sprintf("[optim] Start ID: %.3f | N_sub: %d | Mode: STOCHASTIC\n", id0, n_sub))
  
  max_iter_eff <- if (is.null(max_iter) || !is.finite(max_iter)) 1000L else as.integer(max_iter)
  
  # Optimization Config
  N_SAMPLE_PER_ITER <- 50 
  step_grid <- sort(step_grid, decreasing = FALSE) 
  hot_vars <- integer(0) # Track momentum
  
  id <- id0
  
  for (it in seq_len(max_iter_eff)) {
    
    can_all <- which(allow_update & (w > w_min + 1e-12))
    if (!length(can_all)) break
    
    # 1. Stochastic Sampling: Keep hot vars, fill rest with randoms
    n_rnd <- max(10, N_SAMPLE_PER_ITER - length(hot_vars))
    rnd_vars <- .with_seed(
      as.integer(seed_iter + it),
      sample(can_all, min(length(can_all), n_rnd))
    )
    can_iter <- unique(c(hot_vars, rnd_vars))
    
    # 2. Evaluate Candidates (Greedy)
    best_res <- list(id = Inf)
    
    for (j in can_iter) {
      w_base <- w[j]
      
      for (factor in step_grid) {
        w_new <- max(w_min, w_base * factor)
        if (w_new >= w_base - 1e-12) next 
        
        delta <- w_new - w_base
        val <- calc_id_fast(num_cur + delta * cache$N[[j]], 
                            den_cur + delta * cache$S[[j]], n_sub)
        
        if (is.finite(val) && val < best_res$id) {
          best_res <- list(id = val, j = j, w_new = w_new)
        }
      }
    }
    
    # 3. Commit Best Move
    changed <- FALSE
    if (is.finite(best_res$id) && best_res$id < id - 1e-6) {
      jbest <- best_res$j
      wbest <- best_res$w_new
      delta <- wbest - w[jbest]
      
      num_cur <- num_cur + delta * cache$N[[jbest]]
      den_cur <- den_cur + delta * cache$S[[jbest]]
      w[jbest] <- wbest
      id <- best_res$id
      changed <- TRUE
      
      # Update "Hot" list (Keep this var, drop oldest if too full)
      hot_vars <- unique(c(jbest, hot_vars))
      if (length(hot_vars) > 10) hot_vars <- head(hot_vars, 10)
      
      hist <- rbind(hist, data.frame(iter = it, ID = id, changed = vars[jbest], note = sprintf("%.3f", wbest)))
      if (verbose) cat(sprintf("   iter %d: %s -> %.3f (ID: %.3f)\n", it, vars[jbest], wbest, id))
      
    } else {
      # No improvement? Flush momentum or quit
      if (length(hot_vars) > 0) {
        hot_vars <- integer(0) 
        if (verbose) cat("   [optim] Momentum lost, flushing hot vars.\n")
      } else {
        if (verbose) cat("[optim] No improvement in random subset.\n")
        break
      }
    }

    if (!is.null(progress_fun) && (it == 1L || it %% 10L == 0L)) {
      progress_fun(list(iter = it, ID = id))
    }
    
    # 4. Batch Descent (Every 5th iter)
    if (batch_k > 1 && (it %% 5 == 0)) {
      remain <- setdiff(can_all, best_res$j)
      remain <- .with_seed(
        as.integer(seed_iter + 10000L + it),
        sample(remain, min(length(remain), N_SAMPLE_PER_ITER))
      )
      
      if (length(remain) > 0) {
        scores <- numeric(length(remain))
        for (i in seq_along(remain)) {
          j <- remain[i]
          delta_b <- (max(w_min, w[j] * batch_factor)) - w[j]
          scores[i] <- calc_id_fast(num_cur + delta_b * cache$N[[j]], 
                                    den_cur + delta_b * cache$S[[j]], n_sub)
        }
        
        ord <- order(scores)
        take_idx <- head(ord, min(batch_k, length(ord)))
        take_vars <- remain[take_idx]
        
        if (length(take_vars)) {
          num_b <- num_cur; den_b <- den_cur
          w_b <- w
          for (j in take_vars) {
            wn <- max(w_min, w[j] * batch_factor)
            d  <- wn - w[j]
            num_b <- num_b + d * cache$N[[j]]
            den_b <- den_b + d * cache$S[[j]]
            w_b[j] <- wn
          }
          id_b <- calc_id_fast(num_b, den_b, n_sub)
          
          if (is.finite(id_b) && id_b < id - 1e-6) {
            num_cur <- num_b; den_cur <- den_b
            w <- w_b
            id <- id_b
            changed <- TRUE
            hist <- rbind(hist, data.frame(iter = it, ID = id, changed = "BATCH", note = paste(length(take_vars), "vars")))
            if (verbose) cat(sprintf("   iter %d: [BATCH] x%.2f on %d vars (ID: %.3f)\n", it, batch_factor, length(take_vars), id))
          }
        }
      }
    }
  }
  
  if (plot_progress && requireNamespace("ggplot2", quietly = TRUE)) {
    gp <- ggplot2::ggplot(hist, ggplot2::aes(iter, ID)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(title = sprintf("Stochastic Gower Optim (N=%d, Batch=%d)", n_sub, N_SAMPLE_PER_ITER)) +
      ggplot2::theme_minimal()
    print(gp)
  }
  
  list(weights = w, history = hist, final_ID = id, idx_used = ix_sub_from_reps)
}

# --- Post-Optimization Selection ---

# knee_triangle <- function(w) {
#   y <- sort(as.numeric(w), decreasing = TRUE)
#   n <- length(y)
#   
#   if (n < 3L) {
#     return(list(k = n, thr = if (n) y[n] else NA_real_, curve = data.frame(i = seq_len(n), w = y, d = rep(0, n))))
#   }
#   
#   x <- seq_len(n)
#   num <- abs((y[n] - y[1]) * x - (n - 1) * y + n * y[1] - y[n])
#   den <- sqrt((y[n] - y[1])^2 + (n - 1)^2)
#   d <- num / den
#   k <- which.max(d)
#   
#   list(k = k, thr = y[k], curve = data.frame(i = x, w = y, d = d))
# }

knee_satopaa <- function(w) {
  y <- sort(as.numeric(w), decreasing = TRUE)
  n <- length(y)
  x <- seq_len(n)
  
  if (n < 3L) return(list(k = n, thr = if (n) y[n] else 0.0))
  
  # Normalize to Unit Square
  y_norm <- (y - min(y)) / (max(y) - min(y))
  x_norm <- (x - min(x)) / (max(x) - min(x))
  
  # Calculate Distance from Diagonal 
  # The sensitivity line is y = 1 - x (for decreasing curves).
  # Use abs() to detect 'concave' knees (Toes) that dip below the line.
  d <- abs(y_norm - (1 - x_norm))
  
  # Find the point of maximum deviation
  k <- which.max(d)
  
  list(k = k, thr = y[k])
}

survivors_from_weights <- function(w, 
                                   w_min = W_MIN,
                                   kmin = NULL,
                                   make_plot = TRUE,
                                   plot_file = "FIG_weight_curve_knee.png") {
  
  w_names <- names(w)
  if (is.null(w_names)) w_names <- paste0("V", seq_along(w))
  w <- pmax(w_min, as.numeric(w))
  names(w) <- w_names
  p <- length(w)
  
  w_sorted <- sort(w, decreasing = TRUE)
  
  # Plateau Detection
  plateau_thresh <- 0.95  # Treat anything > 0.99 as part of the ceiling
  idx_plateau_end <- max(which(w_sorted >= plateau_thresh))
  
  if (!is.finite(idx_plateau_end) || idx_plateau_end == p) {
    idx_start <- 1
  } else {
    idx_start <- idx_plateau_end
  }
  
  w_slope <- w_sorted[idx_start:p]
  
  # Kneedle on the Slope
  knee_res <- knee_satopaa(w_slope)
  
  # Map back to global indices
  # The knee index is relative to w_slope, so we add the plateau offset
  cut_idx <- idx_start + knee_res$k - 1
  thr_knee <- w_sorted[cut_idx]
  
  if (is.null(kmin)) kmin <- max(30L, ceiling(0.10 * p))
  cut_idx <- min(p, max(cut_idx, kmin))
  thr_final <- w_sorted[cut_idx]
  
  survivors <- names(w_sorted)[seq_len(cut_idx)]
  active_eps <- GOWER_ACTIVE_EPS
  if (!is.finite(active_eps) || active_eps <= 0) active_eps <- 1e-8
  active_eps <- max(active_eps, .Machine$double.eps)
  active_survivors <- survivors[w_sorted[seq_len(cut_idx)] > (w_min + active_eps)]
  
  cat(sprintf("[Selection] Plateau=%d | Knee Offset=%d | Final: %d/%d (thr=%.4f)\n", 
              idx_start, knee_res$k, length(survivors), p, thr_final))
  
  if (isTRUE(make_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    curve_df <- data.frame(rank = seq_along(w_sorted), weight = w_sorted)
    
    # Points to visualize
    pt_plateau <- curve_df[idx_start, , drop = FALSE]
    pt_final   <- curve_df[cut_idx, , drop = FALSE]
    
    pplt <- ggplot2::ggplot(curve_df, ggplot2::aes(rank, weight)) +
      ggplot2::geom_line() +
      # Show the plateau cutoff
      ggplot2::geom_point(data = pt_plateau, ggplot2::aes(color = "Plateau End"), size = 2) +
      # Show the final cut
      ggplot2::geom_point(data = pt_final, ggplot2::aes(color = "Selected Cut"), shape = 4, size = 4, stroke = 2) +
      ggplot2::scale_color_manual(values = c("Plateau End" = "orange", "Selected Cut" = "blue")) +
      ggplot2::labs(x = "Rank (descending)", y = "Weight", color = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = c(0.8, 0.8))
    
    print(pplt)
    save_plot_gg(plot_file, pplt, width = 6, height = 4, dpi = 300)
  }
  
  list(
    survivors = survivors,
    active_survivors = active_survivors,
    active_eps = active_eps,
    thr_tail = thr_final,
    w_clamped = w
  )
}

# --- PCA and Residuals Helpers ---

design_with_map <- function(X) {
  Xg <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
  if (!ncol(Xg)) stop("[design_with_map] input has 0 columns")
  
  # Filter constant cols
  keep <- vapply(Xg, function(v) length(unique(na.omit(v))) >= 1L, logical(1))
  if (!any(keep)) stop("[design_with_map] all columns are NA-only")
  Xg <- Xg[, keep, drop = FALSE]
  
  # Type standardization for modeling
  for (nm in names(Xg)) {
    v <- Xg[[nm]]
    if (is.ordered(v)) {
      if (isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
        Xg[[nm]] <- factor(v, levels = levels(v), ordered = FALSE, exclude = NULL)
      } else {
        Xg[[nm]] <- as.numeric(v)
      }
      next
    }
    if (is.numeric(v) || is.integer(v)) next
    if (is.logical(v)) {
      Xg[[nm]] <- factor(v, levels = c(FALSE, TRUE))
      next
    }
    if (is.factor(v)) next
    if (is.matrix(v)) {
      Xg[[nm]] <- as.numeric(v)
      next
    }
    Xg[[nm]] <- factor(as.character(v))
  }
  
  fml <- as.formula(paste("~", paste(colnames(Xg), collapse = " + "), "-1"))
  tm <- terms(fml, data = Xg)
  MM <- model.matrix(tm, data = Xg)
  storage.mode(MM) <- "double"
  
  # Filter zero variance in design matrix
  ok <- apply(MM, 2L, function(col) {
    v <- stats::var(as.numeric(col), na.rm = TRUE)
    is.finite(v) && v > 1e-12
  })
  
  if (!any(ok)) stop("[design_with_map] all encoded columns were ~zero-variance")
  
  assign <- attr(MM, "assign")
  tl <- attr(tm, "term.labels")
  varmap <- tl[assign]
  MM <- MM[, ok, drop = FALSE]
  attr(MM, "varmap") <- varmap[ok]
  MM
}

residualise_foldsafe <- function(Xenc, Base, folds, k_gam = 6) {
  n <- nrow(Base)
  V <- colnames(Xenc)
  if (!length(V)) {
    return(matrix(numeric(0), nrow = n, ncol = 0L, dimnames = list(rownames(Base), character(0))))
  }

  Xenc_mat <- as.matrix(Xenc)
  storage.mode(Xenc_mat) <- "double"
  sm_terms <- paste0("s(b", seq_len(ncol(Base)), ",k=", k_gam, ")")
  fml_text <- paste("v ~", paste(sm_terms, collapse = " + "))
  fold_levels <- sort(unique(as.integer(folds)))
  fold_sets <- lapply(fold_levels, function(k) {
    list(tr = which(folds != k), te = which(folds == k))
  })
  base_frames <- lapply(fold_sets, function(ix) {
    list(
      tr = data.frame(Base[ix$tr, , drop = FALSE]),
      te = data.frame(Base[ix$te, , drop = FALSE])
    )
  })
  n_workers <- tryCatch(as.integer(future::nbrOfWorkers()), error = function(e) 1L)
  worker_parent <- new.env(parent = globalenv())
  residualise_worker_core <- residualise_one_col_worker
  environment(residualise_worker_core) <- worker_parent
  residualise_worker_env <- list2env(
    list(
      Xenc_mat = Xenc_mat,
      n = n,
      fold_sets = fold_sets,
      base_frames = base_frames,
      fml_text = fml_text,
      residualise_worker_core = residualise_worker_core
    ),
    parent = worker_parent
  )
  residualise_task <- function(i) {
    residualise_worker_core(i, Xenc_mat, n, fold_sets, base_frames, fml_text)
  }
  environment(residualise_task) <- residualise_worker_env

  cat(sprintf(
    "[Residuals] OOF GAM residualisation over %d encoded columns | folds=%d | workers=%d\n",
    length(V), length(fold_levels), n_workers
  ))

  if (n_workers <= 1L) {
    pb <- utils::txtProgressBar(min = 0, max = length(V), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
    cols <- lapply(seq_along(V), function(i) {
      out <- residualise_task(i)
      utils::setTxtProgressBar(pb, i)
      out
    })
  } else {
    cols <- FUTURE_LAPPLY(
      seq_along(V),
      residualise_task,
      future.packages = c("mgcv", "stats"),
      future.seed = TRUE,
      future.scheduling = 1
    )
  }

  E <- do.call(cbind, cols)
  colnames(E) <- V
  rownames(E) <- rownames(Base)
  E[, colSums(is.finite(E)) > 0, drop = FALSE]
}

residualise_one_col_worker <- function(i, Xenc_mat, n, fold_sets, base_frames, fml_text) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  options(mc.cores = 1)

  z <- as.numeric(Xenc_mat[, i])
  e <- rep(NA_real_, n)

  vz <- stats::var(z, na.rm = TRUE)
  if (!any(is.finite(z)) || !is.finite(vz) || vz <= 1e-12) {
    return(e)
  }

  fml_env <- list2env(list(s = mgcv::s), parent = baseenv())
  fml <- stats::as.formula(fml_text, env = fml_env)
  gam_ctrl <- mgcv::gam.control(nthreads = 1)

  for (j in seq_along(fold_sets)) {
    tr <- fold_sets[[j]]$tr
    te <- fold_sets[[j]]$te
    if (!length(tr) || !length(te)) next

    dftr <- base_frames[[j]]$tr
    dftr$v <- z[tr]

    g <- try(mgcv::gam(fml, data = dftr, method = "REML", control = gam_ctrl), silent = TRUE)
    if (inherits(g, "try-error")) next

    mu <- as.numeric(stats::predict(g, newdata = base_frames[[j]]$te, type = "response"))
    e[te] <- z[te] - mu
  }

  e
}

# --- Diagnostic Metrics Helpers ---

score_item_base <- function(nm, Z, varmap) {
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  
  if (length(idx) == 1L) return(as.numeric(Z[, idx]))
  
  # For multi-column items (factors), use first PC
  sc <- try(suppressWarnings(prcomp(Z[, idx, drop = FALSE], rank. = 1)$x[, 1]), silent = TRUE)
  if (inherits(sc, "try-error")) return(rep(NA_real_, nrow(Z)))
  
  as.numeric(sc)
}

e_from_E <- function(nm, E_scaled, Z, varmap) {
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  
  if (length(idx) == 1L) return(as.numeric(E_scaled[, idx]))
  
  # Project E_scaled columns onto the item's original structure
  pc1 <- try(prcomp(Z[, idx, drop = FALSE], rank. = 1), silent = TRUE)
  if (inherits(pc1, "try-error")) return(rep(NA_real_, nrow(Z)))
  
  as.numeric(as.matrix(E_scaled[, idx, drop = FALSE]) %*% pc1$rotation[, 1])
}

r2_base_linear <- function(Base, v) {
  if (nrow(Base) != length(v)) return(NA_real_)
  fit <- try(lm(v ~ Base[, 1] + Base[, 2]), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  summary(fit)$r.squared
}

r2_residual_cv <- function(e, nb, folds = CV_FOLDS, seed = SEED_GLOBAL, key = NULL) {
  n <- length(e)
  if (n < 10 || var(e, na.rm = TRUE) == 0 || all(is.na(e))) return(NA_real_)

  if (length(folds) > 1L) {
    fold_id <- as.integer(folds)
    if (length(fold_id) != n) stop("fold_id length mismatch: ", length(fold_id), " vs n=", n)
    fold_id[!is.finite(fold_id)] <- 1L
  } else {
    K <- as.integer(folds[1])
    if (!is.finite(K) || K < 2L) K <- 2L
    if (K > n) K <- n

    seed_eff <- if (is.null(key)) as.integer(seed) else .seed_from_key(seed, paste0("r2_resid|", key))
    fold_id <- .with_seed(seed_eff, sample(rep_len(seq_len(K), n)))
  }

  pred <- rep(NA_real_, n)
  
  for (f in sort(unique(fold_id))) {
    te <- which(fold_id == f)
    if (length(te) == 0) next
    
    # KNN imputation prediction on the test fold
    e_mask <- e
    e_mask[te] <- NA 
    nb_vals <- matrix(e_mask[nb], nrow = n) 
    pr <- rowMeans(nb_vals, na.rm = TRUE)
    
    pr[is.na(pr)] <- mean(e[-te], na.rm = TRUE)
    pred[te] <- pr[te]
  }
  
  ve <- stats::var(e, na.rm = TRUE)
  mse <- mean((e - pred)^2, na.rm = TRUE)
  
  if (!is.finite(ve) || ve <= 1e-9) return(0)
  max(0, 1 - mse / ve)
}

# --- Parallel Worker for Interactions ---

interact_worker <- function(i, grid, Zmat, BaseDF, DxDF, map) {
  
  # Prevent over-subscription: Lock internal threads to 1 per worker.
  # This avoids the "exploding threads" issue when parallelizing high-level tasks.
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  options(mc.cores = 1)
  
  dx_name <- grid$dx[i]
  nm      <- grid$var[i]
  
  y <- as.integer(DxDF[[dx_name]] > 0)
  v <- score_item_base(nm, Zmat, map) 
  
  if (all(is.na(v)) || stats::sd(v, na.rm = TRUE) < 1e-6) return(NULL)
  
  # Run OOF GAMs
  res <- oof_R2_two_gams(
    v, BaseDF, y,
    K_target = 5,
    k_gam = 10,
    seed = SEED_GLOBAL,
    seed_key = paste0(dx_name, "::", nm)
  )
  
  if (is.finite(unname(res["dR2"]))) {
    return(data.frame(
      dx = dx_name, var = nm,
      R2_add = unname(res["R2_add"]), R2_int = unname(res["R2_int"]),
      dR2 = unname(res["dR2"]), p_like = unname(res["p_like"]),
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# ==============================================================================
# 2. Data Ingest
# ==============================================================================

df <- readr::read_csv2(
  PSY_CSV,
  progress = FALSE, show_col_types = FALSE
)

id_col <- if ("participant_id" %in% names(df)) "participant_id" else names(df)[1]
ids_all <- as.character(df[[id_col]])

# Clean columns
drop_pattern <- "(?i)diagnosis|SCID|NODIAG"
cols_to_drop <- grep(drop_pattern, names(df), value = TRUE)

X <- dplyr::select(df, -dplyr::all_of(c(id_col, cols_to_drop))) |>
  as.data.frame(stringsAsFactors = FALSE)

# Robust type coercion
is_small_int_scale <- function(v) {
  vn <- suppressWarnings(as.numeric(v))
  if (all(is.na(vn))) return(FALSE)
  u <- sort(unique(na.omit(vn)))
  k <- length(u)
  k >= 3 && k <= 7 && all(abs(u - round(u)) < 1e-8)
}

for (nm in names(X)) {
  v <- X[[nm]]
  if (is.character(v)) {
    vn <- suppressWarnings(as.numeric(v))
    if (!all(is.na(vn))) v <- vn
  }
  if (all(v %in% c(0, 1, NA))) {
    if (isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
      X[[nm]] <- factor(v, ordered = FALSE, exclude = if (isTRUE(MISSING_AS_NOMINAL_LEVEL)) NULL else NA)
    } else {
      X[[nm]] <- factor(v, levels = c(0, 1), ordered = TRUE)
    }
  } else if (is_small_int_scale(v)) {
    if (isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
      X[[nm]] <- factor(
        as.integer(round(as.numeric(v))),
        ordered = FALSE,
        exclude = if (isTRUE(MISSING_AS_NOMINAL_LEVEL)) NULL else NA
      )
    } else {
      X[[nm]] <- factor(as.integer(round(as.numeric(v))), ordered = TRUE)
    }
  } else if (is.numeric(v)) {
    X[[nm]] <- as.numeric(v)
  } else {
    X[[nm]] <- factor(v)
  }
}

keep <- stats::complete.cases(X) & !is.na(ids_all) & ids_all != ""
X <- X[keep, , drop = FALSE]
ids_all <- ids_all[keep]
rownames(X) <- make.unique(ids_all)

degenerate_value <- constant_profile_value(X)
degenerate_mask <- !is.na(degenerate_value)
if (any(degenerate_mask)) {
  degenerate_constant_profiles <- data.frame(
    row = which(degenerate_mask),
    participant_id = rownames(X)[degenerate_mask],
    constant_value = degenerate_value[degenerate_mask],
    stringsAsFactors = FALSE
  )

  if (isTRUE(WRITE_DEGENERATE_CSV)) {
    write_csv(degenerate_constant_profiles, "degenerate_constant_profiles.csv")
  }

  deg_counts <- sort(table(degenerate_constant_profiles$constant_value), decreasing = TRUE)
  top_counts <- head(deg_counts, 5L)
  warning(sprintf(
    "[data] Found %d fully-constant profiles; kept in analysis%s. Top values: %s",
    nrow(degenerate_constant_profiles),
    if (isTRUE(WRITE_DEGENERATE_CSV)) " (see degenerate_constant_profiles.csv)" else "",
    paste(sprintf("%s=%d", names(top_counts), as.integer(top_counts)), collapse = ", ")
  ))
} else {
  degenerate_constant_profiles <- NULL
}

# --- Diagnoses Table Ingest (OPTIONAL) ---
DX_AVAILABLE <- FALSE
diag_wide_full <- NULL

safe_ingest_dx <- function(path, ids_here, optional = TRUE) {
  if (!nzchar(path) || !file.exists(path)) {
    if (optional) {
      message("[dx] Diagnoses file missing; proceeding without diagnoses.")
      return(list(
        available = FALSE,
        diag_wide_full = data.frame(participant_id = ids_here, stringsAsFactors = FALSE)
      ))
    } else {
      stop("Diagnoses file not found: ", path)
    }
  }
  
  dx <- try(
    readr::read_csv2(
      path,
      col_types = readr::cols(),
      progress = FALSE, show_col_types = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(dx, "try-error")) {
    if (optional) {
      warning("[dx] Read failed; disabling diagnoses (optional).")
      return(list(
        available = FALSE,
        diag_wide_full = data.frame(participant_id = ids_here, stringsAsFactors = FALSE)
      ))
    } else {
      stop(dx)
    }
  }
  
  dx <- dx %>% dplyr::mutate(participant_id = as.character(participant_id))
  
  mm <- match(ids_here, as.character(dx$participant_id))
  if (all(is.na(mm))) {
    if (optional) {
      warning("[dx] Join failed for 100% of rows; disabling diagnoses (optional).")
      return(list(
        available = FALSE,
        diag_wide_full = data.frame(participant_id = ids_here, stringsAsFactors = FALSE)
      ))
    } else {
      stop(sprintf("Diagnoses join failed for %d/%d rows.", sum(is.na(mm)), length(mm)))
    }
  }
  
  dxA <- dx[mm, , drop = FALSE]
  dxA$participant_id <- ids_here
  if (ncol(dxA) > 1) {
    dxA <- dxA %>%
      dplyr::mutate(dplyr::across(-participant_id, ~ {
        x <- suppressWarnings(as.integer(.))
        x[is.na(x)] <- 0L
        pmin(pmax(x, 0L), 1L)
      }))
  }
  
  avail <- (ncol(dxA) > 1) && any(colSums(dxA[, -1, drop = FALSE] > 0, na.rm = TRUE) > 0)
  list(available = avail, diag_wide_full = dxA)
}

ids_here <- rownames(X)
dx_ing   <- safe_ingest_dx(DIAG_CSV, ids_here, optional = isTRUE(DX_OPTIONAL))

DX_AVAILABLE   <- isTRUE(dx_ing$available)
diag_wide_full <- dx_ing$diag_wide_full

cat(sprintf("[ingest] X rows=%d, cols=%d | DX_available=%s | DX_cols=%d\n",
            nrow(X), ncol(X), DX_AVAILABLE, max(1L, ncol(diag_wide_full))))

if (DX_AVAILABLE && ncol(diag_wide_full) > 1) {
  DX_wide <- as.data.frame(diag_wide_full[, -1, drop = FALSE])
  rownames(DX_wide) <- diag_wide_full$participant_id
} else {
  DX_wide <- NULL
}

# Align diagnoses to X row order
ids_here <- rownames(X)
mm <- match(ids_here, as.character(diag_wide_full$participant_id))
if (anyNA(mm)) stop(sprintf("Diagnoses join failed for %d/%d rows.", sum(is.na(mm)), length(mm)))

DX_wide <- as.data.frame(diag_wide_full[mm, , drop = FALSE])
DX_wide$participant_id <- NULL

if (isTRUE(DX_DENY_NOS)) {
  nos_cols <- names(DX_wide)[tolower(names(DX_wide)) %in% "nodiag"]
  if (length(nos_cols)) DX_wide[nos_cols] <- NULL
}

cat(sprintf(
  "[ingest] X rows=%d, cols=%d | DX rows=%d, cols=%d (post NOS=%s)\n",
  nrow(X), ncol(X), nrow(DX_wide), ncol(DX_wide),
  if (isTRUE(DX_DENY_NOS)) "drop-nodiag" else "kept"
))

# ==============================================================================
# 3. Deduplication and Coreset Selection
# ==============================================================================

dedup_mode_eff <- if (!isTRUE(DO_DEDUP)) "none" else DEDUP_MODE

if (!dedup_mode_eff %in% c("gower_complete", "hash_exact", "hash_round", "none")) {
  stop("Unknown DEDUP_MODE: ", dedup_mode_eff)
}

# Helper: convert a single column to a stable string representation for hashing.
# Hash modes always use the prepped, non-jittered representation so that ties are
# preserved and the key is deterministic.
dedup_key_col <- function(v, mode = c("hash_exact", "hash_round"),
                          digits = 6L, na_token = "<NA>") {
  mode <- match.arg(mode)
  
  if (is.numeric(v)) {
    x <- as.numeric(v)
    
    if (mode == "hash_exact") {
      out <- ifelse(is.na(x), na_token, sprintf("%.17g", x))
      out[out %in% c("-0", "-0.0")] <- "0"
      return(out)
    }
    
    x <- round(x, digits = digits)
    fmt <- paste0("%.", digits, "f")
    out <- ifelse(is.na(x), na_token, sprintf(fmt, x))
    
    # Normalise negative zero after rounding
    zero_tok <- sprintf(fmt, 0)
    out[grepl(paste0("^-0(?:\\.0+)?$"), out)] <- zero_tok
    return(out)
  }
  
  if (is.factor(v) || is.ordered(v) || is.character(v) || is.logical(v)) {
    out <- as.character(v)
    out[is.na(out)] <- na_token
    return(out)
  }
  
  out <- as.character(v)
  out[is.na(out)] <- na_token
  out
}

# Helper: build duplicate groups from a hashed row key
dedup_groups_from_hash <- function(Xdf, mode = c("hash_exact", "hash_round"),
                                   digits = 6L, na_token = "<NA>") {
  mode <- match.arg(mode)
  
  X_key <- as.data.frame(
    lapply(Xdf, dedup_key_col, mode = mode, digits = digits, na_token = na_token),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Use a separator that is highly unlikely to occur in real data
  key <- do.call(paste, c(X_key, sep = "\r"))
  split(seq_len(nrow(Xdf)), key, drop = TRUE)
}

if (dedup_mode_eff == "gower_complete") {
  
  PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE, seed = SEED_JITTER)
  X_for_id <- PX$X
  Dg <- gower_dist(X_for_id, type_list = PX$type, weights = PX$weights)
  
  # Auto-detect EPS_DEDUP using knee method if NA
  if (is.na(EPS_DEDUP)) {
    d1 <- first_nn_d1(Dg)
    qlo <- as.numeric(stats::quantile(d1[is.finite(d1)], probs = c(0.001, 0.01), na.rm = TRUE))
    eps_grid <- seq(from = max(0, min(qlo, na.rm = TRUE) - 0.02),
                    to = min(0.50, stats::quantile(d1, 0.30, na.rm = TRUE)),
                    by = 0.005)
    cc <- collapse_curve(Dg, eps_grid)
    dprop <- diff(cc$prop_retained) / diff(cc$eps)
    knee_i <- which.min(dprop)
    EPS_DEDUP <- mean(cc$eps[c(knee_i, knee_i + 1)])
  }
  
  gr_all <- complete_groups(Dg, EPS_DEDUP)
  
  # Identify medoids of complete-linkage groups
  Dm_g <- as.matrix(Dg)
  diag(Dm_g) <- 0
  split_idx <- split(seq_len(nrow(Dm_g)), gr_all)
  
  reps <- vapply(split_idx, function(ix) {
    ix[which.min(rowSums(Dm_g[ix, ix, drop = FALSE]))]
  }, integer(1))
  
  mult <- as.integer(lengths(split_idx))
  
  Dg_rep <- stats::as.dist(Dm_g[reps, reps, drop = FALSE])
  core_idx_rep <- twonn_core_by_slope(Dg_rep)
  
  cat(sprintf("[Dedup:gower_complete] eps=%.3f | reps=%d of %d | core_rep=%d\n",
              EPS_DEDUP, length(reps), nrow(X), length(core_idx_rep)))
  
  if (isTRUE(WRITE_DEDUP_CSV)) {
    mult_df <- data.frame(
      rep_row = reps,
      representative_id = rownames(X_for_id)[reps],
      multiplicity = mult
    )
    write_csv(mult_df, sprintf("near_duplicate_groups_complete_eps%g.csv", EPS_DEDUP))
  }
  
} else if (dedup_mode_eff %in% c("hash_exact", "hash_round")) {
  
  # For hash dedup we deliberately disable jitter; otherwise every row becomes unique.
  PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = FALSE, seed = SEED_JITTER)
  X_for_id <- PX$X
  
  split_idx <- dedup_groups_from_hash(
    X_for_id,
    mode = dedup_mode_eff,
    digits = DEDUP_HASH_DIGITS,
    na_token = DEDUP_HASH_NA
  )
  
  # Hash groups do not have an internal distance geometry, so keep the first row
  # in each group as the deterministic representative.
  reps <- vapply(split_idx, `[`, integer(1), 1L)
  mult <- as.integer(lengths(split_idx))
  
  # core_idx_rep is currently passed downstream but not actually used by the optimiser.
  # In hash mode, keep all representatives.
  core_idx_rep <- seq_along(reps)
  
  cat(sprintf("[Dedup:%s] reps=%d of %d | core_rep=%d | hash_digits=%d\n",
              dedup_mode_eff, length(reps), nrow(X), length(core_idx_rep), DEDUP_HASH_DIGITS))
  
  if (isTRUE(WRITE_DEDUP_CSV)) {
    mult_df <- data.frame(
      rep_row = reps,
      representative_id = rownames(X_for_id)[reps],
      multiplicity = mult
    )
    suffix <- if (dedup_mode_eff == "hash_round") {
      paste0("hash_round_d", DEDUP_HASH_DIGITS)
    } else {
      "hash_exact"
    }
    write_csv(mult_df, sprintf("near_duplicate_groups_%s.csv", suffix))
  }
  
} else {
  
  reps <- seq_len(nrow(X))
  mult <- rep(1L, nrow(X))
  core_idx_rep <- seq_along(reps)
  
  cat(sprintf("[Dedup:none] reps=%d of %d | core_rep=%d\n",
              length(reps), nrow(X), length(core_idx_rep)))
  
  if (isTRUE(WRITE_DEDUP_CSV)) {
    mult_df <- data.frame(
      rep_row = reps,
      representative_id = rownames(X)[reps],
      multiplicity = mult
    )
    write_csv(mult_df, "near_duplicate_groups_none.csv")
  }
}

# ==============================================================================
# 4. Gower Weight Optimization
# ==============================================================================

# Drop constant columns and optionally trim high-correlation pairs
drop_constant_cols <- function(X) {
  keep <- vapply(X, function(v) length(unique(na.omit(v))) >= 1L, logical(1))
  X[, keep, drop = FALSE]
}
X_pred <- drop_constant_cols(X)

if (DO_CORR_TRIM) {
  corr_trim <- function(X, thr = 0.95) {
    num <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
    if (!ncol(num)) return(X)
    
    C <- suppressWarnings(cor(num, use = "pairwise.complete.obs"))
    drop <- character(0)
    for (i in seq_len(ncol(C) - 1)) {
      if (colnames(C)[i] %in% drop) next
      j <- which(abs(C[i, (i + 1):ncol(C)]) >= thr) + i
      drop <- union(drop, colnames(C)[j])
    }
    keep <- setdiff(colnames(X), drop)
    X[, keep, drop = FALSE]
  }
  X_pred <- corr_trim(X_pred, CORR_THRESH)
}
cat(sprintf("Start: X_pred has %d columns after constant-drop%s.\n\n", ncol(X_pred), if (DO_CORR_TRIM) " + corr-trim" else ""))

WEIGHTING_MODE <- tolower(WEIGHTING_MODE)
if (!WEIGHTING_MODE %in% c("id_guided", "uniform")) {
  stop("Unknown WEIGHTING_MODE: ", WEIGHTING_MODE)
}

w_init <- setNames(rep(1, ncol(X)), colnames(X))
allow <- setNames(rep(FALSE, ncol(X)), colnames(X))
allow[colnames(X_pred)] <- TRUE
gower_kmin <- max(30L, ceiling(0.10 * length(w_init)))

# --- Gower optimisation: reference run + optional stability replicates ---
wopt_list <- NULL
wopt <- NULL

if (identical(WEIGHTING_MODE, "id_guided")) {
  if (GOWER_MULTI_ENABLE && GOWER_MULTI_RUNS > 1L) {
    cat(sprintf("[Gower-multi] Running %d optimisation replicates (reference run = 1, base seed = %d).\n",
                GOWER_MULTI_RUNS, SEED_GLOBAL))

    run_ids <- seq_len(GOWER_MULTI_RUNS)
    worker_parent <- new.env(parent = globalenv())

    hash32_core <- .hash32
    environment(hash32_core) <- worker_parent

    seed_from_key_core <- .seed_from_key
    environment(seed_from_key_core) <- list2env(
      list(.hash32 = hash32_core),
      parent = worker_parent
    )

    with_seed_core <- .with_seed
    environment(with_seed_core) <- worker_parent

    two_nn_from_distvec_core <- .two_nn_from_distvec
    environment(two_nn_from_distvec_core) <- worker_parent

    twonn_id_from_dist_core <- twonn_id_from_dist
    environment(twonn_id_from_dist_core) <- list2env(
      list(.two_nn_from_distvec = two_nn_from_distvec_core),
      parent = worker_parent
    )

    prep_X_for_gower_core <- prep_X_for_gower
    environment(prep_X_for_gower_core) <- list2env(
      list(
        .seed_from_key = seed_from_key_core,
        .with_seed = with_seed_core
      ),
      parent = worker_parent
    )

    make_NS_cache_core <- make_NS_cache
    environment(make_NS_cache_core) <- worker_parent

    optimise_gower_weights_constrained_core <- optimise_gower_weights_constrained
    environment(optimise_gower_weights_constrained_core) <- list2env(
      list(
        .seed_from_key = seed_from_key_core,
        .with_seed = with_seed_core,
        prep_X_for_gower = prep_X_for_gower_core,
        make_NS_cache = make_NS_cache_core,
        twonn_id_from_dist = twonn_id_from_dist_core,
        RARE_LEVEL_MIN_PROP = RARE_LEVEL_MIN_PROP,
        FIX_REP_SUBSET = FIX_REP_SUBSET
      ),
      parent = worker_parent
    )

    wopt_list <- progressr::with_progress({
      p <- progressr::progressor(steps = length(run_ids))
      gower_multi_run_task_env <- list2env(
        list(
          .seed_from_key = seed_from_key_core,
          X = X,
          w_init = w_init,
          allow = allow,
          SEED_GLOBAL = SEED_GLOBAL,
          SEED_JITTER = SEED_JITTER,
          W_MIN = W_MIN,
          W_STEP_GRID = W_STEP_GRID,
          W_BATCH_K = W_BATCH_K,
          W_BATCH_FACTOR = W_BATCH_FACTOR,
          W_MAX_ITERS = W_MAX_ITERS,
          N_ROWS_SUB = N_ROWS_SUB,
          NCORES_PAR = NCORES_PAR,
          reps = reps,
          core_idx_rep = core_idx_rep,
          GOWER_MULTI_RUNS = GOWER_MULTI_RUNS,
          optimise_gower_weights_constrained_core = optimise_gower_weights_constrained_core,
          p = p
        ),
        parent = worker_parent
      )
      gower_multi_run_task <- function(r) {
        base_seed_r <- if (r == 1L) SEED_GLOBAL else .seed_from_key(SEED_GLOBAL, paste0("gower_run_", r))
        jitter_r <- if (r == 1L) SEED_JITTER else .seed_from_key(SEED_JITTER, paste0("gower_run_", r))
        verbose_r <- (r == 1L)
        plot_r <- FALSE
        progress_cb <- if (r == 1L) {
          function(info) {
            cat(sprintf("[Gower ref] iter %d: ID=%.3f\n", info$iter, info$ID))
          }
        } else {
          NULL
        }

        out <- optimise_gower_weights_constrained_core(
          X,
          init_weights = w_init,
          allow_update = allow,
          objective = "TwoNN_all",
          w_min = W_MIN,
          step_grid = W_STEP_GRID,
          batch_k = W_BATCH_K,
          batch_factor = W_BATCH_FACTOR,
          max_iter = W_MAX_ITERS,
          n_rows_sub = N_ROWS_SUB,
          ncores = NCORES_PAR,
          seed_jitter = jitter_r,
          reps_idx = reps,
          core_idx_rep = core_idx_rep,
          base_seed = base_seed_r,
          verbose = verbose_r,
          plot_progress = plot_r,
          progress_fun = progress_cb
        )

        p(sprintf("replicate %d/%d", r, GOWER_MULTI_RUNS))
        out
      }
      environment(gower_multi_run_task) <- gower_multi_run_task_env

      FUTURE_LAPPLY(
        run_ids,
        gower_multi_run_task,
        future.packages = "progressr",
        future.seed = TRUE
      )
    })

    wopt <- wopt_list[[1L]]
  } else {
    wopt <- optimise_gower_weights_constrained(
      X,
      init_weights = w_init,
      allow_update = allow,
      objective = "TwoNN_all",
      w_min = W_MIN,
      step_grid = W_STEP_GRID,
      batch_k = W_BATCH_K,
      batch_factor = W_BATCH_FACTOR,
      max_iter = W_MAX_ITERS,
      n_rows_sub = N_ROWS_SUB,
      ncores = NCORES_PAR,
      seed_jitter = SEED_JITTER,
      reps_idx = reps,
      core_idx_rep = core_idx_rep,
      base_seed = SEED_GLOBAL,
      verbose = TRUE,
      plot_progress = TRUE
    )
  }

  ix_pre <- wopt$idx_used
  PX_pre <- prep_X_for_gower(
    X[ix_pre, , drop = FALSE],
    rare_prop = RARE_LEVEL_MIN_PROP,
    do_jitter = TRUE,
    seed = SEED_JITTER
  )
  Xsub_pre <- PX_pre$X
  typ_pre <- PX_pre$type

  D_un_pre <- gower_dist(Xsub_pre, type_list = typ_pre, weights = rep(1, ncol(Xsub_pre)))
  w_pre <- wopt$weights[colnames(Xsub_pre)]
  D_op_pre <- gower_dist(Xsub_pre, type_list = typ_pre, weights = w_pre)

  dist_health_console(D_un_pre, "unweighted PRE")
  dist_health_console(D_op_pre, "ID-guided PRE")

  if (!is.null(wopt_list) && length(wopt_list) > 1L) {
    var_names <- names(wopt$weights)
    n_runs <- length(wopt_list)

    surv_mat <- matrix(0L, nrow = length(var_names), ncol = n_runs,
                       dimnames = list(var_names, paste0("run", seq_len(n_runs))))
    sel_mat <- matrix(0L, nrow = length(var_names), ncol = n_runs,
                      dimnames = list(var_names, paste0("run", seq_len(n_runs))))
    final_ID_vec <- rep(NA_real_, n_runs)
    thr_tail_vec <- rep(NA_real_, n_runs)

    sel_ref <- survivors_from_weights(
      w = wopt$weights,
      w_min = W_MIN,
      kmin = gower_kmin,
      make_plot = TRUE
    )
    sel_mat[sel_ref$survivors, 1L] <- 1L
    surv_mat[sel_ref$active_survivors, 1L] <- 1L
    final_ID_vec[1L] <- wopt$final_ID
    thr_tail_vec[1L] <- sel_ref$thr_tail

    if (n_runs > 1L) {
      for (r in 2:n_runs) {
        wr <- wopt_list[[r]]$weights
        sel_r <- survivors_from_weights(
          w = wr,
          w_min = W_MIN,
          kmin = gower_kmin,
          make_plot = FALSE
        )
        sel_mat[sel_r$survivors, r] <- 1L
        surv_mat[sel_r$active_survivors, r] <- 1L
        final_ID_vec[r] <- wopt_list[[r]]$final_ID
        thr_tail_vec[r] <- sel_r$thr_tail
      }
    }

    surv_count <- rowSums(surv_mat)
    surv_prop <- surv_count / n_runs
    sel_count <- rowSums(sel_mat)
    sel_prop <- sel_count / n_runs
    Wmat <- vapply(seq_len(n_runs), function(r) wopt_list[[r]]$weights[var_names], numeric(length(var_names)))
    dimnames(Wmat) <- list(var_names, paste0("run", seq_len(n_runs)))

    stability_tbl <- data.frame(
      var = var_names,
      count = as.integer(surv_count),
      prop = as.numeric(surv_prop),
      count_selected = as.integer(sel_count),
      prop_selected = as.numeric(sel_prop),
      selected_ref = var_names %in% sel_ref$survivors,
      in_ref = var_names %in% sel_ref$active_survivors,
      weight_ref = as.numeric(wopt$weights[var_names]),
      weight_mean = rowMeans(Wmat, na.rm = TRUE),
      weight_sd = apply(Wmat, 1L, stats::sd, na.rm = TRUE),
      weight_min = apply(Wmat, 1L, min, na.rm = TRUE),
      weight_max = apply(Wmat, 1L, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    stability_tbl <- stability_tbl[order(-stability_tbl$prop, -stability_tbl$weight_mean), ]

    write_csv(stability_tbl, "gower_survivor_stability.csv")
    write_csv(
      data.frame(var = rownames(Wmat), Wmat, check.names = FALSE, stringsAsFactors = FALSE),
      "gower_weights_by_run_wide.csv"
    )

    active_floor <- W_MIN + sel_ref$active_eps
    robust_mask <- stability_tbl$prop >= GOWER_MULTI_MIN_PROP & stability_tbl$weight_max > active_floor
    if (!any(robust_mask)) {
      warning(sprintf(
        "[Gower-multi] No active variable met prop >= %.2f; falling back to active reference survivors.",
        GOWER_MULTI_MIN_PROP
      ))
      robust_set <- sel_ref$active_survivors
      if (!length(robust_set)) robust_set <- sel_ref$survivors
      robust_mask <- stability_tbl$var %in% robust_set
    } else {
      robust_set <- stability_tbl$var[robust_mask]
    }

    cat(sprintf(
      "[Gower-multi] Reference selected = %d | reference active = %d | robust survivors = %d (min_prop=%.2f, runs=%d).\n",
      length(sel_ref$survivors), length(sel_ref$active_survivors), length(robust_set), GOWER_MULTI_MIN_PROP, n_runs
    ))

    w_multi_raw <- with(stability_tbl, weight_mean - W_MIN * (1 - prop))
    w_multi_raw[!is.finite(w_multi_raw)] <- 0
    w_multi_raw <- pmax(0, w_multi_raw)
    scale_fac <- max(w_multi_raw[robust_mask], na.rm = TRUE)
    if (!is.finite(scale_fac) || scale_fac <= 0) scale_fac <- 1

    w_full <- setNames(as.numeric(w_multi_raw / scale_fac), stability_tbl$var)
    survivors <- robust_set

    saveRDS(list(
      n_runs = n_runs,
      ref_run = 1L,
      min_prop = GOWER_MULTI_MIN_PROP,
      final_ID = final_ID_vec,
      thr_tail = thr_tail_vec,
      survivors_ref = sel_ref$survivors,
      survivors_ref_active = sel_ref$active_survivors,
      survivors_robust = robust_set,
      weights = lapply(wopt_list, function(z) z$weights)
    ), file = "gower_multi_runs_summary.rds")
  } else {
    sel <- survivors_from_weights(
      w = wopt$weights,
      w_min = W_MIN,
      kmin = gower_kmin,
      make_plot = TRUE
    )
    w_full <- sel$w_clamped
    survivors <- sel$survivors
  }
} else {
  survivors <- colnames(X_pred)
  if (!length(survivors)) {
    stop("[weights] uniform mode retained zero predictors after preprocessing.")
  }

  w_full <- setNames(rep(0, ncol(X)), colnames(X))
  w_full[survivors] <- 1
  wopt <- list(
    idx_used = reps,
    history = data.frame(),
    final_ID = NA_real_,
    weights = w_full
  )

  cat(sprintf(
    "[weights] Uniform mode: skipping Gower optimisation; retaining %d/%d predictors with unit weights.\n",
    length(survivors), ncol(X)
  ))
}

X <- X[, survivors, drop = FALSE]
w_all <- w_full[survivors]
cat(sprintf("[weights] Mode=%s | retained predictors: p=%d\n", WEIGHTING_MODE, ncol(X)))

# Recalculate Distances with final weights on a diagnostic subset only
PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE, seed = SEED_JITTER)
Xg <- PX$X
type_g <- PX$type
w_use <- w_all[colnames(Xg)]

if (!FINAL_DIAG_MODE %in% c("full", "reps", "sample_reps", "sample_all", "none")) {
  stop("Unknown FINAL_DIAG_MODE: ", FINAL_DIAG_MODE)
}

# Build diagnostic index set
diag_idx <- switch(
  FINAL_DIAG_MODE,
  
  "none" = integer(0),
  
  "full" = seq_len(nrow(Xg)),
  
  "reps" = {
    if (exists("reps") && length(reps) > 0L) reps else seq_len(nrow(Xg))
  },
  
  "sample_reps" = {
    pool <- if (exists("reps") && length(reps) > 0L) reps else seq_len(nrow(Xg))
    n_take <- min(length(pool), DIAG_N_MAX)
    if (length(pool) > n_take) {
      warn_diag_subsample("Final diag", n_take, length(pool), DIAG_N_MAX, "representatives")
      .with_seed(.seed_from_key(SEED_GLOBAL, "final_diag_sample_reps"), sample(pool, n_take))
    } else {
      pool
    }
  },
  
  "sample_all" = {
    pool <- seq_len(nrow(Xg))
    n_take <- min(length(pool), DIAG_N_MAX)
    if (length(pool) > n_take) {
      warn_diag_subsample("Final diag", n_take, length(pool), DIAG_N_MAX, "rows")
      .with_seed(.seed_from_key(SEED_GLOBAL, "final_diag_sample_all"), sample(pool, n_take))
    } else {
      pool
    }
  }
)

if (length(diag_idx) >= 3L) {
  Xg_diag <- Xg[diag_idx, , drop = FALSE]
  
  D_final <- cluster::daisy(Xg_diag, metric = "gower", type = type_g, weights = w_use)
  D_un_final <- cluster::daisy(Xg_diag, metric = "gower", type = type_g, weights = rep(1, ncol(Xg_diag)))
  ID_all <- twonn_id_from_dist(D_final)
  core_ix <- twonn_core_by_slope(D_final)
  
  DmF <- as.matrix(D_final)
  diag(DmF) <- Inf
  
  ID_core <- if (length(core_ix) >= 3L) {
    twonn_id_from_dist(as.dist(DmF[core_ix, core_ix, drop = FALSE]))
  } else {
    NA_real_
  }
  
  ID_LB <- if (length(core_ix) >= 15L) {
    lb_mle_id(DmF[core_ix, core_ix, drop = FALSE], 5, 15)
  } else {
    NA_real_
  }
  
  cat(sprintf(
    "[Constrained %s | diag=%s | n=%d] TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%d, p_active=%d)\n",
    toupper(PREF_TARGET), FINAL_DIAG_MODE, nrow(Xg_diag),
    ID_all, ID_core, ID_LB, length(core_ix), ncol(Xg)
  ))

  dist_health_console(D_un_final, sprintf("unweighted FINAL (%s)", FINAL_DIAG_MODE))
  dist_health_console(D_final, sprintf("weighted FINAL (%s)", FINAL_DIAG_MODE))
  for (k in c(10, 15, 30)) {
    nn_un <- knn_from_dist(D_un_final, k = k)
    nn_op <- knn_from_dist(D_final, k = k)
    cat(sprintf("[kNN overlap FINAL] k=%d | mean Jaccard = %.3f\n", k, mean_jaccard_knn(nn_un, nn_op)))
  }
} else {
  cat(sprintf(
    "[Constrained %s | diag=%s] Final distance-based diagnostics skipped.\n",
    toupper(PREF_TARGET), FINAL_DIAG_MODE
  ))
}

# ==============================================================================
# 5. PCA, Whitening, and Residuals
# ==============================================================================

# Expand design matrix
Xenc <- design_with_map(X)
varmap <- attr(Xenc, "varmap")
vars <- unique(varmap)

# Distribute Gower weights to encodings; sqrt-scale so PCA respects them
w_enc <- setNames(rep(1, ncol(Xenc)), colnames(Xenc))
alloc <- table(varmap)
for (nm in names(alloc)) {
  idx <- which(varmap == nm)
  wj <- w_all[nm]
  if (!is.finite(wj)) wj <- 1
  w_enc[idx] <- wj / length(idx)
}
Xenc_w <- sweep(Xenc, 2, sqrt(pmax(w_enc, 0)), "*")

m_star <- as.integer(M_STAR_FIXED)
Z <- scale(Xenc_w, center = TRUE, scale = TRUE)

BASE_DECOMP_METHOD <- tolower(BASE_DECOMP_METHOD)
if (!BASE_DECOMP_METHOD %in% c("robust_pca", "pca")) {
  stop("Unknown BASE_DECOMP_METHOD: ", BASE_DECOMP_METHOD)
}

k_eff <- max(2L, min(m_star, nrow(Z) - 1L, ncol(Z)))

if (identical(BASE_DECOMP_METHOD, "robust_pca")) {
  rp <- rrcov::PcaHubert(
    x        = Z,
    k        = k_eff,
    kmax     = k_eff,
    scale    = FALSE,
    signflip = TRUE
  )
  Base <- rp@scores[, 1:m_star, drop = FALSE]
  base_loadings <- rp@loadings[, 1:m_star, drop = FALSE]
  base_spectrum <- rp@eigenvalues[seq_len(k_eff)]
} else {
  pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE, rank. = k_eff)
  Base <- pc$x[, 1:m_star, drop = FALSE]
  base_loadings <- pc$rotation[, 1:m_star, drop = FALSE]
  base_spectrum <- pc$sdev[seq_len(k_eff)]^2
}

base_total_var <- sum(apply(Z, 2, stats::var))
base_explained_variance_ratio <- base_spectrum / pmax(base_total_var, 1e-12)

colnames(Base) <- paste0("b", seq_len(ncol(Base)))
colnames(base_loadings) <- paste0("b", seq_len(ncol(base_loadings)))

# Base[, 1] <- -Base[, 1]
# base_loadings[, 1] <- -base_loadings[, 1]
# Base[, 2] <- -Base[, 2]
# base_loadings[, 2] <- -base_loadings[, 2]

Base_A <- as.data.frame(Base[, 1:2, drop = FALSE])
colnames(Base_A) <- c("b1", "b2")
ids_base <- rownames(Base)
stopifnot(length(ids_base) == nrow(Base), all(nzchar(ids_base)))

cat(sprintf("[Base] Decomposition=%s. Proceeding to Whitening...\n", BASE_DECOMP_METHOD))

# --- Whitening for Neighbour Search ---
S <- stats::cov(Base)
U <- try(chol(S + diag(1e-8, ncol(Base))), silent = TRUE)
if (inherits(U, "try-error")) {
  eig <- eigen(S, symmetric = TRUE)
  U <- t(eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8))) %*% t(eig$vectors))
}
Base_w <- Base %*% solve(U)

KS_RESIDUAL <- c(6, 8, 10, 12, 16, 20)
nb_list <- setNames(lapply(KS_RESIDUAL, function(k) {
  RANN::nn2(Base_w, Base_w, k = pmin(k + 1L, nrow(Base_w)))$nn.idx[, -1L, drop = FALSE]
}), paste0("k", KS_RESIDUAL))

# --- Residual Extraction  ---
make_strat_folds <- function(y, K, seed = SEED_GLOBAL, key = NULL) {
  y <- as.integer(y)
  n <- length(y)
  K <- as.integer(K[1])
  if (!is.finite(K) || K < 2L) K <- 2L
  if (K > n) K <- n

  seed_eff <- if (is.null(key)) as.integer(seed) else .seed_from_key(seed, paste0("strat|", key))

  .with_seed(seed_eff, {
    folds <- integer(n)
    idx0 <- which(y == 0)
    idx1 <- which(y == 1)
    if (length(idx0)) folds[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
    if (length(idx1)) folds[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
    folds
  })
}

choose_K <- function(y, K_target = CV_FOLDS, min_per_class = 8) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  max(2, min(K_target, floor(n1 / min_per_class), floor(n0 / min_per_class)))
}

# Define Target for CV folds
if (!is.null(DX_wide) && "NODIAG" %in% names(DX_wide)) {
  y_use <- as.integer(DX_wide$NODIAG == 0L)
  stopifnot(length(y_use) == nrow(X))
  
  cat(sprintf("[target] from NODIAG | prev=%.3f | n0=%d (NODIAG=1) | n1=%d (any dx)\n",
              mean(y_use), sum(y_use == 0L), sum(y_use == 1L)))
  
  K_fold <- choose_K(y_use, K_target = CV_FOLDS, min_per_class = 8)
  fold_id <- make_strat_folds(y_use, K = K_fold, seed = SEED_GLOBAL, key = "target_folds")
} else {
  n <- nrow(X)
  K_fold <- max(2L, min(CV_FOLDS, n))
  fold_id <- .with_seed(
    .seed_from_key(SEED_GLOBAL, "target_unsupervised_folds"),
    sample(rep(seq_len(K_fold), length.out = n))
  )
  y_use <- rep(NA_integer_, n)
  
  message("[target] NODIAG unavailable; using unsupervised folds.")
}

E <- residualise_foldsafe(Xenc_w, Base, folds = fold_id, k_gam = 6)
E_scaled <- scale(E, center = TRUE, scale = TRUE)

# XR represents the "Residuals" (high-dimensional, post-OOF)
XR <- E_scaled
if (any(!is.finite(XR))) XR[!is.finite(XR)] <- 0

cat(sprintf("[Residuals] XR matrix: %d rows × %d columns (post-OOF, scaled)\n", nrow(XR), ncol(XR)))

# ==============================================================================
# 6. Item Diagnostics (Base vs Resid Roles)
# ==============================================================================

Z0_std <- scale(Xenc, center = TRUE, scale = TRUE)
vars_diag <- unique(varmap)

# Helper for KNN K selection
choose_k_nb <- function(e, nb_list, folds = CV_FOLDS, seed = SEED_GLOBAL, key = NULL) {
  folds_use <- folds
  if (length(folds) == 1L) {
    K <- as.integer(folds[1])
    if (!is.finite(K) || K < 2L) K <- 2L
    if (K > length(e)) K <- length(e)
    seed_eff <- if (is.null(key)) as.integer(seed) else .seed_from_key(seed, paste0("choose_knn|", key))
    folds_use <- .with_seed(seed_eff, sample(rep_len(seq_len(K), length(e))))
  }

  r2s <- vapply(nb_list, function(nb) {
    r2_residual_cv(e, nb, folds = folds_use, seed = seed)
  }, numeric(1))
  ix <- which.max(r2s)
  list(k = as.integer(sub("^k", "", names(nb_list)[ix])), R2_cv = r2s[ix])
}

RESIDUAL_PERM_B <- N_PERM
MIN_SD_ITEM <- 1e-6

cat(sprintf("[Item Roles] Starting diagnostics on %d items...\n", length(vars_diag)))

roles_rows <- lapply(vars_diag, function(nm) {
  # Score item (Original)
  v <- score_item_base(nm, Z0_std, varmap)
  if (all(is.na(v)) || stats::sd(v, na.rm = TRUE) < MIN_SD_ITEM) return(NULL)
  
  # Score residual
  e_item <- e_from_E(nm, E_scaled, Z0_std, varmap)
  if (all(is.na(e_item)) || stats::sd(e_item, na.rm = TRUE) < MIN_SD_ITEM) return(NULL)
  
  # Metrics
  R2_base <- r2_base_linear(Base, v)
  sel <- choose_k_nb(e_item, nb_list, folds = CV_FOLDS, seed = SEED_GLOBAL, key = nm)
  R2_resid <- sel$R2_cv
  
  data.frame(
    var = nm,
    R2_base = R2_base,
    R2_residual = R2_resid,
    k_residual = sel$k,
    stringsAsFactors = FALSE
  )
})

roles_df <- dplyr::bind_rows(roles_rows)

if (nrow(roles_df) > 0) {
  # Classify items
  ES_BASE <- 0.08      
  ES_RESIDUAL <- 0.05
  
  roles_df$role <- dplyr::case_when(
    roles_df$R2_base >= ES_BASE & roles_df$R2_residual >= ES_RESIDUAL ~ "mixed",
    roles_df$R2_base >= ES_BASE ~ "base-aligned",
    roles_df$R2_residual >= ES_RESIDUAL ~ "residual-structured",
    TRUE ~ "weak"
  )
  
  roles_df$role_final <- factor(roles_df$role, levels = c("base-aligned", "mixed", "residual-structured", "weak"))
  roles_df <- roles_df[order(-pmax(roles_df$R2_base, roles_df$R2_residual), roles_df$var), ]
  
  write_csv(roles_df, "predictive_item_roles_diagnostics.csv")
  cat(sprintf("[roles] wrote: %s  (p=%d items processed)\n", "predictive_item_roles_diagnostics.csv", nrow(roles_df)))
} else {
  warning("[roles] No valid items found for roles analysis.")
}

# ==============================================================================
# 7. Resid-Only Decomposition & Clustering
# ==============================================================================

stopifnot(exists("E_scaled"), is.matrix(E_scaled), nrow(E_scaled) >= 3) 
Ef <- E_scaled
Ef[!is.finite(Ef)] <- 0
n <- nrow(Ef)
pE <- ncol(Ef)
if (pE < 2L || n < 4L) stop("[Resid-only] insufficient columns/rows in E.")

RESIDUAL_BASE_MAX <- min(6L, pE, n - 1L)

pc_f <- prcomp(Ef, rank. = max(2L, RESIDUAL_BASE_MAX))
Bprime_all <- pc_f$x[, 1:RESIDUAL_BASE_MAX, drop = FALSE]
colnames(Bprime_all) <- paste0("f", seq_len(ncol(Bprime_all)))

pick_m_via_tc <- function(Xhigh,
                          Xlow_all,
                          ks = 10:30,
                          mmax = ncol(Xlow_all),
                          lambda = 0.02,
                          max_n = DIAG_N_MAX,
                          pool_idx = NULL,
                          seed = SEED_GLOBAL,
                          key = "pick_m_via_tc") {
  Xhigh <- as.matrix(Xhigh)
  Xlow_all <- as.matrix(Xlow_all)
  stopifnot(nrow(Xhigh) == nrow(Xlow_all))

  mmax <- min(as.integer(mmax), ncol(Xlow_all))
  if (!is.finite(mmax) || mmax < 2L) return(max(1L, mmax))

  n_full <- nrow(Xhigh)
  if (!is.finite(max_n) || max_n <= 0L) max_n <- n_full
  max_n <- max(as.integer(max_n), max(ks) + 2L, 3L)
  if (is.null(pool_idx) || !length(pool_idx)) {
    pool_idx <- seq_len(n_full)
  } else {
    pool_idx <- sort(unique(as.integer(pool_idx)))
    pool_idx <- pool_idx[is.finite(pool_idx) & pool_idx >= 1L & pool_idx <= n_full]
    if (!length(pool_idx)) pool_idx <- seq_len(n_full)
  }

  if (length(pool_idx) > max_n) {
    warn_diag_subsample(
      "pick_m_via_tc",
      max_n,
      length(pool_idx),
      max_n,
      if (length(pool_idx) < n_full) "representatives" else "rows"
    )
    idx <- .with_seed(.seed_from_key(seed, key), sample(pool_idx, max_n))
    Xhigh <- Xhigh[idx, , drop = FALSE]
    Xlow_all <- Xlow_all[idx, , drop = FALSE]
  } else if (length(pool_idx) < n_full) {
    idx <- pool_idx
    Xhigh <- Xhigh[idx, , drop = FALSE]
    Xlow_all <- Xlow_all[idx, , drop = FALSE]
  }

  rank_from_dist <- function(D) {
    n <- nrow(D)
    R <- matrix(0L, n, n)
    for (i in seq_len(n)) {
      ord <- order(D[i, ])
      R[i, ord] <- seq_len(n)
    }
    R
  }

  trust_cont <- function(high, low, ks = 10:30, Rh = NULL) {
    high <- as.matrix(high)
    low <- as.matrix(low)
    stopifnot(nrow(high) == nrow(low))
    n <- nrow(high)
    ks_use <- ks[is.finite(ks) & ks >= 1L & ks < n]
    if (!length(ks_use)) stop("No valid ks for trust/continuity at n=", n)

    if (is.null(Rh)) {
      Dh <- as.matrix(stats::dist(high))
      diag(Dh) <- Inf
      Rh <- rank_from_dist(Dh)
    }

    Dl <- as.matrix(stats::dist(low))
    diag(Dl) <- Inf
    Rl <- rank_from_dist(Dl)

    res <- lapply(ks_use, function(k) {
      H <- t(apply(Rh, 1, function(r) order(r)[1:k]))
      L <- t(apply(Rl, 1, function(r) order(r)[1:k]))
      Tsum <- 0
      Csum <- 0
      for (i in seq_len(n)) {
        U <- setdiff(L[i, ], H[i, ])
        if (length(U)) Tsum <- Tsum + sum(pmax(Rh[i, U] - k, 0))
        V <- setdiff(H[i, ], L[i, ])
        if (length(V)) Csum <- Csum + sum(pmax(Rl[i, V] - k, 0))
      }
      denom <- n * k * (2 * n - 3 * k - 1)
      if (!is.finite(denom) || denom <= 0) {
        return(data.frame(k = k, Trust = NA_real_, Continuity = NA_real_))
      }
      data.frame(k = k, Trust = 1 - (2 / denom) * Tsum, Continuity = 1 - (2 / denom) * Csum)
    })
    do.call(rbind, res)
  }

  Dh <- as.matrix(stats::dist(Xhigh))
  diag(Dh) <- Inf
  Rh <- rank_from_dist(Dh)

  vals <- lapply(2:mmax, function(m) {
    tc <- trust_cont(Xhigh, Xlow_all[, 1:m, drop = FALSE], ks, Rh = Rh)
    data.frame(m = m, T = mean(tc$Trust, na.rm = TRUE), C = mean(tc$Continuity, na.rm = TRUE))
  }) |> dplyr::bind_rows()

  vals$score <- with(vals, (T + C) - lambda * (m - min(m)))
  vals$m[which.max(vals$score)]
}

pick_m_pool <- if (exists("reps") && length(reps) > 0L) reps else NULL
m_f <- pick_m_via_tc(
  Ef,
  Bprime_all,
  ks = 10:30,
  mmax = ncol(Bprime_all),
  lambda = 0.02,
  max_n = DIAG_N_MAX,
  pool_idx = pick_m_pool
)
Bprime <- Bprime_all[, 1:m_f, drop = FALSE]
colnames(Bprime) <- paste0("f", seq_len(ncol(Bprime)))

# OOF residuals of E on B' -> F' (linear)
Kf <- max(2L, min(CV_FOLDS, n))
folds_f <- .with_seed(
  .seed_from_key(SEED_PRED, "fprime_oof_folds"),
  sample(rep(1:Kf, length.out = n))
)

residualise_linear_oof <- function(Y, X, folds) {
  n <- nrow(Y)
  p <- ncol(Y)
  R <- matrix(NA_real_, n, p, dimnames = list(rownames(Y), colnames(Y)))
  Xdf <- as.data.frame(X)
  
  for (j in seq_len(p)) {
    y <- Y[, j]
    if (all(!is.finite(y))) next
    for (k in sort(unique(folds))) {
      tr <- which(folds != k)
      te <- which(folds == k)
      fit <- try(stats::lm(y ~ ., data = cbind.data.frame(y = y[tr], Xdf[tr, , drop = FALSE])), silent = TRUE)
      if (inherits(fit, "try-error")) next
      mu <- as.numeric(predict(fit, newdata = Xdf[te, , drop = FALSE]))
      R[te, j] <- y[te] - mu
    }
  }
  R[, colSums(is.finite(R)) > 0, drop = FALSE]
}
Fprime <- residualise_linear_oof(Ef, Bprime, folds_f)

# ID diagnostics in residual spaces
resid_diag_from_reps <- exists("reps") && length(reps) > 0L
resid_diag_pool <- if (resid_diag_from_reps) reps else seq_len(nrow(Bprime))
resid_diag_n <- min(length(resid_diag_pool), DIAG_N_MAX)
resid_diag_idx <- if (length(resid_diag_pool) > resid_diag_n) {
  warn_diag_subsample(
    "Resid-only",
    resid_diag_n,
    length(resid_diag_pool),
    DIAG_N_MAX,
    if (resid_diag_from_reps) "representatives" else "rows"
  )
  .with_seed(.seed_from_key(SEED_GLOBAL, "resid_space_diag_sample"), sample(resid_diag_pool, resid_diag_n))
} else {
  resid_diag_pool
}

if (length(resid_diag_idx) < 3L) resid_diag_idx <- seq_len(min(nrow(Bprime), 3L))

Bprime_diag <- Bprime[resid_diag_idx, , drop = FALSE]
D_Bprime <- stats::dist(Bprime_diag)
MB <- as.matrix(D_Bprime)
diag(MB) <- Inf
core_B <- core_band_idx(D_Bprime, k = CORE_KNN_K, band = CORE_BAND)
ID_B_all <- twonn_id_from_dist(D_Bprime)
ID_B_core <- if (length(core_B) >= 3) twonn_id_from_dist(stats::as.dist(MB[core_B, core_B])) else NA_real_
ID_B_LB <- if (length(core_B) >= 3) lb_mle_id(MB[core_B, core_B, drop = FALSE], 5, 15) else NA_real_

if (ncol(Fprime) >= 2) {
  Fprime_diag <- Fprime[resid_diag_idx, , drop = FALSE]
  D_Fprime <- stats::dist(Fprime_diag)
  MFp <- as.matrix(D_Fprime)
  diag(MFp) <- Inf
  core_Fp <- core_band_idx(D_Fprime, k = CORE_KNN_K, band = CORE_BAND)
  ID_Fp_all <- twonn_id_from_dist(D_Fprime)
  ID_Fp_core <- if (length(core_Fp) >= 3) twonn_id_from_dist(stats::as.dist(MFp[core_Fp, core_Fp])) else NA_real_
  ID_Fp_LB <- if (length(core_Fp) >= 3) lb_mle_id(MFp[core_Fp, core_Fp, drop = FALSE], 5, 15) else NA_real_
} else {
  ID_Fp_all <- ID_Fp_core <- ID_Fp_LB <- NA_real_
}

cat(sprintf("[Resid-only] ID(B'): TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%d)\n",
            ID_B_all, ID_B_core, ID_B_LB, length(core_B)))
cat(sprintf("[Resid-only] ID(F'): TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%s)\n",
            ID_Fp_all, ID_Fp_core, ID_Fp_LB, if (exists("core_Fp")) length(core_Fp) else "NA"))

# Clustering in B' (Louvain on kNN graph)
knn_graph <- function(X, k = 15) {
  idx <- RANN::nn2(X, X, k = pmin(k + 1L, nrow(X)))$nn.idx[, -1L, drop = FALSE]
  i <- rep(seq_len(nrow(X)), each = ncol(idx))
  j <- as.vector(idx)
  g <- igraph::graph_from_edgelist(cbind(i, j), directed = FALSE)
  igraph::simplify(g)
}
gF <- knn_graph(Bprime, k = 15)
clF <- igraph::cluster_louvain(gF)$membership

saveRDS(list(
  Bprime = Bprime,
  Fprime = Fprime,
  m_f = m_f,
  diag_idx = resid_diag_idx,
  ID_Bprime = c(all = ID_B_all, core = ID_B_core, LB = ID_B_LB),
  ID_Fprime = c(all = ID_Fp_all, core = ID_Fp_core, LB = ID_Fp_LB),
  clusters = clF
), file = "residual_only_summary.rds")

cat("[Resid-only] wrote:", file.path(OUTPUTS_DIR, "residual_only_summary.rds"), "\n")

# ==============================================================================
# 8. Predictive Diagnostics (OOF & Interactions)
# ==============================================================================
if (DX_AVAILABLE) {
  # Align dx labels to Base
  ids_base <- rownames(Base)
  stopifnot(length(ids_base) == nrow(Base), all(nzchar(ids_base)))
  
  DxW <- diag_wide_full %>%
    dplyr::transmute(participant_id = 
                       trimws(as.character(participant_id)),
                     dplyr::across(-participant_id, ~ as.integer(.x))) %>%
    dplyr::right_join(tibble::tibble(participant_id = ids_base, .row = seq_along(ids_base)), by = "participant_id") %>%
    dplyr::arrange(.row) %>%
    dplyr::select(-participant_id, -.row)
  DxW[is.na(DxW)] <- 0L
  
  prev <- colMeans(DxW > 0, na.rm = TRUE)
  cases <- colSums(DxW > 0, na.rm = TRUE)
  keep_dx <- names(prev)[(prev >= DX_PREV_MIN) & (cases >= DX_CASES_MIN)]
  
  if (!length(keep_dx)) {
    warning("[dx] No diagnosis passes thresholds; skipping predictive diagnostics.")
    keep_dx <- intersect(names(prev), names(prev)[cases > 0 & (1 - prev) * nrow(DxW) > 0])
  }
  DxW_A <- as.data.frame(DxW)[, keep_dx, drop = FALSE]
  rownames(DxW_A) <- ids_base
  
  # --- Diagnosis probability fields over Base ---
  dir.create(paste0(OUTPUTS_DIR, "/base_prob_grids"), showWarnings = FALSE)
  
  grid_from_base <- function(Base, nx = 140, ny = 140, pad = 0.05, q = c(0.01, 0.99)) {
    rx <- quantile(Base[, 1], q, na.rm = TRUE)
    ry <- quantile(Base[, 2], q, na.rm = TRUE)
    wx <- diff(rx); wy <- diff(ry)
    xs <- seq(rx[1] - pad * wx, rx[2] + pad * wx, length.out = nx)
    ys <- seq(ry[1] - pad * wy, ry[2] + pad * wy, length.out = ny)
    as.matrix(expand.grid(b1 = xs, b2 = ys))
  }
  
  predict_dx_surface <- function(y, Base, gridXY, k_gam = 30) {
    df <- data.frame(y = as.integer(y), b1 = Base[, 1], b2 = Base[, 2])
    fit <- try(mgcv::gam(y ~ s(b1, b2, k = k_gam), data = df, family = binomial(), method = "REML"), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    
    pg <- try(as.numeric(predict(fit,
                                 newdata = data.frame(
                                   b1 = gridXY[, 1],
                                   b2 = gridXY[, 2]),
                                 type = "response")),
              silent = TRUE)
    
    if (inherits(pg, "try-error")) return(NULL)
    pmin(pmax(pg, 1e-6), 1 - 1e-6)
  }
  
  make_fname_safe <- function(x, max_len = 80) {
    s <- iconv(as.character(x), to = "ASCII//TRANSLIT")
    s <- tolower(s)
    s <- gsub("[^a-z0-9]+", "_", s)
    s <- gsub("^_+|_+$", "", s)
    s <- substr(s, 1, max_len)
    if (!nzchar(s)) s <- "dx"
    s
  }
  
  gridXY <- grid_from_base(Base, nx = as.integer(BF_NGR_UV), ny = as.integer(BF_NGR_UV))
  dx_surface_files <- c()
  
  for (dx in names(DxW_A)) {
    y <- as.integer(DxW_A[[dx]] > 0)
    pg <- predict_dx_surface(y, Base, gridXY, k_gam = 30)
    if (is.null(pg)) {
      warning(sprintf("[base-field] skip dx=%s (fit failed).", dx))
      next
    }
    out <- data.frame(b1 = gridXY[, 1], b2 = gridXY[, 2], p = pg, dx = dx)
    fp <- file.path(paste0(OUTPUTS_DIR, "/base_prob_grids"), paste0("probgrid_", make_fname_safe(dx), ".csv"))
    write_csv(out, fp)
    dx_surface_files <- c(dx_surface_files, fp)
  }
  cat(sprintf("[base-field] wrote %d grid(s) to %s\n", length(dx_surface_files), paste0(OUTPUTS_DIR, "/base_prob_grids/")))
  
  # --- Item x Diagnosis Interactions (Parallelized OOF GAM) ---
  
  vars_item_interact <- unique(varmap)
  n_probe_count <- if (exists("N_TOP_PER_DX")) N_TOP_PER_DX else 50
  vars_probe <- head(vars_item_interact, min(n_probe_count, length(vars_item_interact)))
  
  # Pre-calculate scaled matrix for speed
  if (!exists("Z0_std")) Z0_std <- scale(Xenc, center = TRUE, scale = TRUE)
  
  task_grid <- expand.grid(
    dx = names(DxW_A), 
    var = vars_probe, 
    stringsAsFactors = FALSE
  )
  
  # Filter rare diagnoses for stability
  dx_counts <- colSums(DxW_A > 0, na.rm = TRUE)
  valid_dx <- names(dx_counts)[dx_counts >= 6]
  task_grid <- task_grid[task_grid$dx %in% valid_dx, , drop = FALSE]
  
  cat(sprintf("[Interactions] Parallelizing %d tasks (%d items x %d diagnoses)...\n", 
              nrow(task_grid), length(unique(task_grid$var)), length(unique(task_grid$dx))))
  
  # Helper Helpers
  choose_K_dx <- function(y, K_target = 5L, min_per_class = 6L) {
    y <- as.integer(y > 0)
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    max(2L, min(K_target, floor(n1 / min_per_class), floor(n0 / min_per_class)))
  }
  
  oof_R2_two_gams <- function(v, Base, dx, K_target = 5, k_gam = 10,
                              seed = SEED_GLOBAL, seed_key = NULL) {
    n <- length(v)
    if (n != nrow(Base) || n != length(dx)) {
      return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
    }
    
    K <- choose_K_dx(dx, K_target = K_target, min_per_class = 6L)
    if (K < 2) return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
    
    seed_eff <- if (is.null(seed_key)) {
      as.integer(seed)
    } else {
      .seed_from_key(seed, paste0("oof_gam|", seed_key))
    }
    fold_id <- .with_seed(seed_eff, sample(rep(1:K, length.out = n)))
    lev_all <- levels(factor(dx))
    
    y_add <- rep(NA_real_, n)
    y_int <- rep(NA_real_, n)
    
    # Note: Force single-threaded GAM fitting here since we are inside parallel loop
    ctrl <- mgcv::gam.control(maxit = 100, nthreads = 1)
    
    for (k in sort(unique(fold_id))) {
      tr <- which(fold_id != k)
      te <- which(fold_id == k)
      if (length(tr) < 10 || length(te) == 0) next 
      
      dtr <- data.frame(v = v[tr], b1 = Base[tr, 1], b2 = Base[tr, 2], dx = factor(dx[tr], levels = lev_all))
      dte <- data.frame(b1 = Base[te, 1], b2 = Base[te, 2], dx = factor(dx[te], levels = lev_all))
      
      # 1. Additive Model
      f_add <- try(mgcv::gam(v ~ s(b1, b2, k = k_gam, bs = "tp", m = 2),
                             data = dtr, method = "REML", gamma = 1.4, control = ctrl), silent = TRUE)
      
      # 2. Interaction Model
      f_int <- try(mgcv::gam(v ~ s(b1, b2, k = k_gam, bs = "tp", m = 2) + 
                               dx +
                               s(b1, b2, by = dx, k = k_gam, bs = "tp", m = 2),
                             data = dtr, method = "REML", select = TRUE, gamma = 1.4, control = ctrl), silent = TRUE)
      
      mu <- mean(dtr$v, na.rm = TRUE)
      
      if (!inherits(f_add, "try-error")) {
        pa <- try(predict(f_add, newdata = dte, type = "response"), silent = TRUE)
        if (!inherits(pa, "try-error")) y_add[te] <- as.numeric(pa)
      }
      
      if (!inherits(f_int, "try-error")) {
        pi <- try(predict(f_int, newdata = dte, type = "response"), silent = TRUE)
        if (!inherits(pi, "try-error")) y_int[te] <- as.numeric(pi)
      }
      
      if (any(!is.finite(y_add[te]))) y_add[te][!is.finite(y_add[te])] <- mu
      if (any(!is.finite(y_int[te]))) y_int[te][!is.finite(y_int[te])] <- mu
    }
    
    ve <- stats::var(v, na.rm = TRUE)
    if (!is.finite(ve) || ve <= 0) return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
    
    R2_add <- max(0, 1 - mean((v - y_add)^2, na.rm = TRUE) / ve)
    R2_int <- max(0, 1 - mean((v - y_int)^2, na.rm = TRUE) / ve)
    
    d_sq <- (v - y_add)^2 - (v - y_int)^2
    d_sq <- d_sq[is.finite(d_sq)]
    
    p_like <- if (length(d_sq) < 10L || all(abs(d_sq) < .Machine$double.eps)) NA_real_ else 
      tryCatch(stats::wilcox.test(d_sq, mu = 0, alternative = "greater", exact = FALSE)$p.value,
               error = function(e) NA_real_)
    
    c(R2_add = R2_add, R2_int = R2_int, p_like = p_like, dR2 = R2_int - R2_add)
  }
  
  # Execute Parallel Loop
  n_workers_interact <- tryCatch(as.integer(future::nbrOfWorkers()), error = function(e) 1L)
  worker_parent <- new.env(parent = globalenv())
  score_item_base_core <- score_item_base
  environment(score_item_base_core) <- worker_parent

  choose_K_dx_core <- choose_K_dx
  environment(choose_K_dx_core) <- worker_parent

  hash32_core <- .hash32
  environment(hash32_core) <- worker_parent

  seed_from_key_core <- .seed_from_key
  environment(seed_from_key_core) <- list2env(
    list(.hash32 = hash32_core),
    parent = worker_parent
  )

  with_seed_core <- .with_seed
  environment(with_seed_core) <- worker_parent

  oof_R2_two_gams_core <- oof_R2_two_gams
  environment(oof_R2_two_gams_core) <- list2env(
    list(
      choose_K_dx = choose_K_dx_core,
      .seed_from_key = seed_from_key_core,
      .with_seed = with_seed_core,
      SEED_GLOBAL = SEED_GLOBAL
    ),
    parent = worker_parent
  )

  interact_worker_core <- interact_worker
  environment(interact_worker_core) <- list2env(
    list(
      score_item_base = score_item_base_core,
      oof_R2_two_gams = oof_R2_two_gams_core,
      SEED_GLOBAL = SEED_GLOBAL
    ),
    parent = worker_parent
  )
  interact_worker_env <- list2env(
    list(
      task_grid = task_grid,
      Z0_std = Z0_std,
      Base = Base,
      DxW_A = DxW_A,
      varmap = varmap,
      interact_worker_core = interact_worker_core
    ),
    parent = worker_parent
  )
  interact_task <- function(i) {
    interact_worker_core(i, task_grid, Z0_std, Base, DxW_A, varmap)
  }
  environment(interact_task) <- interact_worker_env

  if (n_workers_interact <= 1L) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(task_grid), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
    rows <- lapply(seq_len(nrow(task_grid)), function(i) {
      out <- interact_task(i)
      utils::setTxtProgressBar(pb, i)
      out
    })
  } else {
    rows <- FUTURE_LAPPLY(
      seq_len(nrow(task_grid)),
      interact_task,
      future.packages = c("mgcv", "stats"),
      future.seed = TRUE,
      future.scheduling = 1
    )
  }
  
  # Aggregate Results
  rows <- Filter(Negate(is.null), rows)
  
  if (length(rows) > 0) {
    int_tbl <- dplyr::bind_rows(rows) %>% dplyr::filter(is.finite(dR2), is.finite(p_like))
    
    if (nrow(int_tbl) > 0) {
      int_tbl <- int_tbl %>%
        dplyr::group_by(dx) %>%
        dplyr::mutate(q_like = p.adjust(p_like, method = "BH")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sig = q_like < SIG_Q) %>%
        dplyr::arrange(dplyr::desc(dR2))
      
      write_csv(int_tbl, "item_dx_interactions.csv")
      cat(sprintf("[interactions] wrote: %s (rows=%d)\n", "item_dx_interactions.csv", nrow(int_tbl)))
    } else {
      warning("[interactions] Table empty after filtering.")
    }
  } else {
    warning("[interactions] No valid interactions computed.")
  }
} else {
  message("[dx] Section 8 (predictive diagnostics & interactions) disabled: DX set as optional (setup.R)")
}

# ==============================================================================
# 9. Outputs and Session
# ==============================================================================

if (DX_AVAILABLE && !is.null(DX_wide) && ncol(DX_wide) > 0) {

  # Cluster Summaries
  kF <- length(unique(clF))
  cat(sprintf("[clusters] using %d clusters.\n", kF))
  
  clusters <- sort(unique(clF))
  cl_fac <- factor(clF, levels = clusters)
  
  mean_by_cluster <- function(vec, clf) {
    out <- tapply(vec, clf, function(z) mean(z, na.rm = TRUE))
    as.numeric(out)
  }
  
  n_by_cluster <- as.integer(tabulate(cl_fac))
  base_means <- cbind(b1_mean = mean_by_cluster(Base[, 1], cl_fac), b2_mean = mean_by_cluster(Base[, 2], cl_fac))
  
  f_means <- sapply(seq_len(ncol(Bprime)), function(j) mean_by_cluster(Bprime[, j], cl_fac))
  if (!is.null(dim(f_means))) {
    colnames(f_means) <- paste0("f", seq_len(ncol(Bprime)), "_mean")
  } else {
    f_means <- cbind(`f1_mean` = as.numeric(f_means))
  }
  
  sum_tbl <- data.frame(cluster = clusters, n = n_by_cluster, base_means, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(f_means)) {
    f_means <- as.matrix(f_means)
    storage.mode(f_means) <- "double"
    sum_tbl <- cbind(sum_tbl, f_means, stringsAsFactors = FALSE)
  }
  write_csv(sum_tbl, "cluster_summary.csv")
  
  # Cluster Enrichment Analysis
  DxW2 <- DxW_A
  DxW2[is.na(DxW2)] <- 0L
  prev2 <- colSums(DxW2 > 0L, na.rm = TRUE)
  keepc <- prev2 >= max(DX_CASES_MIN, ceiling(0.01 * nrow(DxW2)))
  DxW2 <- DxW2[, keepc, drop = FALSE]
  
  if (ncol(DxW2) > 0) {
    tab <- sapply(colnames(DxW2), function(dn) {
      tapply(DxW2[[dn]] > 0L, factor(clF, levels = clusters), sum, na.rm = TRUE)
    })
    tab <- as.matrix(tab)
    storage.mode(tab) <- "double"
    rownames(tab) <- as.character(clusters)
    
    suppressWarnings({
      chi <- try(chisq.test(tab), silent = TRUE)
    })
    
    if (inherits(chi, "try-error")) {
      Eexp_used <- outer(rowSums(tab), colSums(tab), function(r, c) r * c / max(sum(tab), 1))
      Z_enrich <- (tab - Eexp_used) / sqrt(pmax(Eexp_used, 1e-9))
      pmat <- 2 * pnorm(-abs(Z_enrich))
    } else {
      Z_enrich <- chi$stdres
      Eexp_used <- chi$expected
      pmat <- 2 * pnorm(-abs(Z_enrich))
    }
    
    # BH corrections
    q_by_dx <- apply(pmat, 2, function(p) p.adjust(p, method = "BH"))
    q_global <- p.adjust(as.vector(pmat), method = "BH")
    q_global_mat <- matrix(
      q_global,
      nrow = nrow(pmat), ncol = ncol(pmat),
      byrow = FALSE, dimnames = dimnames(pmat)
    )
    
    # Prevalences
    ncl <- as.integer(tabulate(factor(clF, levels = clusters)))
    names(ncl) <- as.character(clusters)
    prev_in_cluster <- sweep(tab, 1, ncl, "/")
    prev_overall <- colSums(DxW2 > 0L) / nrow(DxW2)
    
    # Export wide matrices
    write_csv(data.frame(cluster = rownames(Z_enrich), as.data.frame(Z_enrich), check.names = FALSE), "dx_cluster_enrichment_Z.csv")
    write_csv(data.frame(cluster = rownames(pmat), as.data.frame(pmat), check.names = FALSE), "dx_cluster_enrichment_p.csv")
    write_csv(data.frame(cluster = rownames(q_by_dx), as.data.frame(q_by_dx), check.names = FALSE), "dx_cluster_enrichment_q_by_dx.csv")
    write_csv(data.frame(cluster = rownames(q_global_mat), as.data.frame(q_global_mat), check.names = FALSE), "dx_cluster_enrichment_q_global.csv")
    write_csv(data.frame(cluster = rownames(prev_in_cluster), as.data.frame(prev_in_cluster), check.names = FALSE), "dx_cluster_prev_in_cluster.csv")
    
    # Export Long table
    grid <- expand.grid(
      cluster = rownames(tab),
      dx = colnames(tab),
      stringsAsFactors = FALSE
    )
    idx_r <- match(grid$cluster, rownames(tab))
    idx_c <- match(grid$dx, colnames(tab))
    
    grid$n_cluster <- ncl[grid$cluster]
    grid$cases_cluster <- tab[cbind(idx_r, idx_c)]
    grid$prev_cluster <- prev_in_cluster[cbind(idx_r, idx_c)]
    grid$prev_overall <- prev_overall[grid$dx]
    grid$Z <- Z_enrich[cbind(idx_r, idx_c)]
    grid$p <- pmat[cbind(idx_r, idx_c)]
    grid$q_by_dx <- q_by_dx[cbind(idx_r, idx_c)]
    grid$q_global <- q_global_mat[cbind(idx_r, idx_c)]
    grid$expected_cases <- Eexp_used[cbind(idx_r, idx_c)]
    grid$enriched <- (grid$Z > 0) & (grid$q_by_dx < SIG_Q)
    
    grid <- grid[order(-grid$Z, grid$q_by_dx, grid$dx, as.integer(grid$cluster)), ]
    write_csv(grid, "dx_cluster_enrichment_long.csv")
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Cap Z for display
      Zcap <- pmax(pmin(Z_enrich, 5), -5)
      z_long <- data.frame(
        cluster = rep(rownames(Zcap), times = ncol(Zcap)),
        dx = rep(colnames(Zcap), each = nrow(Zcap)),
        Z = as.vector(Zcap),
        stringsAsFactors = FALSE
      )
      gp <- ggplot2::ggplot(z_long, ggplot2::aes(x = dx, y = factor(cluster))) +
        ggplot2::geom_tile(ggplot2::aes(fill = Z)) +
        ggplot2::scale_fill_gradient2() +
        ggplot2::labs(
          x = "Diagnosis", y = "Cluster", fill = "Z (std. resid.)",
          title = "Diagnosis enrichment by cluster (std. residuals)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      print(gp)
      save_plot_gg("FIG_dx_cluster_enrichment_heatmap.png",
                   gp,
                   width = 10, height = 6, dpi = 150
      )
    }
  } else {
    warning("[clusters] No Dx columns after filtering for enrichment; skipped.")
  }
} else {
  message("[dx] Enrichment by diagnosis disabled: DX set as optional.")
}

# ------------------------------------------------------------------------------
# Encoding & geometry
# ------------------------------------------------------------------------------

as_base2_numeric <- function(Base_in) {
  A <- as.data.frame(Base_in[, 1:2, drop = FALSE])
  for (j in 1:2) {
    v <- suppressWarnings(as.numeric(A[[j]]))
    m <- mean(v, na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    A[[j]] <- ifelse(is.finite(v), v, m)
  }
  B <- as.matrix(A)
  storage.mode(B) <- "double"
  colnames(B) <- c("b1", "b2")
  rownames(B) <- rownames(Base_in)
  B
}

# whiten base to unit disk (pure linear algebra)
standardise_to_circle <- function(Base2, cover = 0.95) {
  X <- as.matrix(Base2[, 1:2, drop = FALSE]); storage.mode(X) <- "double"
  mu <- colMeans(X); S <- stats::cov(X)
  if (!all(is.finite(S)) || det(S) <= 0) S <- S + diag(1e-8, 2)
  
  eig <- eigen(S, symmetric = TRUE)
  V   <- eig$vectors
  L   <- pmax(eig$values, 1e-12)
  S_half <- V %*% diag(sqrt(L), 2) %*% t(V)
  S_hi   <- V %*% diag(1/sqrt(L), 2) %*% t(V)
  
  Xc <- sweep(X, 2, mu, "-")
  U0 <- Xc %*% t(S_hi)
  r  <- sqrt(rowSums(U0^2))
  s  <- as.numeric(stats::quantile(r, probs = cover, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) s <- max(r, na.rm = TRUE)
  
  fwd <- function(xb1, xb2) {
    xb1 <- as.numeric(xb1); xb2 <- as.numeric(xb2)
    Xnew <- cbind(xb1, xb2)
    Uq   <- (sweep(Xnew, 2, mu, "-") %*% t(S_hi)) / s
    colnames(Uq) <- c("u1", "u2"); Uq
  }
  inv <- function(u1, u2) {
    Unew <- cbind(as.numeric(u1), as.numeric(u2)) * s
    Xq   <- Unew %*% t(S_half) + matrix(mu, nrow(Unew), 2, byrow = TRUE)
    colnames(Xq) <- c("b1", "b2"); Xq
  }
  
  list(mu = mu, S_half = S_half, S_half_inv = S_hi, s = s, fwd = fwd, inv = inv)
}

make_unitdisk_square <- function(nu) {
  u1 <- seq(-1, 1, length.out = nu)
  u2 <- seq(-1, 1, length.out = nu)
  g <- expand.grid(u1 = u1, u2 = u2)
  mask <- with(g, sqrt(u1^2 + u2^2) <= 1 + 1e-9)
  list(u1 = u1, u2 = u2, grid = g, mask_vec = mask)
}

inside_hull <- function(px, py, poly, eps = 1e-12) {
  n <- nrow(poly)
  j <- n
  inside <- rep(FALSE, length(px))
  on_edge <- rep(FALSE, length(px))
  for (i in seq_len(n)) {
    xi <- poly$b1[i]; yi <- poly$b2[i]
    xj <- poly$b1[j]; yj <- poly$b2[j]
    dx <- xj - xi; dy <- yj - yi
    seg_len2 <- dx * dx + dy * dy + eps
    t <- pmin(1, pmax(0, ((px - xi) * dx + (py - yi) * dy) / seg_len2))
    projx <- xi + t * dx; projy <- yi + t * dy
    on_edge <- on_edge | ((px - projx)^2 + (py - projy)^2 <= (1e-9)^2)
    cross <- ((yi > py) != (yj > py)) & (px < (xj - xi) * (py - yi) / (yj - yi + eps) + xi)
    inside <- xor(inside, cross)
    j <- i
  }
  inside | on_edge
}

build_geometry <- function(Base_A, GRID_N_B = 600, GRID_N_U = 600) {
  B1 <- Base_A[, 1]; B2 <- Base_A[, 2]
  df_base <- data.frame(b1 = as.numeric(B1), b2 = as.numeric(B2))
  rownames(df_base) <- rownames(Base_A)
  
  Xstd <- standardise_to_circle(Base_A, cover = 0.95)
  U <- as.data.frame(Xstd$fwd(B1, B2))
  names(U) <- c("u1", "u2")
  rownames(U) <- rownames(Base_A)
  
  UD <- make_unitdisk_square(GRID_N_B)
  gridX_sq <- as.data.frame(Xstd$inv(UD$grid$u1, UD$grid$u2))
  mask_sq <- UD$mask_vec
  
  hidx <- grDevices::chull(df_base$b1, df_base$b2)
  hpoly <- df_base[hidx, , drop = FALSE]
  
  qx <- range(df_base$b1); qy <- range(df_base$b2)
  gridB_full <- expand.grid(
    b1 = seq(qx[1], qx[2], length.out = GRID_N_U),
    b2 = seq(qy[1], qy[2], length.out = GRID_N_U)
  )
  mask_hull <- inside_hull(gridB_full$b1, gridB_full$b2, hpoly)
  
  U$rin  <- sqrt(U$u1^2 + U$u2^2) <= 1 + 1e-9
  U$insq <- abs(U$u1) <= 1 + 1e-9 & abs(U$u2) <= 1 + 1e-9
  
  list(
    Xstd = Xstd, U = U, UD = UD,
    gridX_sq = gridX_sq, mask_sq = mask_sq,
    df_base = df_base, gridB_full = gridB_full, mask_hull = mask_hull
  )
}

# Basic objects
varmap <- attr(Z, "varmap")
Z_A    <- Z

# Geometry & grids
geom <- build_geometry(Base_A)

# Export PC score tables
pc_scores_base <- data.frame(
  participant_id = ids_base,
  b1 = Base[, 1], b2 = Base[, 2],
  cluster = as.integer(clF),
  target = as.integer(y_use),
  stringsAsFactors = FALSE
)
write_csv(pc_scores_base, "pc_scores_base_b1b2.csv")

pc_scores_residual <- data.frame(
  participant_id = ids_base,
  Bprime,
  cluster = as.integer(clF),
  stringsAsFactors = FALSE
)
write_csv(pc_scores_residual, "pc_scores_residual_Bprime.csv")

saveRDS(list(
  participant_id = ids_base, 
  XR = E_scaled, 
  Fprime = Fprime
),
file = "Fprime_matrix.rds"
)

# Export Weights and Maps
w_tbl <- data.frame(
  var = names(w_full),
  weight = as.numeric(w_full),
  selected = names(w_full) %in% survivors,
  stringsAsFactors = FALSE
)
write_csv(w_tbl, "gower_weights_id_guided.csv")

enc_map <- data.frame(
  mm_col = colnames(Xenc),
  source_var = as.character(varmap),
  weight_share = as.numeric(w_enc[colnames(Xenc)]),
  stringsAsFactors = FALSE
)
write_csv(enc_map, "encoding_map_and_weight_share.csv")

saveRDS(list(
  seeds = list(
    rng_kind = RNG_KIND,
    seed_global = SEED_GLOBAL,
    seed_pred = SEED_PRED,
    seed_jitter = SEED_JITTER,
    seed_boot = SEED_BOOT,
    bundle_seed = BUNDLE_SEED
  ),
  folds = list(
    target = fold_id,
    fibre = folds_f
  ),
  data_quality = list(
    constant_profile_n = if (is.null(degenerate_constant_profiles)) 0L else nrow(degenerate_constant_profiles),
    constant_profile_ids = if (is.null(degenerate_constant_profiles)) character(0) else degenerate_constant_profiles$participant_id,
    constant_profile_values = if (is.null(degenerate_constant_profiles)) character(0) else degenerate_constant_profiles$constant_value
  ),
  decomposition = list(
    method = BASE_DECOMP_METHOD,
    m_star = m_star,
    k_eff = k_eff
  ),
  gower = list(
    mode = WEIGHTING_MODE,
    reps = reps,
    core_idx_rep = core_idx_rep,
    idx_used = wopt$idx_used,
    history = wopt$history,
    final_id = wopt$final_ID,
    survivors = survivors,
    multi_run_enabled = isTRUE(GOWER_MULTI_ENABLE && GOWER_MULTI_RUNS > 1L),
    multi_run_runs = if (!is.null(wopt_list)) length(wopt_list) else 1L
  )
), file = "qc_manifest.rds")

# Session Info
zz <- file("sessionInfo.txt", open = "wt")
sink(zz)
print(sessionInfo())
cat("\n\nConfig snapshot:\n")
print(str(list(
  outputs_dir = OUTPUTS_DIR,
  weighting_mode = WEIGHTING_MODE,
  base_decomp_method = BASE_DECOMP_METHOD,
  rng = list(
    kind = RNG_KIND,
    seed_global = SEED_GLOBAL,
    seed_pred = SEED_PRED,
    seed_jitter = SEED_JITTER,
    seed_boot = SEED_BOOT
  ),
  compute = list(
    bam_threads = BAM_THREADS,
    omp_threads = OMP_THREADS,
    blas_threads = BLAS_THREADS,
    ncores_par = NCORES_PAR
  ),
  palette = list(
    engine = PALETTE_ENGINE,
    name = PALETTE_NAME,
    name_div = PALETTE_NAME_DIV,
    direction = PALETTE_DIRECTION
  ),
  toggles = list(
    do_plots = DO_PLOTS,
    do_diagnostics = DO_DIAGNOSTICS,
    do_surface = DO_SURFACE,
    do_sweep = DO_SWEEP
  ),
  pc_score_map = list(
    knn_k = KNN_K,
    knn_variant = KNN_VARIANT,
    local_scale = LOCAL_SCALE,
    multiplicity_weight = MULT_WEIGHT,
    m_star_fixed = M_STAR_FIXED,
    id_k_range = K_ID_LO_HI,
    pc_max = PC_MAX,
    m_default = M_DEFAULT
  )
), max.level = 1))
sink()
close(zz)

cat("[export] wrote score tables, enrichment tables, weights, and session info.\n")

# Check alignment
cat("[Orientation] Base head:\n")
print(head(Base[, 1:2]))
cat("[Orientation] Base_A head:\n")
print(head(Base_A))

# ------------------------------------------------------------------------------
# Base-space Plots
# ------------------------------------------------------------------------------

density_circle_outline <- function(radius = 1, n = 361) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(u1 = radius * cos(th), u2 = radius * sin(th))
}

density_plots <- function(geom) {
  stopifnot(all(c("U") %in% names(geom)))
  U_df <- geom$U
  U_df <- U_df[is.finite(U_df$u1) & is.finite(U_df$u2), , drop = FALSE]
  if (!nrow(U_df)) return(NULL)
  
  r <- sqrt(U_df$u1^2 + U_df$u2^2)
  lim <- max(1.05, as.numeric(stats::quantile(r, probs = 0.995, na.rm = TRUE)))
  lim <- 1.05 * lim
  
  ggplot() +
    stat_density_2d(
      data = U_df,
      aes(u1, u2, colour = after_stat(level)),
      bins = 6, linewidth = 0.45, alpha = 0.9
    ) +
    { scico::scale_colour_scico(palette = PALETTE_NAME, direction = 1) } +
    geom_point(
      data = U_df, aes(u1, u2),
      shape = 16, size = 0.9, colour = scales::alpha("black", 0.35)
    ) +
    geom_path(
      data = density_circle_outline(),
      aes(u1, u2),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.35
    ) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = "u1 (whitened/scaled b1,b2)", y = "u2") +
    theme_pub(12)
}

density_plots_base <- function(Base_A, pad_frac = 0.06) {
  dfB <- data.frame(
    b1 = as.numeric(Base_A[, 1]),
    b2 = as.numeric(Base_A[, 2]),
    stringsAsFactors = FALSE
  )
  dfB <- dfB[is.finite(dfB$b1) & is.finite(dfB$b2), , drop = FALSE]
  if (!nrow(dfB)) return(NULL)
  
  rx <- range(dfB$b1, na.rm = TRUE)
  ry <- range(dfB$b2, na.rm = TRUE)
  wx <- diff(rx)
  wy <- diff(ry)
  if (!is.finite(wx) || wx <= 0) wx <- 1
  if (!is.finite(wy) || wy <= 0) wy <- 1
  
  xlim <- c(rx[1] - pad_frac * wx, rx[2] + pad_frac * wx)
  ylim <- c(ry[1] - pad_frac * wy, ry[2] + pad_frac * wy)
  
  ggplot() +
    stat_density_2d(
      data = dfB,
      aes(b1, b2, colour = after_stat(level)),
      bins = 6, linewidth = 0.45, alpha = 0.9
    ) +
    { scico::scale_colour_scico(palette = PALETTE_NAME, direction = 1) } +
    geom_point(
      data = dfB, aes(b1, b2),
      shape = 16, size = 0.9, colour = scales::alpha("black", 0.35)
    ) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = "b1", y = "b2") +
    theme_pub(12)
}

direction_wheel_plot <- function(geom) {
  U_disk <- subset(geom$U, rin)
  U_disk <- U_disk[is.finite(U_disk$u1) & is.finite(U_disk$u2), , drop = FALSE]
  if (!nrow(U_disk)) return(NULL)
  
  nu <- 500
  pad <- 0.75
  gx <- seq(-1 - pad, 1 + pad, length.out = nu)
  gy <- seq(-1 - pad, 1 + pad, length.out = nu)
  
  kd <- with(U_disk, MASS::kde2d(u1, u2, n = nu,
                                 lims = c(-1 - pad, 1 + pad, -1 - pad, 1 + pad)))
  D <- kd$z
  D <- log1p(D / max(D, na.rm = TRUE))
  D <- D / quantile(D, 0.99, na.rm = TRUE)
  D[D > 1] <- 1; D[D < 0] <- 0
  ALPHA <- D^0.70
  
  G <- expand.grid(u1 = gx, u2 = gy)
  theta <- atan2(G$u2, G$u1)
  r <- sqrt(G$u1^2 + G$u2^2)
  
  H0 <- 170; L0 <- 60; Cmax <- 110; betaC <- 0.90
  H <- (H0 + theta * 180 / pi) %% 360
  C <- pmin(Cmax * (pmin(r, 1)^betaC), Cmax)
  L <- pmax(0, pmin(100, L0 - 6 * (pmin(r, 1)^1.1)))
  
  G$fill <- grDevices::hcl(H, C, L)
  G$alpha <- as.vector(ALPHA)
  
  feather_disk <- function(r, radius = 1, w = 0.08) {
    a <- rep(0, length(r))
    a[r <= (radius - w)] <- 1
    edge <- r > (radius - w) & r <= radius
    a[edge] <- (cos(((r[edge] - (radius - w)) / w) * pi / 2))^2
    a
  }
  G$alpha <- G$alpha * feather_disk(r)
  
  anchor <- data.frame(
    x = c(0.85, -0.85, 0.02, 0.02),
    y = c(0.02, 0.02, 0.85, -0.85),
    lab = c("+b1", "-b1", "+b2", "-b2"),
    col = grDevices::hcl((H0 + c(0, 180, 90, -90)) %% 360, Cmax, L0)
  )
  
  ggplot() +
    geom_raster(data = G, aes(u1, u2, fill = I(fill), alpha = alpha), interpolate = TRUE) +
    scale_alpha(range = c(0, 1), guide = "none") +
    geom_point(data = U_disk, aes(u1, u2),
               shape = 16, size = 0.8, colour = scales::alpha("black", 0.32)) +
    geom_path(
      data = density_circle_outline(),
      aes(u1, u2),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.35
    ) +
    coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5),
                expand = FALSE, clip = "on") +
    geom_point(data = anchor, aes(x, y), shape = 15, size = 3, colour = anchor$col) +
    geom_text(data = anchor, aes(x, y, label = lab), nudge_x = 0.07, size = 3.2) +
    labs(x = "u1 (whitened/scaled b1,b2)", y = "u2") +
    theme_pub(12)
}

direction_wheel_plot_base <- function(Base_A, cover = 1, pad_frac = 0.06) {
  dfB <- data.frame(
    b1 = as.numeric(Base_A[, 1]),
    b2 = as.numeric(Base_A[, 2]),
    stringsAsFactors = FALSE
  )
  dfB <- dfB[is.finite(dfB$b1) & is.finite(dfB$b2), , drop = FALSE]
  if (!nrow(dfB)) return(NULL)
  
  Xstd <- standardise_to_circle(Base_A, cover = cover)
  mu <- Xstd$mu
  S_hi <- Xstd$S_half_inv
  s <- Xstd$s
  
  nu <- 500
  
  rx <- range(dfB$b1, na.rm = TRUE)
  ry <- range(dfB$b2, na.rm = TRUE)
  wx <- diff(rx)
  wy <- diff(ry)
  if (!is.finite(wx) || wx <= 0) wx <- 1
  if (!is.finite(wy) || wy <= 0) wy <- 1
  
  lims <- c(
    rx[1] - pad_frac * wx, rx[2] + pad_frac * wx,
    ry[1] - pad_frac * wy, ry[2] + pad_frac * wy
  )
  
  kd <- with(dfB, MASS::kde2d(b1, b2, n = nu, lims = lims))
  D <- kd$z
  D <- log1p(D / max(D, na.rm = TRUE))
  D <- D / stats::quantile(D, 0.99, na.rm = TRUE)
  D[D > 1] <- 1
  D[D < 0] <- 0
  ALPHA <- D^0.70
  
  gx <- kd$x
  gy <- kd$y
  G <- expand.grid(b1 = gx, b2 = gy)
  
  Uq <- (sweep(as.matrix(G), 2, mu, "-") %*% t(S_hi)) / s
  theta <- atan2(Uq[, 2], Uq[, 1])
  r <- sqrt(Uq[, 1]^2 + Uq[, 2]^2)
  
  H0 <- 170
  L0 <- 60
  Cmax <- 110
  betaC <- 0.90
  H <- (H0 + theta * 180 / pi) %% 360
  C <- pmin(Cmax * (pmin(r, 1)^betaC), Cmax)
  L <- pmax(0, pmin(100, L0 - 6 * (pmin(r, 1)^1.1)))
  
  G$fill <- grDevices::hcl(H, C, L)
  G$alpha <- as.vector(ALPHA)
  
  ggplot() +
    geom_raster(
      data = G,
      aes(b1, b2, fill = I(fill), alpha = alpha),
      interpolate = TRUE
    ) +
    scale_alpha(range = c(0, 1), guide = "none") +
    geom_point(
      data = dfB,
      aes(b1, b2),
      shape = 16, size = 0.8,
      colour = scales::alpha("black", 0.32)
    ) +
    coord_equal(
      xlim = lims[1:2], ylim = lims[3:4],
      expand = FALSE, clip = "on"
    ) +
    labs(x = "b1", y = "b2") +
    theme_pub(12)
}

# Participant density plot
p_dens <- density_plots(geom)
save_plot_gg("FIG_dens_unitsquare_scatter", p_dens, width = 8.0, height = 7.0)

p_dens_base <- density_plots_base(Base_A)
if (!is.null(p_dens_base)) {
  save_plot_gg("FIG_dens_base_scatter", p_dens_base, width = 8.0, height = 7.0)
}

# Direction plot
p_dir <- direction_wheel_plot(geom)
if (!is.null(p_dir)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth", p_dir, width = 8.0, height = 7.0)
}

p_dir_base <- direction_wheel_plot_base(Base_A)
if (!is.null(p_dir_base)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth_BASE", p_dir_base, width = 8.0, height = 7.0)
}

# ------------------------------------------------------------------------------
# Biplots
# ------------------------------------------------------------------------------

score_item_1d <- function(nm, Z, varmap) {
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  v <- if (length(idx) == 1L) {
    as.numeric(Z[, idx])
  } else {
    sc <- try(suppressWarnings(prcomp(Z[, idx, drop = FALSE], rank. = 1)$x[, 1]), silent = TRUE)
    if (inherits(sc, "try-error")) rep(NA_real_, nrow(Z)) else as.numeric(sc)
  }
  as.numeric(scale(v))
}

pearson_r <- function(a, b) {
  if (length(a) != length(b)) return(NA_real_)
  suppressWarnings(cor(as.numeric(a), as.numeric(b), use = "complete.obs", method = "pearson"))
}

build_biplot_data <- function(Z_A, varmap, Base_A, U) {
  B1 <- Base_A[, 1]; B2 <- Base_A[, 2]
  items <- unique(varmap)
  Rtab <- dplyr::bind_rows(lapply(items, function(nm) {
    v <- score_item_1d(nm, Z_A, varmap)
    if (!any(is.finite(v))) return(NULL)
    data.frame(
      item = nm,
      r_b1 = pearson_r(v, B1), r_b2 = pearson_r(v, B2),
      r_u1 = pearson_r(v, U$u1), r_u2 = pearson_r(v, U$u2),
      stringsAsFactors = FALSE
    )
  }))
  if (!nrow(Rtab)) stop("[biplot] No item correlations computed.")
  Rtab |>
    dplyr::mutate(
      across(c(r_b1, r_b2, r_u1, r_u2), ~ suppressWarnings(as.numeric(.))),
      r_b1 = ifelse(is.finite(r_b1), r_b1, 0),
      r_b2 = ifelse(is.finite(r_b2), r_b2, 0),
      r_u1 = ifelse(is.finite(r_u1), r_u1, 0),
      r_u2 = ifelse(is.finite(r_u2), r_u2, 0),
      mag_r_base = sqrt(r_b1^2 + r_b2^2),
      mag_r_disk = sqrt(r_u1^2 + r_u2^2)
    )
}

plot_biplots <- function(Rtab, Base_A, U) {
  suppressPackageStartupMessages(requireNamespace("ggrepel", quietly = TRUE))
  
  B1n <- suppressWarnings(as.numeric(Base_A[, 1]))
  B2n <- suppressWarnings(as.numeric(Base_A[, 2]))
  B1n[!is.finite(B1n)] <- mean(B1n, na.rm = TRUE)
  B2n[!is.finite(B2n)] <- mean(B2n, na.rm = TRUE)
  
  Hidx <- grDevices::chull(B1n, B2n)
  H <- data.frame(b1 = B1n[Hidx], b2 = B2n[Hidx])
  
  cx <- mean(B1n); cy <- mean(B2n)
  Rscale <- 0.80 * min(diff(range(B1n)), diff(range(B2n)))
  
  # Helper to assign octants based on angle (-pi to pi)
  # Octant 1: -22.5 to +22.5 (East)
  # Octant 2: +22.5 to +67.5 (North-East) ... etc.
  get_octant <- function(x, y) {
    theta_deg <- atan2(y, x) * 180 / pi
    # Rotate so -22.5 becomes 0 for easier flooring
    # Result: 0=East, 1=NE, 2=N, 3=NW, 4=W, 5=SW, 6=S, 7=SE
    idx <- floor((theta_deg + 22.5 + 360) %% 360 / 45) + 1
    paste0("Octant_", idx)
  }
  
  # Take top 4 items from each of the 8 octants (32 items total)
  S_base <- Rtab |>
    dplyr::mutate(
      octant = get_octant(r_b1, r_b2)
    ) |>
    dplyr::group_by(octant) |>
    dplyr::arrange(dplyr::desc(mag_r_base)) |>
    dplyr::slice_head(n = 4) |>  # Top 4 per octant
    dplyr::ungroup() |>
    dplyr::mutate(
      x0 = cx, 
      y0 = cy,
      x1 = cx + Rscale * r_b1,
      y1 = cy + Rscale * r_b2
    )
  
  p_base <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = H, ggplot2::aes(b1, b2),
                          fill = NA, colour = "black", linewidth = 0.4) +
    ggplot2::geom_segment(
      data = S_base,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_label_repel(
      data = S_base,
      ggplot2::aes(x = x1, y = y1, label = item),
      size = 3.1, max.overlaps = Inf, label.size = 0,
      label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12)
  
  draw_disk_outline <- function() {
    th <- seq(0, 2 * pi, length.out = 361)
    data.frame(x = cos(th), y = sin(th))
  }
  
  Rdisk <- 0.85
  
  S_disk <- Rtab |>
    dplyr::mutate(
      octant = get_octant(r_u1, r_u2)
    ) |>
    dplyr::group_by(octant) |>
    dplyr::arrange(dplyr::desc(mag_r_disk)) |>
    dplyr::slice_head(n = 4) |>  # Top 4 per octant in U-space
    dplyr::ungroup() |>
    dplyr::transmute(
      item, 
      u0 = 0, v0 = 0,
      u1 = as.numeric(Rdisk * r_u1),
      v1 = as.numeric(Rdisk * r_u2)
    )
  
  p_disk <- ggplot2::ggplot() +
    ggplot2::geom_path(data = draw_disk_outline(), ggplot2::aes(x, y)) +
    ggplot2::geom_segment(
      data = S_disk,
      ggplot2::aes(x = u0, y = v0, xend = u1, yend = v1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_text_repel(
      data = S_disk,
      ggplot2::aes(x = u1, y = v1, label = item),
      size = 3.1, max.overlaps = Inf
    ) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    ggplot2::labs(x = "u1 (whitened b1,b2)", y = "u2") +
    theme_pub(12)
  
  list(p_base = p_base, p_disk = p_disk)
}

# Biplots
if (!is.null(varmap)) {
  Rtab <- build_biplot_data(Z_A, varmap, Base_A, geom$U)
  write_csv(Rtab, "items_vs_base_and_unitdisk_correlations.csv")
  
  BIP <- plot_biplots(Rtab, Base_A, geom$U)
  save_plot_gg("FIG_biplot_items_BASE", BIP$p_base, width = 8.0, height = 7.0)
  save_plot_gg("FIG_biplot_items_UNITDISK", BIP$p_disk, width = 8.0, height = 7.0)
  item_component_correlations_base <- as.matrix(Rtab[, c("r_b1", "r_b2")])
  rownames(item_component_correlations_base) <- Rtab$item
  colnames(item_component_correlations_base) <- c("u1", "u2")
} else {
  msgf("[biplot] varmap missing; skipping biplots.")
  item_component_correlations_base <- matrix(numeric(0), nrow = 0L, ncol = 2L)
  colnames(item_component_correlations_base) <- c("u1", "u2")
}

pc_scores_2d <- as.matrix(Base_A)
colnames(pc_scores_2d) <- c("u1", "u2")

saveRDS(list(
  participant_id = ids_base,
  pc_scores_2d = pc_scores_2d,
  spectrum = base_spectrum,
  explained_variance_ratio = base_explained_variance_ratio,
  selected_items = survivors,
  weights = w_full,
  item_component_correlations = item_component_correlations_base,
  weighting_mode = WEIGHTING_MODE,
  decomp_method = BASE_DECOMP_METHOD
), file = "method_sensitivity_fit_bundle.rds")
