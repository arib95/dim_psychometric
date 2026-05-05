# Shared cache and stochastic weight-search core for item-level Gower optimisation.

gower_opt_is_nominal_fast_path <- function(Xdf, type = NULL) {
  type_empty <- is.null(type) || !length(type) || all(lengths(type) == 0L)
  type_empty && all(vapply(Xdf, function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1)))
}

gower_opt_nominal_codes <- function(x) {
  factor_x <- if (is.factor(x) || is.logical(x)) {
    factor(as.character(x), exclude = NULL)
  } else {
    factor(x, exclude = NULL)
  }
  as.integer(factor_x)
}

gower_opt_make_nominal_mismatch_cache <- function(Xdf) {
  n <- nrow(Xdf)
  p <- ncol(Xdf)
  n_pairs <- n * (n - 1L) / 2L

  if (n < 2L) {
    return(list(N = vector("list", p), S = NULL, n = n, n_pairs = n_pairs, nominal_raw = TRUE))
  }

  pair_j <- rep.int(seq_len(n - 1L), times = (n - 1L):1L)
  pair_i <- sequence((n - 1L):1L) + pair_j

  N_list <- vector("list", p)
  for (j in seq_len(p)) {
    z <- gower_opt_nominal_codes(Xdf[[j]])
    N_list[[j]] <- as.raw(as.integer(z[pair_i] != z[pair_j]))
  }

  list(N = N_list, S = NULL, n = n, n_pairs = n_pairs, nominal_raw = TRUE)
}

gower_opt_make_ns_cache <- function(Xdf, type = NULL, gower_fun = NULL) {
  if (gower_opt_is_nominal_fast_path(Xdf, type)) {
    return(gower_opt_make_nominal_mismatch_cache(Xdf))
  }

  if (is.null(gower_fun)) {
    gower_fun <- function(Xdf, type, weights) {
      cluster::daisy(Xdf, metric = "gower", type = type, weights = weights)
    }
  }

  p <- ncol(Xdf)
  N_list <- vector("list", p)
  S_list <- vector("list", p)

  for (j in seq_len(p)) {
    w1 <- rep(0, p)
    w1[j] <- 1
    Dj <- gower_fun(Xdf, type, w1)
    vD <- as.numeric(Dj)
    ok <- as.numeric(is.finite(vD))
    vD[!is.finite(vD)] <- 0
    N_list[[j]] <- vD
    S_list[[j]] <- ok
  }

  list(N = N_list, S = S_list, n = nrow(Xdf), nominal_raw = FALSE)
}

gower_opt_cache_N <- function(cache, j) {
  if (isTRUE(cache$nominal_raw)) as.integer(cache$N[[j]]) else cache$N[[j]]
}

gower_opt_cache_S <- function(cache, j) {
  if (isTRUE(cache$nominal_raw)) 1 else cache$S[[j]]
}

gower_opt_initial_state <- function(cache, w) {
  if (isTRUE(cache$nominal_raw)) {
    num <- numeric(cache$n_pairs)
    for (j in seq_along(w)) {
      if (is.finite(w[[j]]) && w[[j]] != 0) {
        num <- num + w[[j]] * as.integer(cache$N[[j]])
      }
    }
    return(list(num = num, den = sum(w, na.rm = TRUE)))
  }

  list(
    num = Reduce(`+`, Map(`*`, cache$N, as.list(w))),
    den = Reduce(`+`, Map(`*`, cache$S, as.list(w)))
  )
}

gower_opt_dist_from_num_den <- function(num, den, n_rows) {
  Dvec <- num / pmax(den, .Machine$double.eps)
  attr(Dvec, "Size") <- n_rows
  attr(Dvec, "Diag") <- FALSE
  attr(Dvec, "Upper") <- FALSE
  class(Dvec) <- "dist"
  Dvec
}

gower_opt_calc_id <- function(num, den, n_rows, id_fun) {
  id_fun(gower_opt_dist_from_num_den(num, den, n_rows))
}

gower_opt_sample <- function(x, size, seed = NULL, seed_fun = NULL) {
  size <- min(length(x), size)
  if (!length(x) || size <= 0L) return(x[0])

  if (!is.null(seed) && !is.null(seed_fun)) {
    return(seed_fun(seed, function() x[sample.int(length(x), size)]))
  }

  x[sample.int(length(x), size)]
}

gower_opt_stochastic_weights <- function(cache,
                                         vars,
                                         init_weights,
                                         allow_update = rep(TRUE, length(init_weights)),
                                         id_fun,
                                         w_min = 0.05,
                                         step_grid = c(0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.05),
                                         batch_k = 3L,
                                         batch_factor = 0.75,
                                         max_iter = 1000L,
                                         eval_per_iter = 50L,
                                         seed_iter = NULL,
                                         seed_fun = NULL,
                                         progress_fun = NULL,
                                         verbose = FALSE,
                                         n_rows = cache$n) {
  w <- init_weights
  names(w) <- vars
  allow_update <- as.logical(allow_update)
  allow_update[is.na(allow_update)] <- FALSE
  w[!allow_update] <- pmax(w_min, w[!allow_update])

  state <- gower_opt_initial_state(cache, w)
  num_cur <- state$num
  den_cur <- state$den

  id <- gower_opt_calc_id(num_cur, den_cur, n_rows, id_fun)
  hist <- data.frame(iter = 0L, ID = id, changed = NA_character_, note = NA_character_)

  if (verbose) {
    cat(sprintf("[optim] Start ID: %.3f | N_sub: %d | Mode: STOCHASTIC\n", id, n_rows))
    if (isTRUE(cache$nominal_raw)) {
      cat(sprintf("[optim] Cache: nominal mismatch raw | pairs=%d | vars=%d\n", cache$n_pairs, length(cache$N)))
    }
  }

  max_iter_eff <- if (is.null(max_iter) || !is.finite(max_iter)) 1000L else as.integer(max_iter)
  eval_per_iter <- if (is.null(eval_per_iter) || !is.finite(eval_per_iter)) 50L else as.integer(eval_per_iter)
  eval_per_iter <- max(1L, eval_per_iter)
  step_grid <- sort(unique(step_grid), decreasing = FALSE)
  hot_vars <- integer(0)

  for (it in seq_len(max_iter_eff)) {
    can_all <- which(allow_update & (w > w_min + 1e-12))
    if (!length(can_all)) break

    n_rnd <- max(10L, eval_per_iter - length(hot_vars))
    rnd_vars <- gower_opt_sample(
      can_all,
      min(length(can_all), n_rnd),
      seed = if (!is.null(seed_iter)) as.integer(seed_iter + it) else NULL,
      seed_fun = seed_fun
    )
    can_iter <- unique(c(hot_vars, rnd_vars))

    best <- list(id = Inf)
    for (j in can_iter) {
      w_base <- w[[j]]
      for (factor in step_grid) {
        w_new <- max(w_min, w_base * factor)
        if (w_new >= w_base - 1e-12) next

        delta <- w_new - w_base
        id_try <- gower_opt_calc_id(
          num_cur + delta * gower_opt_cache_N(cache, j),
          den_cur + delta * gower_opt_cache_S(cache, j),
          n_rows,
          id_fun
        )

        if (is.finite(id_try) && id_try < best$id) {
          best <- list(id = id_try, j = j, w_new = w_new)
        }
      }
    }

    if (is.finite(best$id) && best$id < id - 1e-6) {
      jbest <- best$j
      wbest <- best$w_new
      delta <- wbest - w[[jbest]]

      num_cur <- num_cur + delta * gower_opt_cache_N(cache, jbest)
      den_cur <- den_cur + delta * gower_opt_cache_S(cache, jbest)
      w[[jbest]] <- wbest
      id <- best$id

      hot_vars <- unique(c(jbest, hot_vars))
      if (length(hot_vars) > 10L) hot_vars <- head(hot_vars, 10L)

      hist <- rbind(hist, data.frame(iter = it, ID = id, changed = vars[[jbest]], note = sprintf("%.3f", wbest)))
      if (verbose) cat(sprintf("   iter %d: %s -> %.3f (ID: %.3f)\n", it, vars[[jbest]], wbest, id))
    } else {
      if (length(hot_vars) > 0L) {
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

    if (batch_k > 1L && (it %% 5L == 0L)) {
      remain <- setdiff(can_all, best$j)
      remain <- gower_opt_sample(
        remain,
        min(length(remain), eval_per_iter),
        seed = if (!is.null(seed_iter)) as.integer(seed_iter + 10000L + it) else NULL,
        seed_fun = seed_fun
      )

      if (length(remain) > 0L) {
        scores <- numeric(length(remain))
        for (i in seq_along(remain)) {
          j <- remain[[i]]
          delta_b <- max(w_min, w[[j]] * batch_factor) - w[[j]]
          scores[[i]] <- gower_opt_calc_id(
            num_cur + delta_b * gower_opt_cache_N(cache, j),
            den_cur + delta_b * gower_opt_cache_S(cache, j),
            n_rows,
            id_fun
          )
        }

        take_vars <- remain[head(order(scores), min(batch_k, length(scores)))]

        if (length(take_vars)) {
          num_b <- num_cur
          den_b <- den_cur
          w_b <- w
          for (j in take_vars) {
            wn <- max(w_min, w[[j]] * batch_factor)
            d <- wn - w[[j]]
            num_b <- num_b + d * gower_opt_cache_N(cache, j)
            den_b <- den_b + d * gower_opt_cache_S(cache, j)
            w_b[[j]] <- wn
          }

          id_b <- gower_opt_calc_id(num_b, den_b, n_rows, id_fun)
          if (is.finite(id_b) && id_b < id - 1e-6) {
            num_cur <- num_b
            den_cur <- den_b
            w <- w_b
            id <- id_b
            hist <- rbind(hist, data.frame(iter = it, ID = id, changed = "BATCH", note = paste(length(take_vars), "vars")))
            if (verbose) cat(sprintf("   iter %d: [BATCH] x%.2f on %d vars (ID: %.3f)\n", it, batch_factor, length(take_vars), id))
          }
        }
      }
    }
  }

  list(weights = w, history = hist, final_ID = id)
}
