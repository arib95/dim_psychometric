# =============================== Unified Setup ================================
suppressPackageStartupMessages({
  # Base + utils
  library(stats)
  library(utils)
  library(grDevices)
  library(Matrix)
  
  # Data
  library(readr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(zoo)
  library(tibble)
  
  # Maths / ML
  library(cluster)
  library(RANN)
  library(FNN)
  library(RSpectra)
  library(dbscan)
  library(glmnet)
  library(MASS)
  library(expm)
  library(clue)
  library(rrcov)
  
  # Graphs / viz
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(webshot2)
  library(patchwork)
  
  # Domain / extras
  library(igraph)
  library(aricode)
  library(vegan)
  library(princurve)
  library(FactoMineR)
  library(reticulate)
  library(R.utils)
  library(pROC)
  library(PRROC)
  library(robustbase)
  library(isotone)
  library(grid)
  library(metR)
  
  # parallel
  library(withr)
  library(future)
  library(future.apply)
  library(parallelly)
  library(progressr)
  library(RhpcBLASctl)
})

# =========================== User-defined defaults ============================

# General / logging / outputs
VERBOSE        <- FALSE
OUTPUTS_DIR    <- "out"
PNG_DPI        <- 300L
SAVE_IMAGES    <- TRUE
SAVE_PLOT_RDS  <- TRUE

# RNG / reproducibility
RNG_KIND     <- "L'Ecuyer-CMRG"
SEED         <- 42L
SEED_GLOBAL  <- 42L
BUNDLE_SEED  <- SEED_GLOBAL
SEED_PRED    <- SEED_GLOBAL
SEED_JITTER  <- SEED_GLOBAL
SEED_BOOT    <- SEED_GLOBAL
OOF_SEED     <- SEED_GLOBAL

# Compute / threading
default_setup_cfg <- list(
  BAM_THREADS = "all_but_one",
  SET_ENV_THREADS = TRUE,
  OMP_THREADS = 1L,
  BLAS_THREADS = 1L,
  NCORES_PAR = "all_but_one",
  FUTURE_GLOBALS_MAXSIZE = 1024^3
)
setup_cfg <- utils::modifyList(default_setup_cfg, SETUP_CFG)
BAM_THREADS <- setup_cfg$BAM_THREADS
SET_ENV_THREADS <- setup_cfg$SET_ENV_THREADS
OMP_THREADS <- setup_cfg$OMP_THREADS
BLAS_THREADS <- setup_cfg$BLAS_THREADS
NCORES_PAR <- setup_cfg$NCORES_PAR
FUTURE_GLOBALS_MAXSIZE <- setup_cfg$FUTURE_GLOBALS_MAXSIZE

# Palette / plotting defaults
PALETTE_ENGINE     <- "scico"
PALETTE_NAME       <- "lipari"   # sequential
PALETTE_NAME_DIV   <- "vik"      # diverging
PALETTE_DIRECTION  <- 1L
PALETTE_COLOURS    <- NULL

# Global toggles
DO_PLOTS        <- TRUE
DO_DIAGNOSTICS  <- TRUE
DO_SURFACE      <- TRUE
DO_SWEEP        <- TRUE
RUN_PRINSURF    <- TRUE
USE_SURFACE_CV  <- TRUE
DO_DELTAQ       <- TRUE
DO_CORR_TRIM    <- FALSE

# Embedding (clusters)
KNN_K        <- 9L
KNN_VARIANT  <- "mutual"   # mutual | union
LOCAL_SCALE  <- TRUE
MULT_WEIGHT  <- TRUE

# Embedding (psychometric)
M_STAR_FIXED  <- 2L
K_ID_LO_HI    <- c(8L, 20L)
PC_MAX        <- 6L
M_DEFAULT     <- 4L
BASE_DECOMP_METHOD <- "robust_pca"  # robust_pca | pca

# Community detection
COMMUNITY_ALGO      <- "leiden"
LEIDEN_OBJECTIVE    <- "modularity"
LEIDEN_GAMMA        <- 0.75
LEIDEN_ITERS        <- -1L
MIN_CLUSTER_SIZE    <- 2L   # on unique profiles
MIN_CLUSTER_WEIGHT  <- 8L   # expanded weight

# Surface
K_BASIS_UV        <- 40L
SURFACE_CV_FOLDS  <- 3L

# Bundle-fibre
BF_NGR_UV      <- 35L
BF_K_BASE      <- 100L
BF_K_RANK_CAP  <- 4L
BF_MP_SHRINK   <- 0.90
BF_MIN_N_NEIGH <- 50L
BF_DENS_K      <- 30L

# Nulls / resampling
N_PERM            <- 200L
DEG_REPS          <- N_PERM
NULL_SCOPE        <- "kept"    # kept | full
N_BOOT            <- 200L
SIL_REPS          <- 200L
COL_SHUFFLE_REPS  <- 200L
CV_FOLDS          <- 5L
SIG_METHOD        <- "perm"
SIG_FWER          <- 0.95
PRED_BOOT         <- 200L

# Features / guards
INCLUDE_NODIAG   <- TRUE
DX_MIN_PREV      <- 0.00
DX_MAX_PREV      <- 0.99
MIN_CASES_TOTAL  <- 10L
MIN_CASES_IN     <- 5L
MIN_CASES_OUT    <- 5L
SHOW_NS_LABELS   <- TRUE

# Targets / features (additional)
PREF_TARGET         <- "ANY_DX"
CORR_THRESH         <- 0.95
DX_DENY_NOS         <- FALSE
DX_PREV_MIN         <- 0.00
DX_CASES_MIN        <- 10L
N_TOP_PER_DX        <- 80L
RARE_LEVEL_MIN_PROP <- 0.01

# Pillars (A, B, C)
ALPHA_FDR      <- 0.05
MIN_PREV_IN_CL <- 0.00
MIN_OR         <- 2.0
ALPHA_LOCALIZE <- 0.05
AUC_MIN        <- 0.70
PREV_MIN       <- 0.03
NCASE_MIN      <- 5L
NIN_MIN        <- 5L
NOUT_MIN       <- 5L
MAJORS_BOOT_B  <- 200L
PILLAR_B_REPS  <- 200L

# Clusters
K_GRID       <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
VARIANT_GRID <- c("union", "mutual")
SWEEP_BOOT_B <- 200L

# Weights
WEIGHTING_MODE        <- "id_guided"  # id_guided | uniform
TREAT_ORDINALS_AS_NOMINAL <- TRUE
MISSING_AS_NOMINAL_LEVEL  <- TRUE
W_MIN                <- 0.05
W_STEP_GRID          <- c(0.95, 0.90, 0.75, 0.5, 0.25, 0.10, 0.05)
W_BATCH_K            <- 3L
W_BATCH_FACTOR       <- 0.75
W_MAX_ITERS          <- NA_integer_
N_ROWS_SUB           <- 200L
FIX_REP_SUBSET       <- TRUE
GOWER_MULTI_RUNS     <- 10L
GOWER_MULTI_MIN_PROP <- 0.35
GOWER_MULTI_ENABLE   <- TRUE
GOWER_ACTIVE_EPS     <- 1e-8

# Dedup / diagnostics
DO_DEDUP         <- TRUE
EPS_DEDUP        <- NA_real_
WRITE_DEDUP_CSV  <- TRUE
WRITE_DEGENERATE_CSV <- TRUE
CORE_BAND        <- c(0.20, 0.70)  # effective default in original script
CORE_KNN_K       <- 10L            # effective default in original script
SIG_Q            <- 0.01
K_FIBRE_CAP      <- 3L
MIN_N_DX         <- 25L

# Dedup mode:
#   - "gower_complete": current exact complete-linkage dedup in full Gower space
#   - "hash_exact": exact hash on the prepped, non-jittered representation
#   - "hash_round": rounded hash on numeric columns of the prepped, non-jittered representation
#   - "none": no deduplication
DEDUP_MODE        <- "none"
DEDUP_HASH_DIGITS <- 4L
DEDUP_HASH_NA     <- "NA"

# Final Gower diagnostics after variable selection:
#   - "full": use all retained rows (expensive; not suitable for large n)
#   - "reps": use all deduplicated representatives
#   - "sample_reps": sample from deduplicated representatives
#   - "sample_all": sample from all retained rows
#   - "none": skip final distance-based diagnostics entirely
FINAL_DIAG_MODE <- "sample_reps"
DIAG_N_MAX      <- 4000L

# Core sweep / optimisation
PCL_K_MIN       <- 2L
PCL_K_MAX       <- 4L
PCL_LAMBDA_GRID <- c(0, 0.05, 0.1, 0.2, 0.4, 0.8)
PCL_KNN_BASE    <- 12L
PCL_LOCAL_SCALE <- TRUE
PCL_N_START     <- 3L
PCL_N_ITER      <- 700L
PCL_TOL         <- 1e-5
PCL_EPS         <- 1e-10
PCL_MIN_POS_DEF <- 10L
PCL_MIN_NEG_DEF <- 10L

# Usage / regularisation
PCL_MASS_FLOOR_FRAC <- 0.05
PCL_BALANCE_EVERY   <- 8L
PCL_RESEED_MAX      <- 3L
PCL_ENTROPY_PUSH    <- 1e-3

# Selection weights / penalties
PCL_W_Q       <- 0.45
PCL_W_SNN     <- 0.50
PCL_W_REC     <- 0.20
PCL_W_SIL     <- 0.15
PCL_PEN_EMPTY <- 0.60
PCL_PEN_IMB   <- 0.15
PCL_PEN_RED   <- 0.30

# Leiden
LEIDEN_MEMBERSHIP_CSV <- "cluster_membership_all_participants.csv"
USE_LEIDEN_K          <- FALSE

# Plot behaviour
USE_RAGG      <- FALSE
MM_UNITS      <- FALSE
GRID_N_B      <- 400L
GRID_N_U      <- 400L
CONTOUR_PROBS <- c(0.15, 0.25, 0.40, 0.60, 0.80)
CALIB_SHOWPTS <- TRUE
CL_DRAW_HULLS <- TRUE

# k-NN stability plot
STAB_K_RANGE <- c(10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L)
STAB_SD_GRID <- c(0, 0.05, 0.10, 0.15, 0.20)
STAB_REPS    <- 600L

# Outcome plots
behaviour_csv <- "data/wide_diagnoses.csv"
OUT_SUBDIR    <- "wide_diagnoses"
# behaviour_csv <- "data/cluster_membership_all_participants.csv"
# OUT_SUBDIR    <- "clusters"
# behaviour_csv <- "data/hopkins.csv"
# OUT_SUBDIR    <- "hopkins"

vars_keep   <- NULL
write_plots <- DO_PLOTS

MAX_K       <- 8L
MIN_POS     <- 15L
FIT_TIMEOUT <- 60.0

FDR_METHOD <- "BH"   # BH | BY
DO_BOOT    <- TRUE
BOOT_B     <- 199L

# OOF significance config (AUC/C with paired bootstrap)
OOF_BOOT_B    <- 400L
OOF_MAX_PAIRS <- 2e6
OOF_PERM_B    <- 200L

# Inputs
PSY_CSV     <- "data/psychometric_matrix.csv"
DX_OPTIONAL <- TRUE
DIAG_CSV    <- "data/wide_diagnoses.csv"

# ==================== Build mgcv with OpenMP in Apple Silicon =================
# There is experimental support for building mgcv with OpenMP, inc. for Apple
# Silicon Macs. A useful Stack Overflow thread can be found at:
# https://stackoverflow.com/questions/78590558/how-can-i-enable-openmp-in-mac-os-for-mgcv-r-package
library(mgcv)

# =============================== Helpers =====================================
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

.is_relpath  <- function(p) !grepl("^(/|[A-Za-z]:)", p) && dirname(p) == "."
.is_abs_path <- function(p) grepl("^(/|~|[A-Za-z]:)", p)
.is_output_relative <- function(p, out_dir = OUTPUTS_DIR) {
  p_norm <- gsub("\\\\", "/", sub("^\\./", "", p))
  out_norm <- gsub("\\\\", "/", sub("^\\./", "", out_dir))
  identical(p_norm, out_norm) || startsWith(p_norm, paste0(out_norm, "/"))
}
.resolve_output_path <- function(file) {
  if (.is_abs_path(file) || .is_output_relative(file)) file else file.path(OUTPUTS_DIR, file)
}

safe_file <- function(s) gsub("[^A-Za-z0-9_]+", "_", s)

.timeit <- function(expr) {
  t0 <- proc.time()[[3]]
  val <- eval.parent(substitute(expr))
  list(value = val, elapsed = proc.time()[[3]] - t0)
}

.is_binary <- function(x) is.numeric(x) && all(x %in% c(0, 1, NA))
.check_binary_df <- function(df) stopifnot(all(vapply(df, .is_binary, logical(1))))

# =============================== Outputs & I/O ================================
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  if (isTRUE(VERBOSE)) message(sprintf(...))
}

write_csv <- function(x, file, ...) {
  file <- .resolve_output_path(file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv2(x, file, ...)
}

saveRDS <- function(object, file, ...) {
  file <- .resolve_output_path(file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  base::saveRDS(object, file = file, ...)
}

read_csv <- function(file, ...) {
  file <- .resolve_output_path(file)
  as.data.frame(
    readr::read_csv2(file, progress = FALSE, show_col_types = FALSE, ...),
    stringsAsFactors = FALSE
  )
}

readRDS <- function(file, ...) {
  file <- .resolve_output_path(file)
  base::readRDS(file, ...)
}

clean_ggplot_env <- function(p) {
  if (inherits(p, "ggplot")) {
    p$plot_env <- new.env(parent = emptyenv())
  }
  p
}

save_plot_gg <- function(name, plot, width, height, dpi = PNG_DPI, save_rds = FALSE) {
  if (!SAVE_IMAGES && !SAVE_PLOT_RDS) return(invisible())
  
  base_path <- file.path(OUTPUTS_DIR, name)
  png_path  <- paste0(base_path, ".png")
  pdf_path  <- paste0(base_path, ".pdf")
  rds_path  <- paste0(base_path, ".rds")
  
  if (SAVE_IMAGES) {
    if (isTRUE(USE_RAGG) && requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(
        png_path, plot = plot,
        width = width, height = height, dpi = dpi,
        device = ragg::agg_png
      )
    } else {
      ggplot2::ggsave(
        png_path, plot = plot,
        width = width, height = height, dpi = dpi
      )
    }
    
    ggplot2::ggsave(
      pdf_path, plot = plot,
      width = width, height = height,
      device = grDevices::cairo_pdf
    )
  }
  
  if (save_rds) {
    plot_clean <- clean_ggplot_env(plot)
    base::saveRDS(plot_clean, file = rds_path)
  }
  
  invisible(plot)
}

save_plot_base <- function(name, width, height, dpi = PNG_DPI, expr) {
  if (!SAVE_IMAGES) return(invisible())
  png_path <- file.path(OUTPUTS_DIR, paste0(name, ".png"))
  pdf_path <- file.path(OUTPUTS_DIR, paste0(name, ".pdf"))
  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE); eval.parent(substitute(expr)); grDevices::dev.off()
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE); eval.parent(substitute(expr)); grDevices::dev.off()
}

# =============================== RNG & Threads ================================
do.call(RNGkind, as.list(strsplit(RNG_KIND, ",")[[1]]))
set.seed(SEED_GLOBAL)

.hash32 <- function(s, mod = 2147483647) {
  s <- enc2utf8(as.character(s))
  b <- as.integer(charToRaw(s))
  h <- 0
  mul <- 65599
  for (x in b) h <- (h * mul + x) %% mod
  as.integer(h)
}

.seed_from_key <- function(seed, key, mod = 2147483647) {
  seed <- as.integer(seed)
  as.integer((seed + .hash32(key, mod = mod)) %% mod)
}

.with_seed <- function(seed, expr) {
  expr <- substitute(expr)
  if (is.null(seed) || !is.finite(seed)) return(eval(expr, parent.frame()))
  seed <- as.integer(seed)
  
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  
  set.seed(seed)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  eval(expr, parent.frame())
}

if (is.character(BAM_THREADS) && identical(BAM_THREADS, "all_but_one")) {
  BAM_THREADS <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
} else {
  BAM_THREADS <- as.integer(BAM_THREADS)
}
options(mc.cores = BAM_THREADS)

if (isTRUE(SET_ENV_THREADS)) {
  Sys.setenv(
    OMP_NUM_THREADS        = as.character(as.integer(OMP_THREADS)),
    OPENBLAS_NUM_THREADS   = as.character(as.integer(BLAS_THREADS)),
    MKL_NUM_THREADS        = as.character(as.integer(BLAS_THREADS)),
    VECLIB_MAXIMUM_THREADS = as.character(as.integer(BLAS_THREADS)),
    BLAS_NUM_THREADS       = as.character(as.integer(BLAS_THREADS))
  )
}

if (is.character(NCORES_PAR) && identical(NCORES_PAR, "all_but_one")) {
  NCORES_PAR <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
} else {
  NCORES_PAR <- as.integer(NCORES_PAR %||% max(1L, parallel::detectCores(logical = FALSE) - 1L))
}

# ==================================== Plots ===================================
scale_shared_div <- function(name = NULL, limits = NULL) {
  if (requireNamespace("scico", quietly = TRUE)) {
    scico::scale_fill_scico(
      palette   = PALETTE_NAME_DIV,
      direction = PALETTE_DIRECTION,
      limits    = limits,
      oob       = scales::squish,
      name      = name
    )
  } else {
    ggplot2::scale_fill_gradient2(
      limits = limits,
      oob    = scales::squish,
      name   = name
    )
  }
}

scale_prob_fill <- function(limits = NULL, name = NULL) {
  eng <- tolower(PALETTE_ENGINE)
  dir <- PALETTE_DIRECTION
  pal <- PALETTE_NAME
  
  if (eng == "scico" && requireNamespace("scico", quietly = TRUE)) {
    scico::scale_fill_scico(
      palette = pal,
      direction = dir,
      limits = limits,
      oob = scales::squish,
      name = name
    )
  } else if (eng == "colorspace" && requireNamespace("colorspace", quietly = TRUE)) {
    colorspace::scale_fill_continuous_sequential(
      palette = pal,
      rev = (dir == -1),
      limits = limits,
      name = name
    )
  } else if (eng == "brewer") {
    ggplot2::scale_fill_distiller(
      palette = pal,
      type = "seq",
      direction = if (dir == 1) 1 else -1,
      limits = limits,
      oob = scales::squish,
      name = name
    )
  } else if (eng == "paletteer" && requireNamespace("paletteer", quietly = TRUE)) {
    paletteer::scale_fill_paletteer_c(
      pal,
      direction = dir,
      limits = limits,
      name = name
    )
  } else if (eng == "manual" && !is.null(PALETTE_COLOURS)) {
    ggplot2::scale_fill_gradientn(
      colours = PALETTE_COLOURS,
      limits = limits,
      name = name
    )
  } else {
    ggplot2::scale_fill_gradient(
      low = "#f7fbff",
      high = "#08306b",
      limits = limits,
      name = name
    )
  }
}

cluster_colours <- function(levels_vec) {
  if (requireNamespace("scico", quietly = TRUE)) {
    vals <- scico::scico(
      length(levels_vec),
      palette = PALETTE_NAME,
      direction = PALETTE_DIRECTION
    )
    stats::setNames(vals, levels_vec)
  } else {
    cols <- grDevices::hcl.colors(length(levels_vec), "Dark 3")
    stats::setNames(cols, levels_vec)
  }
}

theme_pub <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 1)
    )
}

theme_marginal <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linetype = "dotted",
        linewidth = 0.25
      ),
      panel.grid.minor = ggplot2::element_line(
        colour = "grey92",
        linetype = "dotted",
        linewidth = 0.25
      ),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(3, "pt"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0, size = ggplot2::rel(1.4)
      ),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      legend.position = "none",
      plot.margin = grid::unit(c(5.5, 5.5, 5.5, 5.5), "pt")
    )
}

# ============================= Parallel Computing =============================
options(
  future.rng.onMisuse    = "warn",
  future.globals.maxSize = as.numeric(FUTURE_GLOBALS_MAXSIZE)
)

options(progressr.enable = TRUE)
progressr::handlers("txtprogressbar")
progressr::handlers(global = TRUE)

NWORKERS <- as.integer(NCORES_PAR %||% 1L)

ok <- try({
  if (NWORKERS <= 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = NWORKERS)
  }
}, silent = TRUE)

if (inherits(ok, "try-error")) {
  message("[future] fallback to sequential: ", conditionMessage(attr(ok, "condition")))
  future::plan(future::sequential)
}

.drop_future_args <- function(dots) {
  if (is.null(names(dots))) return(dots)
  dots[!grepl("^future\\.", names(dots), perl = TRUE)]
}

FUTURE_LAPPLY <- function(X, FUN, ...) {
  dots <- list(...)
  ok <- try(do.call(future.apply::future_lapply, c(list(X, FUN), dots)), silent = TRUE)
  if (!inherits(ok, "try-error")) return(ok)
  message("[parallel] future_lapply failed → base::lapply: ",
          conditionMessage(attr(ok, "condition")))
  do.call(base::lapply, c(list(X, FUN), .drop_future_args(dots)))
}

FUTURE_SAPPLY <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  dots <- list(...)
  ok <- try(
    do.call(
      future.apply::future_sapply,
      c(list(X, FUN, simplify = simplify, USE.NAMES = USE.NAMES), dots)
    ),
    silent = TRUE
  )
  if (!inherits(ok, "try-error")) return(ok)
  message("[parallel] future_sapply failed → base::sapply: ",
          conditionMessage(attr(ok, "condition")))
  do.call(
    base::sapply,
    c(list(X, FUN, simplify = simplify, USE.NAMES = USE.NAMES), .drop_future_args(dots))
  )
}

FUTURE_MAPPLY <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE) {
  dots <- list(...)
  ok <- try(
    do.call(
      future.apply::future_mapply,
      c(list(FUN), dots,
        list(MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
    ),
    silent = TRUE
  )
  if (!inherits(ok, "try-error")) return(ok)
  message("[parallel] future_mapply failed, falling back to mapply: ",
          conditionMessage(attr(ok, "condition")))
  do.call(
    base::mapply,
    c(list(FUN), .drop_future_args(dots),
      list(MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
  )
}

# ========================= Resave plot as RDS =================================
resave_plot <- function(name,
                        width,
                        height,
                        dpi = PNG_DPI,
                        types = c("png", "pdf"),
                        suffix = NULL) {
  p <- readRDS(file.path(paste0(name, ".rds")))
  
  if (!is.null(suffix) && nzchar(suffix)) {
    name_out <- paste0(name, "_", suffix)
  } else {
    name_out <- paste0(name, "_", width, "x", height)
  }
  
  base_path <- file.path(OUTPUTS_DIR, name_out)
  
  if ("png" %in% types) {
    ggplot2::ggsave(
      paste0(base_path, ".png"),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  if ("pdf" %in% types) {
    ggplot2::ggsave(
      paste0(base_path, ".pdf"),
      plot = p,
      width = width,
      height = height,
      device = grDevices::cairo_pdf
    )
  }
  
  invisible(base_path)
}

# resave_plot("plotXYZ", width = 4, height = 6, suffix = "tall/wide/square")
