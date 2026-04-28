default_run_cfg <- list(
  root_dir = NULL,
  out_dir = "out/method_sensitivity_16pf/pca_uniform",
  weighting_mode = "uniform",
  decomp_method = "pca",
  psych_csv = "data/psychometric_matrix.csv",
  diag_csv = "data/wide_diagnoses.csv",
  setup_script = "0_setup.R",
  main_script = "1_dimension_psychometric.R",
  bam_threads = "all_but_one",
  ncores_par = "all_but_one",
  omp_threads = 1L,
  blas_threads = 1L
)

`%||%` <- function(x, y) if (is.null(x)) y else x

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

normalise_decomp_method <- function(x) {
  x <- tolower(x)
  if (identical(x, "robpca")) return("robust_pca")
  x
}

run_dimension_variant_main <- function(cfg = default_run_cfg) {
  cfg <- utils::modifyList(default_run_cfg, cfg)
  cfg$weighting_mode <- tolower(cfg$weighting_mode)
  cfg$decomp_method <- normalise_decomp_method(cfg$decomp_method)

  if (!cfg$weighting_mode %in% c("uniform", "id_guided")) {
    stop("Unknown weighting_mode: ", cfg$weighting_mode)
  }
  if (!cfg$decomp_method %in% c("pca", "robust_pca")) {
    stop("Unknown decomp_method: ", cfg$decomp_method)
  }

  root_dir <- cfg$root_dir %||% infer_root_dir()
  setup_path <- normalise_path(cfg$setup_script, root_dir, must_work = TRUE)
  main_path <- normalise_path(cfg$main_script, root_dir, must_work = TRUE)
  out_dir <- normalise_path(cfg$out_dir, root_dir, must_work = FALSE)
  psych_csv <- normalise_path(cfg$psych_csv, root_dir, must_work = TRUE)
  diag_csv <- normalise_path(cfg$diag_csv, root_dir, must_work = FALSE)

  SETUP_CFG <- list(
    BAM_THREADS = cfg$bam_threads,
    NCORES_PAR = cfg$ncores_par,
    OMP_THREADS = cfg$omp_threads,
    BLAS_THREADS = cfg$blas_threads,
    SET_ENV_THREADS = TRUE
  )
  source(setup_path, chdir = TRUE, local = TRUE)

  OUTPUTS_DIR <- out_dir
  WEIGHTING_MODE <- cfg$weighting_mode
  BASE_DECOMP_METHOD <- cfg$decomp_method
  PSY_CSV <- psych_csv
  DIAG_CSV <- diag_csv
  dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(OUTPUTS_DIR)

  message(sprintf(
    "[run_dimension_variant] out=%s | weighting=%s | decomp=%s | bam_threads=%d | workers=%d",
    OUTPUTS_DIR,
    WEIGHTING_MODE,
    BASE_DECOMP_METHOD,
    BAM_THREADS,
    NWORKERS
  ))

  source(main_path, chdir = FALSE, local = TRUE)
  invisible(OUTPUTS_DIR)
}
