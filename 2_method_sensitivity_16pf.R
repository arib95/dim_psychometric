default_method_sensitivity_cfg <- list(
  out_dir = "out/method_sensitivity_16pf",
  data_path = "data/psychometric_matrix_without_missings.csv",
  diag_path = "data/wide_diagnoses.csv",
  qbank_path = "data/QBank.csv",
  scale_spec_path = "data/16pf_scale_spec.csv",
  setup_script = "0_setup.R",
  runner_script = "tools/run_dimension_variant.R",
  summary_script = "tools/post_method_sensitivity_16pf.R",
  bam_threads = "all_but_one",
  ncores_par = "all_but_one",
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
  post_parallel_diagnostics = TRUE,
  spectrum_rank = 10L,
  baseline_reps = 500L,
  baseline_subset_modes = c("unstratified", "stratified_by_scale"),
  distance_eval_n = 1200L,
  distance_eval_k = 15L,
  future_scheduling = 1L,
  save_plots = TRUE,
  run_fits = TRUE,
  run_post = TRUE,
  overwrite = FALSE,
  run_residual_diagnostics = FALSE
)

messagef <- function(...) {
  message(sprintf(...))
}

infer_root_dir <- function() {
  cwd <- normalizePath(".", mustWork = TRUE)
  if (file.exists(file.path(cwd, "0_setup.R"))) {
    return(cwd)
  }
  stop("Run from the project root or pass root_dir explicitly.")
}

normalise_path <- function(path, root_dir, must_work = TRUE) {
  if (!grepl("^(/|~|[A-Za-z]:)", path)) {
    path <- file.path(root_dir, path)
  }
  normalizePath(path, mustWork = must_work)
}

source_script_env <- function(path) {
  env <- new.env(parent = globalenv())
  source(path, local = env, chdir = FALSE)
  env
}

method_specs <- tibble::tribble(
  ~weighting_mode, ~decomp_method, ~run_dir,
  "uniform", "pca", "pca_uniform",
  "uniform", "robust_pca", "robpca_uniform",
  "id_guided", "pca", "pca_id_guided",
  "id_guided", "robust_pca", "robpca_id_guided"
)

run_one_variant <- function(spec, cfg, root_dir, runner_env) {
  run_dir <- file.path(cfg$out_dir, spec$run_dir[[1]])
  bundle_path <- file.path(run_dir, "method_sensitivity_fit_bundle.rds")
  label <- spec$run_dir[[1]]
  
  if (!isTRUE(cfg$overwrite) && file.exists(bundle_path)) {
    messagef("Skipping %s; bundle already exists at %s", label, bundle_path)
    return(invisible(bundle_path))
  }
  
  messagef("Running %s -> %s", label, run_dir)
  run_cfg <- list(
    root_dir = root_dir,
    out_dir = run_dir,
    weighting_mode = spec$weighting_mode[[1]],
    decomp_method = spec$decomp_method[[1]],
    psych_csv = cfg$data_path,
    diag_csv = cfg$diag_path,
    setup_script = file.path(root_dir, "0_setup.R"),
    main_script = file.path(root_dir, "1_dimension_psychometric.R"),
    bam_threads = cfg$bam_threads,
    ncores_par = cfg$ncores_par,
    omp_threads = cfg$omp_threads,
    blas_threads = cfg$blas_threads,
    optimisation_subsample_n = cfg$optimisation_subsample_n,
    optimisation_w_min = cfg$optimisation_w_min,
    optimisation_multi_runs = cfg$optimisation_multi_runs,
    optimisation_multi_min_prop = cfg$optimisation_multi_min_prop,
    optimisation_step_grid = cfg$optimisation_step_grid,
    optimisation_batch_k = cfg$optimisation_batch_k,
    optimisation_batch_factor = cfg$optimisation_batch_factor,
    optimisation_max_iter = cfg$optimisation_max_iter,
    optimisation_eval_per_iter = cfg$optimisation_eval_per_iter,
    run_residual_diagnostics = cfg$run_residual_diagnostics
  )
  runner_env$run_dimension_variant_main(run_cfg)
  invisible(bundle_path)
}

run_method_sensitivity_16pf <- function(cfg = default_method_sensitivity_cfg) {
  cfg <- utils::modifyList(default_method_sensitivity_cfg, cfg)
  root_dir <- cfg$root_dir %||% infer_root_dir()
  cfg$out_dir <- normalise_path(cfg$out_dir, root_dir, must_work = FALSE)
  cfg$data_path <- normalise_path(cfg$data_path, root_dir, must_work = TRUE)
  cfg$diag_path <- normalise_path(cfg$diag_path, root_dir, must_work = FALSE)
  cfg$qbank_path <- normalise_path(cfg$qbank_path, root_dir, must_work = TRUE)
  cfg$scale_spec_path <- normalise_path(cfg$scale_spec_path, root_dir, must_work = TRUE)
  cfg$setup_script <- normalise_path(cfg$setup_script, root_dir, must_work = TRUE)
  cfg$runner_script <- normalise_path(cfg$runner_script, root_dir, must_work = TRUE)
  cfg$summary_script <- normalise_path(cfg$summary_script, root_dir, must_work = TRUE)
  
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  
  runner_env <- source_script_env(cfg$runner_script)
  summary_env <- source_script_env(cfg$summary_script)
  
  if (isTRUE(cfg$run_fits)) {
    for (i in seq_len(nrow(method_specs))) {
      run_one_variant(method_specs[i, , drop = FALSE], cfg, root_dir, runner_env)
    }
  }
  
  if (isTRUE(cfg$run_post)) {
    messagef("Running method summaries in %s", cfg$out_dir)
    summary_env$run_method_sensitivity_summary(list(
      root_dir = root_dir,
      out_dir = cfg$out_dir,
      data_path = cfg$data_path,
      qbank_path = cfg$qbank_path,
      scale_spec_path = cfg$scale_spec_path,
      setup_script = cfg$setup_script,
      bam_threads = cfg$bam_threads,
      ncores_par = cfg$ncores_par,
      omp_threads = cfg$omp_threads,
      blas_threads = cfg$blas_threads,
      seed = cfg$seed,
      min_scale_completion_prop = cfg$min_scale_completion_prop,
      cv_folds = cfg$cv_folds,
      optimisation_subsample_n = cfg$optimisation_subsample_n,
      optimisation_w_min = cfg$optimisation_w_min,
      optimisation_multi_runs = cfg$optimisation_multi_runs,
      selection_stability_optimisation_multi_runs = cfg$selection_stability_optimisation_multi_runs,
      optimisation_multi_min_prop = cfg$optimisation_multi_min_prop,
      optimisation_step_grid = cfg$optimisation_step_grid,
      optimisation_batch_k = cfg$optimisation_batch_k,
      optimisation_batch_factor = cfg$optimisation_batch_factor,
      optimisation_max_iter = cfg$optimisation_max_iter,
      optimisation_eval_per_iter = cfg$optimisation_eval_per_iter,
      stability_reps = cfg$stability_reps,
      stability_fixed_selection_reps = cfg$stability_fixed_selection_reps,
      stability_half_n = cfg$stability_half_n,
      split_half_full_refit_parallel = cfg$split_half_full_refit_parallel,
      split_half_data_cache = cfg$split_half_data_cache,
      post_parallel_diagnostics = cfg$post_parallel_diagnostics,
      spectrum_rank = cfg$spectrum_rank,
      baseline_reps = cfg$baseline_reps,
      baseline_subset_modes = cfg$baseline_subset_modes,
      distance_eval_n = cfg$distance_eval_n,
      distance_eval_k = cfg$distance_eval_k,
      future_scheduling = cfg$future_scheduling,
      save_plots = cfg$save_plots
    ))
  }
  
  messagef("Finished. Outputs written to %s", cfg$out_dir)
  invisible(cfg$out_dir)
}

run_method_sensitivity_16pf(default_method_sensitivity_cfg)
