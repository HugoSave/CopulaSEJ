source("dev/dev_utils.R")

compare_decoupler_margin_estimations <- function(studies, settings, seed=42) {
  # We don't need to filter for questions with more than 10 questions because we are not doing any copula fitting.
  combinations <- tidyr::expand_grid(study=studies, settings)
  set.seed(seed)
  # combinations[19,]
  res <- purrr::pmap(combinations, \(study, error_estimation_settings, decoupler) {
    get_obj_func <- purrr::partial(fit_and_construct_posterior_indep,
                   error_estimation_settings = error_estimation_settings,
                   copula_model = "indep",
                   q_support_restriction=NULL,
                   q_support_overshoot=0.1,
                   rejection_threshold=NULL)
    res <- evalute_marginal_fit(study, decoupler, get_obj_func)
    res$settings <- decoupler_margin_estimation_settings_to_shortname(error_estimation_settings)
    res$decoupler <- decoupler_name_with_settings(decoupler)
    res
  }, .progress="studies&margin")
  res <- res |> purrr::list_rbind()
  res
}

decoupler_name_with_settings <- function(decoupler) {
  if (inherits(decoupler, "relative_decoupler")) {
    decoupler_name <- paste0("Relative_mu_G", ":", decoupler$k)
  } else if (inherits(decoupler, "CDF_decoupler")) {
    decoupler_name <-"CDF"
  } else {
    stop("Unknown decoupler type.")

  }

}

decoupler_margin_estimation_settings_to_shortname <- function(settings) {
  # settings =list(
  #  method="beta_MAP",
  #  beta_mean=0.5,
  #  beta_var=1/12,
  #  prior_std=1
  #)
  settings$method <- switch(
    settings$method,
    "beta_MAP" = "MAP",
    "beta_MLE" = "MLE",
    settings$method # default value
  )
  #settings$prior_std <- switch(
  #  settings$prior_std
  #)

  prior_std <- if (!is.null(settings$prior_std)) paste0(":", round(settings$prior_std, digits=4)) else ""
  glue::glue(
    "{settings$method}{prior_std}"
  )
}

decoupler_and_margin_estimation_settings_study <- function() {
  decouplers <- list(
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=1.5),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=1),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.5),
    #get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.5),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.1),
    #get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.05),
    #get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.01),
    get_CDF_decoupler()
  )
  error_estimation_settings_list <- list(
      list(
        method="beta_MAP",
        prior_std=1
      ),
      list(
        method="beta_MAP",
        prior_std=0.5
      ),
      list(
        method="beta_MAP",
        prior_std=0.1
      ),
      list(
        method="beta_MAP",
        prior_std=0.05
      ),
     list(
       method="beta_MAP",
       prior_std=0.01
     ),
    list(
      method="uniform"
    ),
    list(
      method="beta_MLE"
    ),
    list(
      method="beta_prior"
    )
  )
  combinations <- tidyr::expand_grid(error_estimation_settings=error_estimation_settings_list, decoupler=decouplers )
  combinations
}

evalute_marginal_fit <-  function(study_data, decoupler, get_posterior_obj, k_percentiles=c(5,50,95)) {
  fold_combinations <- create_cross_validation_sets(study_data)

  stats <- purrr::pmap(fold_combinations, function(training_set, test_set) {
    arr_format <- df_format_to_array_format(training_set, test_set, get_three_quantiles_summarizing_function()$f, k_percentiles)
    post_obj <- get_posterior_obj(arr_format$training_summaries, arr_format$training_realizations, arr_format$test_summaries, decoupler)
    decoupler_margins <- post_obj$decoupler_margins
    test_question_id <- test_set$question_id |> unique()
    stopifnot(length(test_question_id) == 1)
    q_realization <- test_set$realization |> unique()
    test_decouple_values <- decoupler$f(q_realization, arr_format$test_summaries) # 1xEx\tilde{D}
    flattened <- abind::adrop(test_decouple_values, drop=1) |> flatten_matrix_row_by_row()

    cdf_values <- purrr::imap_dbl(decoupler_margins, \(margin, i) margin$cdf(flattened[i]))
    pdf_values <- purrr::imap_dbl(decoupler_margins, \(margin, i) margin$pdf(flattened[i]))

    tibble::tibble(
      likelihoods = pdf_values,
      cdf_values = cdf_values,
      test_question_id = test_question_id
    )
  }) |> purrr::list_rbind()
  stats$study_id <- study_data$study_id |> unique()
  stats
}


run_decoupler_margin_comparison <- function() {
  library(devtools)
  devtools::load_all(".")
  studies <- load_data_49(relative_dev_folder = FALSE) # |> filter_studies_few_questions(min_questions = 11)
  #studies <- filter_study_remove_ids(studies,7)
  settings <- decoupler_and_margin_estimation_settings_study()
  res <- compare_decoupler_margin_estimations(studies[1], settings)
  print("Saving to file dev/output/margin_estimation_comparison.rds")
  saveRDS(res, file = "dev/output/margin_estimation_comparison.rds")
}

if (!interactive()) {
  run_decoupler_margin_comparison()
}
