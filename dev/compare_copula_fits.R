source("dev/dev_utils.R")

get_posterior_obj <- function(training, realization, test, decoupler) {
  fit_and_construct_posterior_indep(
    training, realization, test, decoupler,
    copula_model = "vine",
    vine_fit_settings = list(
      family_set = c("onepar","indep"),
      selcrit="mbicv", psi0=0.9,
      threshold=0.0
    ),
    error_estimation_settings=list(
      method="beta_hierarchical",
      k=1
    ),
    q_support_restriction=NULL,
    q_support_overshoot=0.1,
    rejection_threshold=NULL,
    rejection_min_experts=1,
    rejection_test="distance_correlation",
    connection_threshold=NULL
  )
}

compare_decoupler_margin_estimations <- function(studies, settings) {
  # We don't need to filter for questions with more than 10 questions because we are not doing any copula fitting.
  combinations <- tidyr::expand_grid(study=studies, settings)
  res <- purrr::pmap(combinations[577,], \(study, error_estimation_settings, decoupler) {
    get_obj_func <- purrr::partial(fit_and_construct_posterior_indep,
                   error_estimation_settings = error_estimation_settings,
                   copula_model = "indep",
                   q_support_restriction=NULL,
                   q_support_overshoot=0.1,
                   rejection_threshold=NULL)
    res <- evalute_marginal_fit(study, decoupler, get_obj_func)
    res$settings <- decoupler_margin_estimation_settings_to_shortname(error_estimation_settings)
    res$decoupler <- decoupler$short_name
    res
  })
  res <- res |> purrr::list_rbind()
  res
}

decoupler_margin_estimation_settings_to_shortname <- function(settings) {
  # settings =list(
  #  method="beta_hierarchical",
  #  beta_mean=0.5,
  #  beta_var=1/12,
  #  prior_var=1
  #)
  settings$method <- switch(
    settings$method,
    "beta_hierarchical" = "betaH",
    "beta" = "beta",
  )

  glue::glue(
    "{settings$method}:{round(settings$prior_var, digits=4)}"
  )

}

decoupler_and_margin_estimation_settings_study <- function() {
  decouplers <- list(
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=1),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.5),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.1),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.05),
    get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k=0.01),
    get_CDF_decoupler()
  )
  error_estimation_settings_list <- list(
    list(
      method="beta_hierarchical",
      prior_var=1
    ),
    list(
      method="beta_hierarchical",
      prior_var=0.1
    ),
    list(
      method="beta_hierarchical",
      prior_var=0.01
    ),
    list(
      method="beta_hierarchical",
      prior_var=0.001
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

evalute_copula_fit <- function(study_data, decoupler, get_posterior_obj, k_percentiles=c(5,50,95)) {

    fold_combinations <- create_cross_validation_sets(study_data)

    stats <- list(
      prediction = vector(mode = "numeric"),
      realization = vector(mode = "numeric"),
      test_question = vector(mode = "numeric"),
      posterior = list(),
      accepted_experts = vector(mode = "numeric"),
      experts = vector(mode = "numeric")
    )
    warnings <- c()
    single_expert_warning <- FALSE
    below_alpha_warning <- FALSE
    total_nr_experts <- length(unique(study_data$expert_id))

    set.seed(1234)
    stats <- purrr::pmap(fold_combinations, function(training_set, test_set) {
      arr_format <- df_format_to_array_format(training_set, test_set, get_three_quantiles_summarizing_function()$f, k_percentiles)
      post_obj <- get_posterior_obj(arr_format$training_summaries, arr_format$training_realizations, arr_format$test_summaries, decoupler)
      copula <- post_obj$decoupler_copula
      decoupler_margins <- post_obj$decoupler_margins
      stopifnot((test_set$question_id |> unique() |> length()) == 1)
      q_realization <- test_set$realization |> unique()
      test_decouple_values <- decoupler$f(q_realization, arr_format$test_summaries) # 1xEx\tilde{D}
      flattened <- abind::adrop(test_decouple_values, drop=1) |> flatten_matrix_row_by_row()

      cdf_values <- purrr::imap(decoupler_margins, \(margin, i) margin$cdf(flattened[i])) |> do.call(what=cbind)
      #pdf_values <- purrr::imap(decoupler_margins, \(margin, i) margin$pdf(flattened[i])) |> do.call(what=cbind)
      #cdf_values_training <- purrr::imap(decoupler_margins, \(margin, i) margin$cdf(post_obj$flattened_errors[,i])) |> do.call(what=cbind)

      tibble::tibble_row(
        likelihood = copula$density(cdf_values),
        cum_prob = copula$distribution(cdf_values)
      )
    }) |> purrr::list_rbind()
}

copula_fits <- list(
list(copula_model = "vine",
         vine_fit_settings = list(
           family_set = c("onepar","indep"),
           selcrit="mbicv", psi0=0.9,
           threshold=0.0
         ),
     connection_threshold=NULL
     ),
list(copula_model = "vine",
     vine_fit_settings = list(
       family_set = c("onepar","indep"),
       selcrit="mbicv", psi0=0.9,
       threshold=0.8
     ),
     connection_threshold=NULL
     ),
list(copula_model = "vine",
     vine_fit_settings = list(
       family_set = c("onepar","indep"),
       selcrit="mbicv", psi0=0.9,
       threshold=0.8
     ),
     connection_threshold=0.05)
)

if (!interactive()) {
  library(devtools)
  devtools::load_all(".")
  studies <- load_data_49(relative_dev_folder = FALSE) # |> filter_studies_few_questions(min_questions = 11)
  settings <- decoupler_and_margin_estimation_settings_study()
  res <- compare_decoupler_margin_estimations(studies, settings)
  print("Saving to file dev/output/margin_estimation_comparison.rds")
  saveRDS(res, file = "dev/output/margin_estimation_comparison.rds")
}
