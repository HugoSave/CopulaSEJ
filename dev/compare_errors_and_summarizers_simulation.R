
run_study <- function(studies, error_metric, summarizing_function, prediction_method="copula", copula_model="vine", sim_params = list()){
  print(glue::glue("Analyzing method: {prediction_method} with copula model: {copula_model} and error metric: {error_metric$name} and summarising function: {summarizing_function$name}"))

  #params <- default_simulation_params(copula_model = "vine", prediction_method = prediction_method, error_metric = error_metric, summarizing_function = summarizing_function)
  params <- rlang::inject(
    default_simulation_params(copula_model = copula_model, prediction_method = prediction_method,
                              error_metric = error_metric, summarizing_function = summarizing_function,
                              !!!sim_params
                              ))


  analysis_res <- run_analysis_per_study(studies, params)
  analysis_res$results["prediction_method"] <- glue::glue("{prediction_method}:{copula_model}:{error_metric$name}:{summarizing_function$name}")
  analysis_res
}


if (!interactive()) {
  library(devtools)
  load_all(".")
  source("dev/dev_utils.R")

  file_name <- "dev/output/data49_nov24.rds"
  numerical_col_names <- c(k_percentiles_to_colname(c(5, 50, 95)), "realization")
  studies <- readRDS(file_name)
  studies <- filter_questions_in_studies_non_negative(studies)
  #studies <- filter_studies_support_magnitude_diff(studies, numerical_col_names, max_mag_diff = 10000)
  #studies <- filter_studies_support_diff(studies, numerical_col_names, max_abs_diff = 10000)
  studies <- filter_studies_few_questions(studies, min_questions=11)
  studies <- filter_study_remove_ids(studies, study_ids=7)
  #list_mod <- change_value_in_study_list(studies, numerical_col_names, 0.001, 0.0015)$study_list
  #list_mod <- change_small_value_in_study_list(list_mod, numerical_col_names, 0.001)$study_list
  # list_mod <- change_value_in_study_list(list_mod, numerical_col_names, 0, 0.001)$study_list
  #list_mod <- filter_studies_support_magnitude_diff(list_mod, numerical_col_names, max_mag_diff = 1000)

  data_list_short <-studies[18]

  res_list <- list(
    run_study(data_list_short, get_CDF_decoupler(), get_three_quantiles_summarizing_function(), prediction_method="copula",
              copula_model="vine",
              sim_params = list(
                error_estimation_settings = list(out_of_boundary="clamp", method="kde", bw=1),
                vine_fit_settings = list(family_set = c("onepar","indep"),
                                         selcrit="mbicv", threshold=0.7),
                q_support_restriction = 'non_negative'
              )),
    run_study(data_list_short, get_CDF_decoupler(), get_three_quantiles_summarizing_function(), prediction_method="copula",
              copula_model = "indep",
              sim_params = list(
                error_estimation_settings = list(out_of_boundary="clamp", method="kde", bw=1),
                q_support_restriction = 'non_negative'
              )),
    run_study(data_list_short, get_linear_decoupler(), get_median_summarizing_function(), prediction_method="copula",
              copula_model = "vine",
              sim_params = list(
                q_support_restriction = 'non_negative',
                vine_fit_settings = list(family_set = c("onepar","indep"),
                                          selcrit="mbicv", threshold = 0.7)
              )),
    run_study(data_list_short, get_linear_decoupler(), get_median_summarizing_function(), prediction_method="copula",
              copula_model = "indep",
              sim_params = list(
                q_support_restriction = 'non_negative'
              )),

    run_study(data_list_short, get_linear_error_metric(), get_median_summarizing_function(), prediction_method="perfect_expert",
              copula_model = "frank",
              sim_params = list(
                q_support_restriction = 'non_negative'
              ))
    # run_study(data_list_short, get_linear_error_metric(), get_median_summarizing_function(), prediction_method="copula",
    #           copula_model = "vine",
    #           sim_params = list(
    #             q_support_restriction = 'non_negative',
    #             vine_fit_settings = list(family_set = c("all"),
    #                                      selcrit="mbicv", threshold = 0.7)
    #           ))
  )
  # vine_fit_settings = list(threshold = 0.8),
  res_combined <- res_list |> purrr::transpose()
  post_df = res_combined$results |> dplyr::bind_rows()
  warnings = res_combined$warning

  # error_metrics = list(get_sigmoid_relative_error_metric(), get_sigmoid_linear_error_metric())
  # #error_metrics = list(get_ratio_error_metric(), get_linear_error_metric())
  # summarizing_functions = list(get_median_summarizing_function())
  # for (error_metric in error_metrics) {
  #   for (summarizing_function in summarizing_functions) {
  #     res <- run_study(data_list_short, error_metric, summarizing_function)
  #     results <- append(results, list(res$results))
  #     warnings <- append(warnings, list(res$warnings))
  #   }
  # }

  print("Sampling predictions...")

  #preds <- results_df$posterior |> purrr::map_dbl(\(post) mean(sample_log_unnormalized_density(post$logDM, post$support, 1000)))
  list_of_metrics <- purrr::map2(post_df$posterior, post_df$realization, \(post, realization) {
    if (is.null(post)) {
      return(performance_metrics_list())
    } else {
      return(posterior_performance_metrics(post$logDM, post$support, realization, num_samples=1000))
    }

  })

  metric_df <- list_of_metrics |> purrr::list_transpose() |> tibble::as_tibble()
  results <- dplyr::bind_cols(post_df, metric_df)
  results$prediction <- results$median
  results$error = results$prediction - results$realization
  results$rel_error = results$error / results$realization
  results$support = results$posterior |> purrr::map(\(x) x$support)
  results$posterior = NULL  # we drop this column

  saveRDS(list(results_df=results, warnings=warnings), "dev/output/compare_errors_and_summarizers_simulation.rds")

  print("Results saved to dev/output/compare_errors_and_summarizers_simulation.rds")
}


