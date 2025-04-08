
run_study <- function(studies, error_metric, summarizing_function, prediction_methd="copula", sim_params = list()){
  print(paste("Analyzing error metric:", error_metric$name, ", and summarizing function:", summarizing_function$name, ", and prediction method:", prediction_methd))

  #params <- default_simulation_params(copula_model = "vine", prediction_method = prediction_methd, error_metric = error_metric, summarizing_function = summarizing_function)
  params <- rlang::inject(
    default_simulation_params(copula_model = "vine", prediction_method = prediction_methd,
                              error_metric = error_metric, summarizing_function = summarizing_function,
                              !!!sim_params
                              ))

  analysis_res <- run_analysis_per_study(studies, params)
  analysis_res$results["prediction_method"] <- glue::glue("copula={error_metric$name}:{summarizing_function$name} summarizer:prediction method={prediction_methd}")
  analysis_res
}


if (!interactive()) {
  library(devtools)
  load_all(".")
  source("dev/dev_utils.R")

  file_name <- "dev/output/data49_nov24.rds"
  studies <- readRDS(file_name)
  studies <- filter_questions_in_studies_non_negative(studies)
  studies <- filter_studies_few_questions(studies, min_questions=11)
  studies <- filter_study_remove_ids(studies, study_ids=7)
  numerical_col_names <- c(k_percentiles_to_colname(c(5, 50, 95)), "realization")
  list_mod <- change_value_in_study_list(studies, numerical_col_names, 0.001, 0.0015)$study_list
  list_mod <- change_value_in_study_list(list_mod, numerical_col_names, 0, 0.001)$study_list


  data_list_short <- list_mod
  results <- list()
  warnings <- list()

  res <- run_study(data_list_short, get_CDF_decoupler(), get_three_quantiles_summarizing_function(),
                   sim_params = list(
                     error_estimation_settings = list(out_of_boundary="clamp", method="kde", bw=1)
                   ))
  results <- append(results, list(res$results))
  warnings <- append(warnings, list(res$warnings))

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

  res <- run_study(data_list_short, get_linear_error_metric(), get_median_summarizing_function(), "copula_assumption")
  results <- append(results, list(res$results))
  warnings <- append(warnings, list(res$warnings))

  results_df <- dplyr::bind_rows(results)

  print("Sampling predictions...")

  preds <- results_df$posterior |> purrr::map_dbl(\(post) mean(sample_log_unnormalized_density(post$logDM, post$support, 1000)))
  results_df$prediction = preds
  results_df$error = results_df$prediction - results_df$realization
  results_df$rel_error = results_df$error / results_df$realization

  saveRDS(list(results_df=results_df, warnings=warnings), "dev/output/compare_errors_and_summarizers_simulation.rds")

  print("Results saved to dev/output/compare_errors_and_summarizers_simulation.rds")
}


