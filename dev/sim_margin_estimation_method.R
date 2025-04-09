library(CopulaSEJ)
source("dev/dev_utils.R")


produce_posteriors <- function(){
  studies <- load_data_49_filtered()
  studies <- studies[1:2]

  error_metrics = list(get_ratio_error_metric(), get_linear_error_metric())
  summarizing_functions = list(get_three_quantiles_summarizing_function(), get_median_summarizing_function())

  params <- default_simulation_params(copula_model = "vine",
                                      prediction_method = "copula",
                                      vine_fit_settings = list(family_set=c("one_param", "indep",
                                                                            selcrit="mbicv",
                                                                            psi0=0.6)))

  results <- list()
  warnings <- list()
  for (error_metric in error_metrics) {
    for (summarizing_function in summarizing_functions) {
      print(paste("Analyzing error metric:", error_metric$name, ", and summarizing function:", summarizing_function$name))

      params <- default_simulation_params(copula_model = "vine", prediction_method = "copula", error_metric = error_metric, summarizing_function = summarizing_function)

      analysis_res <- run_analysis_per_study(studies, params)
      analysis_res$results["prediction_method"] <- glue::glue("perfect_expert:{error_metric$name}:{summarizing_function$name} summarizer")
      results <- append(results, list(analysis_res$results))
      warnings <- append(warnings, list(analysis_res$warnings))
    }
  }

  results_df <- dplyr::bind_rows(results)
  settings <- list(
    params,
    error_metrics,
    summarizing_functions,
    "load_data_49_filtered"
  )
  saveRDS(list(results=results, settings=settings), "dev/output/copula_vine_metric_comparison.rds")
}


produce_posterior_estimates <- function(res_df) {
  row_list <- res_df |> split(seq_len(nrow(res_df)))
  results <- purrr::map(row_list, function(row) {
    log_post <- row$posterior[[1]]$logDM
    support <- row$posterior[[1]]$support
    realization <- row$realization
    posterior_performance_metrics(log_post, support, realization, num_samples=500)
  })
}
