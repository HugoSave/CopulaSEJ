
copula_shortname <- function(copula_family) {
  # Create a short name for the copula model
  if (copula_family == "vine") {
    return("Vi")
  } else if (copula_family == "indep") {
    return("In")
  } else if (copula_family == "frank") {
    return("Fr")
  } else if (copula_family == "hierarchical") {
    return("Hi")
  } else {
    return(copula_family)
  }
}

rejection_shortname <- function(rejection_test, rejection_threshold, rejection_min_experts) {
  if (is.null(rejection_threshold)) {
    return("NoR")
  }
  # Create a short name for the rejection test
  if (rejection_test == "distance_correlation") {
    prefix <- "DC"
  } else if (rejection_test == "classical_calibration") {
    prefix <- "CS"
  } else if (rejection_test == "kruskal") {
    prefix <- "K"
  } else {
    return(rejection_test)
  }
  return(glue::glue("{prefix}({rejection_threshold},{rejection_min_experts})"))
}

prediction_method_shortname <- function(prediction_method) {
  # Create a short name for the prediction method
  if (prediction_method == "copula") {
    return("Co")
  } else if (prediction_method == "perfect_expert") {
    return("PE")
  } else if (prediction_method == "density_product") {
    return("DP")
  }
  else if (prediction_method == "median_average") {
    return("MA")
  } else if (prediction_method == "equal_weights") {
    return("EW")
  } else if (prediction_method == "classical_global_opt") {
    return("CGO")
  } else if (prediction_method == "uniform") {
    return("Unif")
  } else {
    return(prediction_method)
  }
}

error_metric_and_summ_shortname <- function(error_metric, summarizing_function) {
  # Create a short name for the error metric and summarizing function
  if (error_metric$name == "CDF") {
    stopifnot(summarizing_function$name == "three_quantiles")
    return("CDF")
  } else {
    return(glue::glue("{error_metric$name}_{summarizing_function$name}"))
  }
}

connection_threshold_shortname <- function(connection_threshold) {
  # Create a short name for the connection threshold
  if (is.null(connection_threshold)) {
    return("NoCT")
  } else {
    return(glue::glue("CT({connection_threshold})"))
  }
}

beta_std_shortname <- function(error_estimation_settings) {
  if (is.null(error_estimation_settings) || is.null(error_estimation_settings$prior_std)) {
    return("")
  } else {
    return(glue::glue("BS({error_estimation_settings$prior_std})"))
  }
}


parameter_shortname <- function(sim_params, sep="") {
  method = sim_params$prediction_method

  # For copula method, use copula model name instead of prediction method
  if (method == "copula") {
    # Build name components for copula method
    components <- c(
      copula_shortname(sim_params$copula_model),
      sim_params$error_metric$short_name,
      sim_params$summarizing_function$short_name,
      rejection_shortname(sim_params$rejection_test,
                          sim_params$rejection_threshold,
                          sim_params$rejection_min_experts),
      connection_threshold_shortname(sim_params$connection_threshold),
      beta_std_shortname(sim_params$error_estimation_settings)
    )
    # Remove empty components
    components <- components[components != ""]
    return(paste(components, collapse = sep))
  }
  # For non-parameter methods (no copula involved)
  else if (method %in% c("density_product", "uniform", "median_average", "equal_weights", "classical_global_opt")) {
    return(prediction_method_shortname(method))
  }
  # For other methods that might use parameters
  else {
    components <- c(
      prediction_method_shortname(sim_params$prediction_method),
      copula_shortname(sim_params$copula_model),
      sim_params$error_metric$short_name,
      sim_params$summarizing_function$short_name,
      rejection_shortname(sim_params$rejection_test,
                          sim_params$rejection_threshold,
                          sim_params$rejection_min_experts),
      connection_threshold_shortname(sim_params$connection_threshold),
      beta_std_shortname(sim_params$error_estimation_settings)
    )
    # Remove empty components
    components <- components[components != ""]
    return(paste(components, collapse = sep))
  }
}

get_latest_sim_nr <- function() {
  # Get the latest simulation number for the given simulation parameters
  date_tag <- format(Sys.time(), "%Y-%m-%d")
  folder <- "dev/output/"
  file_list <- list.files(folder, pattern = "sim_.*")
  # Extract the simulation numbers from the file names on current date
  sim_numbers <- grep(paste0(date_tag,"_sim_([0-9]+)"), file_list, perl = TRUE, value=TRUE)

  if (length(sim_numbers) == 0) {
    return(NULL)
  } else {
    # Extract the simulation numbers from the file names
    sim_numbers <-as.numeric(sim_numbers)
    # Return the maximum simulation number
    return(max(sim_numbers))
  }
}

get_simulation_group_files <- function(group="tmp", date=NULL, rel_dev=FALSE) {
  # Get the simulation group based on the provided group and date
  if (is.null(date)) {
    date <- format(Sys.time(), "%Y-%m-%d")
  }
  if (rel_dev) {
    folder <- "output/"
  } else {
    folder <- "dev/output/"
  }
  file_list <- list.files(folder, pattern = paste0(date, "_", group, "_"), full.names=TRUE)
  # Extract the simulation numbers from the file names
  return(file_list)
}

load_files <- function(files) {
  files |> purrr::map(readRDS)
}

combine_simulation_results <- function(list_results) {
  # Combine the simulation results from the provided files
  combined_results <- purrr::map_dfr(list_results, \(x) x$results)
  combined_warnings <- purrr::map(list_results, \(x) x$warnings)
  combined_warnings <- purrr::reduce(combined_warnings, \(x,y) c(x,y))
  return(list(results=combined_results, warnings=combined_warnings))
}

create_file_name <- function(sim_params, sim_group="tmp", sim_nr=NULL) {
  # Create a file name based on the simulation parameters
  date_tag <- format(Sys.time(), "%Y-%m-%d")
  if (is.null(sim_nr)) {
    latest_nr <- get_latest_sim_nr()
    sim_nr <- if (is.null(latest_nr)) { 1 } else { latest_nr + 1 }
  }
  file_name <- paste0("dev/output/",
                      date_tag,
                      "_",
                      sim_group,
                      "_",
                      parameter_shortname(sim_params, "_"),
                      ".rds")

  return(file_name)
}

run_study_find_posterior <- function(studies, params, sim_group="tmp"){
  print(glue::glue("Analyzing method: {params$prediction_method} with copula model: {params$copula_model} and error metric: {params$error_metric$name} and summarising function: {params$summarizing_function$name}"))

  analysis_res <- run_analysis_per_study(studies, params)
  analysis_res$results["prediction_method"] <- parameter_shortname(params)
  #file_name <- create_file_name(params, sim_group=sim_group)
  analysis_res$params <- params
  #print(glue::glue("Saving results to {file_name}"))
  #saveRDS(analysis_res, file_name)
  #analysis_res$file_name <- file_name
  analysis_res
}

sample_df_parallel <- function(df, num_samples=5000, n_cores=4) {

  # Extract only what you need first (avoids copying the whole df)
  post_list <- df$posterior
  realization_list <- df$realization
  sample_prior_list <- df$sample_prior

  # Create a minimal task list
  tasks <- lapply(1:length(post_list), function(i) {
    list(
      post = post_list[[i]],
      realization = realization_list[[i]],
      sample_prior = sample_prior_list[[i]]
    )
  })

  list_of_metrics <- mclapply(tasks, function(task) {
    if (!is.list(task$post)) {
      return(performance_metrics_list())
    } else {
      return(posterior_performance_metrics(
        task$post$logDM, task$post$support, task$realization,
        num_samples = num_samples, mean_value = task$post$mean,
        median_value = task$post$median, sample_prior = task$sample_prior
      ))
    }
  }, mc.cores = n_cores)
  return(list_of_metrics)

}

sample_and_add_metrics <- function(analys_res, run_parallel=TRUE, num_samples=5000, n_cores=4) {

  df <- analys_res$results
  if (run_parallel) {
    list_of_metrics <- sample_df_parallel(df, num_samples=num_samples, n_cores=n_cores)
  } else {
    # Fallback to sequential processing
    list_of_metrics <- purrr::pmap(list(df$posterior, df$realization, df$sample_prior), \(post, realization, sample_prior) {
      if (!is.list(post)) {
        return(performance_metrics_list()) # defaults to NA for all other metrics
      } else {
        return(posterior_performance_metrics(post$logDM, post$support, realization, num_samples=num_samples, mean_value=post$mean, median_value=post$median, sample_prior = sample_prior))
      }
    })
  }

  metric_df <- list_of_metrics |> purrr::list_transpose() |> tibble::as_tibble()
  results <- dplyr::bind_cols(df, metric_df)
  results$prediction <- dplyr::if_else(is.finite(results$prediction), results$prediction, results$median)
  # results$prediction <- results$median
  results$error = results$prediction - results$realization
  results$rel_error = results$error / results$realization
  results$support = results$posterior |> purrr::map(\(x) if(is.list(x)) x$support else NA)
  results
}

push_list <- function(list, item) {
  list[[length(list) + 1]] <- item
  return(list)
}

run_benchmarking_methods <- function() {
  devtools::load_all(".")
  source("dev/dev_utils.R")

  file_name <- "dev/output/data49_nov24.rds"
  studies <- load_data_49(relative_dev_folder = FALSE)
  #studies <- filter_studies_few_questions(studies, min_questions=11)
  studies <- filter_study_remove_ids(studies, study_ids=7)

  data_list_short <- studies
  param_list <- list()

  param_list <- push_list(param_list,
                          default_simulation_params(
                            prediction_method = "density_product",
                            q_support_overshoot=0.1,
                            summarizing_function = get_three_quantiles_summarizing_function(),
                            q_support_restriction = NULL,
                          ))

  param_list <- push_list(param_list,
                          default_simulation_params(
                            prediction_method = "median_average",
                            summarizing_function = get_three_quantiles_summarizing_function(),
                          ))

  param_list <- push_list(param_list,
                          default_simulation_params(
                            prediction_method = "equal_weights"
                          ))

  param_list <- push_list(param_list,
                          default_simulation_params(
                            prediction_method = "classical_global_opt"
                          ))


  param_list <- push_list(param_list,
                          default_simulation_params(
                            prediction_method = "uniform"
                          ))

  result_list <- purrr::map(param_list, \(x) {
    analys_res <- run_study_find_posterior(data_list_short, x, "benchmark")
    results_with_metrics <- sample_and_add_metrics(analys_res, run_parallel = TRUE)
    results_with_metrics$posterior <- NULL
    file_name <- create_file_name(x, "benchmark")
    print(glue::glue("Saving results to {file_name}"))
    saveRDS(results_with_metrics, file_name)
    analys_res$results <- results_with_metrics
    analys_res
  })

  res_combined <- result_list |> combine_simulation_results()

  saveRDS(res_combined, "dev/output/benchmarking_methods_performance.rds")

  print("Results saved to dev/output/benchmarking_methods_performance.rds")
}

run_performance_test <- function() {
  library(devtools)
  devtools::load_all(".")
  source("dev/dev_utils.R")

  studies <- load_data_47(relative_dev_folder = FALSE)
  #studies <- filter_study_remove_ids(studies, study_ids=7)

  data_list_short <-studies
  param_list <- list()

  decoupler_estimation_setting_list <- list(
    list(
      method = "beta_prior"
    ),
    list(
      method = "beta_MAP",
      prior_std = 0.1
    ),
    list(
      method = "beta_MAP",
      prior_std = 0.5
    ),
    list(
      method = "beta_MAP",
      prior_std = 1
    ),
    list(
      method = "beta_MAP",
      prior_std = 1.5
    ),
    list(
      method = "beta_MAP"
    )
  )
  copula_models <- list("indep", "hierarchical")
  decouplers <- list(
    get_CDF_decoupler(),
    get_relative_decoupler(D_tilde=1, compose_sigmoid=TRUE, m_preprocess="median", k=0.05)
  )

  connection_thresholds <- list(0.7)

  settings <- expand.grid(
    copula_model = copula_models,
    decoupler = decouplers,
    decoupler_estimation_settings = decoupler_estimation_setting_list,
    ct = connection_thresholds
  ) |> as_tibble()

  for (i in seq_len(nrow(settings))) {
    copula_model <- settings$copula_model[[i]]
    decoupler <- settings$decoupler[[i]]
    decoupler_estimation_settings <- settings$decoupler_estimation_settings[[i]]
    connection_threshold <- settings$ct[[i]]

    param_list <- push_list(param_list,
                            default_simulation_params(
                              prediction_method = "copula",
                              copula_model = copula_model,
                              error_metric = decoupler,
                              summarizing_function = get_three_quantiles_summarizing_function(),
                              q_support_restriction = NULL,
                              q_support_overshoot = 0.1,
                              rejection_threshold = 0.05,
                              rejection_min_experts = 1,
                              rejection_test = "distance_correlation",
                              vine_fit_settings = list(
                                eta=10,
                                recover_numerical_failure = TRUE
                              ),
                              connection_threshold = connection_threshold,
                              error_estimation_settings = decoupler_estimation_settings
                            ))
  }

  result_list <- purrr::map(param_list, \(x) {
    analys_res <- run_study_find_posterior(data_list_short, x, "main")
    results_with_metrics <- sample_and_add_metrics(analys_res, run_parallel = TRUE, num_samples=5000, n_cores=8)
    results_with_metrics$posterior <- NULL
    file_name <- create_file_name(x, "main")
    print(glue::glue("Saving results to {file_name}"))
    saveRDS(results_with_metrics, file_name)
    analys_res$results <- results_with_metrics
    analys_res
  }, .progress="Parameter choice")
  #set.seed(42, "L'Ecuyer")
  #result_list <- mclapply(param_list[1:2], function(x) {
  #  analys_res <- run_study_find_posterior(data_list_short, x, "main")
  #  results_with_metrics <- sample_and_add_metrics(analys_res)
  #  results_with_metrics$posterior <- NULL
  #  file_name <- create_file_name(x, "main")
  #  print(glue::glue("Saving results to {file_name}"))
  #  saveRDS(results_with_metrics, file_name)
  #  analys_res$results <- results_with_metrics
  #  analys_res
  #}, mc.cores = 2)

  res_combined <- result_list |> combine_simulation_results()

  file_name = "dev/output/compare_performance_simulation.rds"
  saveRDS(res_combined, file_name)

  print(glue::glue("Results saved to {file_name}"))

}
