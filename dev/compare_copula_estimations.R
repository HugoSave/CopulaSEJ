source("dev/dev_utils.R")
library(parallel)

compose_copula_settings_with_posterior_fitting <- function(copula_setting) {
  purrr::partial(
    fit_and_construct_posterior_indep,
    copula_model = copula_setting$copula_model,
    vine_fit_settings = copula_setting$vine_fit_settings,
    connection_threshold = copula_setting$connection_threshold,
    q_support_restriction = NULL,
    q_support_overshoot = 0.1,
    rejection_threshold = 0.05,
    use_pseudo_obs = FALSE
  )
}


# What we will do is evaluate how well the new point is predicted in a leave one out cross validation
evalute_copula_fit <- function(study_data,
                               decoupler,
                               get_posterior_obj,
                               k_percentiles = c(5, 50, 95)) {
  fold_combinations <- create_cross_validation_sets(study_data)

  study_id <- study_data$study_id |> unique()
  stats <- purrr::pmap(fold_combinations, function(training_set, test_set) {
    arr_format <- df_format_to_array_format(
      training_set,
      test_set,
      get_three_quantiles_summarizing_function()$f,
      k_percentiles
    )
    test_question_id <- test_set$question_id |> unique()
    stopifnot(length(test_question_id) == 1)

    # post_obj <- withCallingHandlers(
    #   get_posterior_obj(
    #     arr_format$training_summaries,
    #     arr_format$training_realizations,
    #     arr_format$test_summaries,
    #     decoupler
    #   ),
    #   stan_optimization_error = function(e) {
    #     browser()
    #     message("Stan optimization error: ", conditionMessage(e), "for study_id: ", study_id, "and question_id: ", test_question_id)
    #     return(NULL)
    #   },
    #   stan_optimization_warning = function(w) {
    #     message("Stan optimization warning: ", conditionMessage(w), "for study_id: ", study_id, "and question_id: ", test_question_id)
    #     return(NULL)
    #   }
    # )

    post_obj <- tryCatch(
      get_posterior_obj(
        arr_format$training_summaries,
        arr_format$training_realizations,
        arr_format$test_summaries,
        decoupler
      ),
      stan_optimization_error = function(e) {
        message("Stan optimization error: ", conditionMessage(e), "for study_id: ", study_id, "and question_id: ", test_question_id)
        return(NULL)
      },
      stan_optimization_warning = function(w) {
        message("Stan optimization warning: ", conditionMessage(w), "for study_id: ", study_id, "and question_id: ", test_question_id)
        return(NULL)
      },
      error = function(e) { # for some reason not all 'stan_optimization_error' errors are caught by the stan_optimization_error handle
        message("Error: ", conditionMessage(e), "for study_id: ", study_id, " and question_id: ", test_question_id)
        return(NULL)
      }
    )

    if (is.null(post_obj)) {
      return(tibble::tibble_row(
          likelihood = NA,
          cum_prob = NA,
          test_question_id = test_set$question_id |> unique(),
        )
      )
    }
    copula <- post_obj$decoupler_copula
    decoupler_margins <- post_obj$decoupler_margins
    q_realization <- test_set$realization |> unique()
    test_decouple_values <- decoupler$f(q_realization, arr_format$test_summaries) # 1xEx\tilde{D}
    flattened <- abind::adrop(test_decouple_values, drop = 1) |> flatten_matrix_row_by_row()

    cdf_values <- purrr::imap(decoupler_margins, \(margin, i) margin$cdf(flattened[i])) |> do.call(what =
                                                                                                     cbind)

    tibble::tibble_row(
      likelihood = copula$density(cdf_values),
      cum_prob = copula$distribution(cdf_values),
      test_question_id = test_question_id,
    )
  }) |> purrr::list_rbind()
  stats$study_id <- study_data$study_id |> unique()
  stats
}

run_copula_estimation_comparision <- function(seed=42) {
  devtools::load_all(".")
  studies <- load_data_47(relative_dev_folder = FALSE)
  studies_few <- filter_studies_few_questions(studies, min_questions = 11)
  copula_settings <- get_copula_estimation_settings()
  decoupler = get_CDF_decoupler(overshoot = 0.1)
  decoupler = get_relative_decoupler(D_tilde=1, compose_sigmoid=TRUE, m_preprocess="median", name="Rel.Md.k(0.05)", k=0.05)


  combinations = tidyr::expand_grid(
    study = studies_few,
    setting = copula_settings
  )

  ## Use this for sequential
  # set.seed(seed)
  # res <- purrr::pmap(combinations, \(study, setting) {
  #   get_posterior_func = compose_copula_settings_with_posterior_fitting(setting)
  #   evaluation <- evalute_copula_fit(
  #     study,
  #     decoupler,
  #     get_posterior_func,
  #     k_percentiles = c(5, 50, 95)
  #   )
  #   evaluation$settings <- get_setting_name(setting)
  #   evaluation
  # }, .progress = "Evaluating copula settings")

  # Parallel if you are on linux/mac
  res <- run_mclapply_copula(
    combinations,
    decoupler,
    seed = seed,
    mc.cores = 4
  )

  res_df <- purrr::list_rbind(res)
  # save
  print("Saving results to output/copula_estimation_comparison.rds")
  saveRDS(res_df, file = "dev/output/copula_estimation_comparison.rds")
}


# Function to run mclapply with memory management
run_mclapply_copula <- function(combinations, decoupler, seed, mc.cores = 2) {

  # Check if we're on Windows (mclapply won't work)
  if (.Platform$OS.type == "windows") {
    warning("mclapply doesn't work on Windows. Will run sequentially")
  }

  cat("Starting mclapply with", mc.cores, "cores...\n")
  cat("Total combinations:", nrow(combinations), "\n")

  # Generate reproducible seeds
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, nrow(combinations))

  # Clean up environment before forking to reduce memory footprint
  # Remove any large objects you don't need
  gc(verbose = FALSE)  # Force garbage collection

  # Create a list of tasks (more memory efficient than accessing data frame repeatedly)
  tasks <- vector("list", nrow(combinations))
  for (i in seq_len(nrow(combinations))) {
    tasks[[i]] <- list(
      study = combinations$study[[i]],
      setting = combinations$setting[[i]],
      seed = seeds[i],
      index = i
    )
  }

  # Define the worker function
  process_task <- function(task) {
    # Set seed for reproducibility
    set.seed(task$seed)

    # Your original computation
    get_posterior_func <- compose_copula_settings_with_posterior_fitting(task$setting)
    evaluation <- evalute_copula_fit(
      task$study,
      decoupler,
      get_posterior_func,
      k_percentiles = c(5, 50, 95)
    )
    evaluation$settings <- get_setting_name(task$setting)

    # Optional: print progress (comment out if too noisy)
    if (task$index %% 10 == 0) {
      cat("Completed task", task$index, "\n")
    }

    return(evaluation)
  }
  # Run mclapply with conservative settings
  start_time <- Sys.time()

  res <- mclapply(
    tasks,
    process_task,
    mc.cores = mc.cores,
    mc.silent = FALSE,        # Show errors
    mc.cleanup = TRUE,        # Clean up child processes
    mc.preschedule = TRUE     # Pre-schedule tasks (can be more efficient)
  )

  end_time <- Sys.time()

  cat("mclapply completed in",
      round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")

  # Check for errors
  errors <- sapply(res, function(x) inherits(x, "try-error"))
  if (any(errors)) {
    cat("Warning:", sum(errors), "tasks produced errors.\n")
    cat("First error:", res[[which(errors)[1]]], "\n")
  }

  return(res)
}

get_setting_name <- function(setting) {
  copula_model <- setting$copula_model
  vfs <- setting$vine_fit_settings
  ct <- setting$connection_threshold

  name_parts <- c()

  # Base: copula model
  name_parts <- c(name_parts, copula_model)

  # Specifics for 'hierarchical'
  if (copula_model == "hierarchical") {
    if (!is.null(vfs$eta)) {
      name_parts <- c(name_parts, paste0("eta", vfs$eta))
    }
  }

  # Specifics for 'vine'
  if (copula_model == "vine") {
    if (!is.null(vfs$family_set)) {
      fs <- paste(vfs$family_set, collapse = "+")
      name_parts <- c(name_parts, paste0("fs=", fs))
    }
    if (!is.null(vfs$selcrit)) {
      name_parts <- c(name_parts, vfs$selcrit)
    }
    if (!is.null(vfs$psi0)) {
      name_parts <- c(name_parts, paste0("psi0", vfs$psi0))
    }
    if (!is.null(vfs$threshold)) {
      name_parts <- c(name_parts, paste0("thr", vfs$threshold))
    }
  }
  if (!is.null(ct)) {
    name_parts <- c(name_parts, paste0("ct", ct))
  }

  # Specifics for 'frank'
  if (copula_model == "frank") {
    # do nothing
    #name_parts <- c(name_parts, "frank_default")
  }

  # Combine into a single string
  return(paste(name_parts, collapse = "_"))
}

get_copula_estimation_settings <- function() {
  list(
    list(
      copula_model = "frank",
      vine_fit_settings = list(
        family_set = c("onepar", "indep"),
        selcrit = "aic"
      ),
      connection_threshold = NULL
    ),
    ## Hierarchical copula settings
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=1
      ),
      connection_threshold = 0.7,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=10
      ),
      connection_threshold = 0.7,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=50
      ),
      connection_threshold = 0.7,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=1
      ),
      connection_threshold = 0.5,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=10
      ),
      connection_threshold = 0.5,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=50
      ),
      connection_threshold = 0.5,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=1
      ),
      connection_threshold = NULL
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=10
      ),
      connection_threshold = NULL
    ),
    list(
      copula_model = "hierarchical",
      vine_fit_settings = list(
        eta=50
      ),
      connection_threshold = NULL
    ),
    # list(
    #   copula_model = "hierarchical",
    #   vine_fit_settings = list(
    #     eta=10
    #   ),
    #   connection_threshold = 0.8,
    #   connection_metric = "kendall"
    # ),
    #list(
    #  copula_model = "frank",
    #  vine_fit_settings = list(),
    #  connection_threshold = NULL
    #),
    # list(
    #   copula_model = "indep",
    #   vine_fit_settings = list(),
    #   connection_threshold = NULL
    # ),
    list(
      copula_model = "vine",
      vine_fit_settings = list(
        family_set = c("onepar", "indep"),
        selcrit = "aic"
      ),
      connection_threshold = NULL
    ),
    list(
      copula_model = "vine",
      vine_fit_settings = list(
        family_set = c("onepar", "indep"),
        selcrit = "aic"
      ),
      connection_threshold = 0.5,
      connection_metric = "kendall"
    ),
    list(
      copula_model = "vine",
      vine_fit_settings = list(
        family_set = c("onepar", "indep"),
        selcrit = "aic"
      ),
      connection_threshold = 0.7,
      connection_metric = "kendall"
    )
    #list(
    #  copula_model = "vine",
    #  vine_fit_settings = list(
    #    family_set = c("onepar", "indep"),
    #    selcrit = "aic"
    #  )
    #),
    # list(
    #   copula_model = "vine",
    #   vine_fit_settings = list(
    #     family_set = c("onepar", "indep"),
    #     selcrit = "mbicv",
    #     psi0 = 0.9,
    #     threshold = 0.8
    #   )
    # ),
    # list(
    #   copula_model = "vine",
    #   vine_fit_settings = list(
    #     family_set = c("onepar", "indep"),
    #     selcrit = "mbic",
    #     psi0 = 0.9
    #   ),
    #   connection_threshold = 0.7,
    #   connection_metric = "kendall"
    # )
    #list(
    #  copula_model = "vine",
    #  vine_fit_settings = list(
    #    family_set = c("onepar", "indep"),
    #    selcrit = "mbicv",
    #    psi0 = 0.9,
    #    threshold = 0.8
    #  ),
    #  connection_threshold = 0.05
    #),
  )
}
