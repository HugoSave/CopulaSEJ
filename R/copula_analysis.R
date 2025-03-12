# Read formatted data
# #######################################################
# read formatted data




add_rel_error_data <- function(data) {
  # Get the relative error of an experts median assessment relative to the realization.
  data |> dplyr::mutate(rel_error = (`50th percentile` - realization) / realization)
}

add_error_data <- function(data) {
  # Get the error of an experts median assessment relative to the realization.
  data |> dplyr::mutate(error = (`50th percentile` - realization))
  # arguably this is how it is defined in the paper Jouini and Clemen but
  # it does not seem to work well
  #data |> mutate(error = (realization - `50th percentile`))
}

calculate_error <- function(data, err_metric="error") {
  # Calculate the error of the median assessment relative to the realization.
  if (err_metric == "error") {
    data$`50th percentile` - data$realization
  } else if (err_metric == "rel_error") {
    (data$`50th percentile` - data$realization)/data$realization
  } else {
    stop("Unknown error metric")
  }
}

create_error_observations <- function(data, err_metric="error") {
  # if we have E experts we get a E dimensional return
  # if we have N questions then we have N observations
  data$error <- calculate_error(data, err_metric)
  t(data |> split(data$question_id) |> map(~.$error) |> simplify2array())
}

widen_error_per_expert <- function(study_data, err_name="error") {
  # Flatten the data per expert_id
  df <- study_data |> dplyr::select(dplyr::all_of(err_name), expert_id, question_id)
  df |> dplyr::pivot_wider(names_from = expert_id, names_prefix = "Expert ", values_from = dplyr::all_of(err_name))
}

plot_error_between_experts <- function(study_data, err_name="rel_error") {
  # So if we have E experts we have E * (E-1) / 2 pairs of experts
  # Flatten the data per expert_id
  wide_df <- widen_error_per_expert(study_data, err_name)
  pairs(wide_df, main = paste("Error between experts using column name:", err_name))
}

# error_obs is nx(d*E) matrix (n questions, d * E errors)
find_copula <- function(error_obs, copula_model="joe", fit_method="ml") {
  # Find
  error_length = ncol(error_obs)
  lower=NULL
  upper=NULL

  pseudo_obs <- rvinecopulib::pseudo_obs(error_obs)

  # fit a copula to the data
  if (copula_model == "joe") {
    copula_model <- copula::joeCopula(dim=error_length)
    lower = 1
  } else if (copula_model == "indep") {
    return(copula::indepCopula(dim=error_length))
  } else if (copula_model == "frank") {
    copula_model <- copula::frankCopula(dim=error_length)
    lower = 0
  } else if (copula_model == "clayton") {
    copula_model <- copula::claytonCopula(dim=error_length)
  } else if (copula_model == "gumbel") {
    copula_model <- copula::gumbelCopula(dim=error_length)
  } else if (copula_model == "t") {
    copula_model <- copula::tCopula(dim=error_length)
  } else if(copula_model == "normal"){
    copula_model <- copula::normalCopula(dim=error_length, dispstr = "un")
  } else if(copula_model == "vine") {
    return(rvinecopulib::vinecop(pseudo_obs, family_set = "onepar", cores=2))
    #

    # res <- try(rvinecopulib::vinecop(pseudo_obs), silent = TRUE)
    # if (inherits(res, "try-error")) {
    #   if (attr(res, "condition")$message == "copula has not been fitted from data or its parameters have been modified manually") {
    #     warning("vine not fitted using familiy_set 'all'. Trying 'onepar' instead.")
    #     res <- rvinecopulib::vinecop(pseudo_obs, family_set = "onepar")
    #     return(res)
    #   }
    #   stop("Error in vinecopulib::vinecop")
    # } else {
    #   return(rvinecopulib::vinecop(pseudo_obs))
    # }
  }
  else {
    stop("Unknown copula model")
  }

  # check if copula is normal or tcopula class
  if (is(copula_model, "acopula")) {
    # mpl does not work well for higher dimensions and gives the same output as ml
    # regardless. Only the estimation of the standard error is different.
    kendall_guess <- mean(copula::P2p(copula::corKendall(pseudo_obs)))
    if (kendall_guess < 0) {
      kendall_guess = 0
    }
    alpha_start <- copula::iTau(copula_model, kendall_guess)

    # optim.control=list(factr=1e8), makes so that the tolerance is about 1e-7.
    fitted_params <- copula::fitCopula(copula_model, pseudo_obs, method = "mpl",
                               start=alpha_start,
                               estimate.variance=FALSE,
                               lower=lower,
                               upper=upper,
                               optim.control=list(factr=1e8),
                               optim.method="L-BFGS-B")
  } else{
    fitted_vine_cop <- rvinecopulib::vinecop(pseudo_obs)
    fiited_params <- copula::fitCopula(copula_model, pseudo_obs, method = "mpl",
                                       optim.control=list(factr=1e8),
                                       optim.method="L-BFGS-B")
  }
  #fitted_param <- optim(alpha_start, loglikCopula, lower=lower, upper=upper,
  #      method = "L-BFGS-B", copula = copula_model, u = pseudo_obs)

  # Sometimes a specific starting value is needed for convergence, sometimes not.
  # fitted_params <- fitCopula(copula_model, psuedo_obs, method = "ml", start=1)

  # copula_model@parameters <- coef(fitted_params)
  #tau(copula_model)
  #print(cor(wide_df, method = "kendall"))
  fitted_params@copula
}

# Widen percentile columns
widen_percentiles <- function(data) {
  # convert columns like `xth percentile` to a percentile column
  data |> dplyr::pivot_longer(cols = dplyr::contains("percentile"), names_to = "k", values_to = "k_percentile") |>
    dplyr::mutate(k = readr::parse_number(k))
}

# Returns a data frame with the L, U, L* and U* for the assessmenet given. Must be usedd for a single study and question.
calculate_assessment_support <- function(assessments, k_percentiles=c(5,50,95), overshoot=0.01) {
  # Get the support of the assessments
  min_k = min(k_percentiles)
  min_k_colname = k_percentiles_to_colname(min_k)
  max_k = max(k_percentiles)
  max_k_colname = k_percentiles_to_colname(max_k)
  support_df <- assessments |> summarise(
    across(all_of(min_k_colname), min, .names="L"), across(all_of(max_k_colname), max, .names="U"),
    .groups="drop") |> mutate(L_star = L - (U - L) * overshoot, U_star = U + (U - L) * overshoot)
  # single row df to vector
  support_df |> as.list()
}


add_0_and_100_percentiles <- function(data, k_percentiles=c(5,50,95), overshoot=0.1) {
  support_df <- calculate_assessment_support(data, k_percentiles, overshoot)

  group_by_vars = c("expert_id", "study_id")
  data |> dplyr::group_by(dplyr::across(any_of(group_by_vars))) |>
    dplyr::mutate(`0th percentile` = support_df$L_star, `100th percentile` = support_df$U_star) |> dplyr::ungroup()
}



#' Constructs the error distribution from JC assumption. I actually don't think
#' this is currently needed however.
#'
#' @param distribution A single distribution list for one expert
#' @param error_metric A error metric class object
#' @param expert_m A vector of summarized properties for a single expert. Length d.
#'
#' @export
#'
construct_error_distribution_JC_assumption <- function(distribution, error_metric, expert_m) {
  error_increasing = error_metric$f_increasing_q(expert_m)
  force(distribution)
  force(error_metric)
  d = length(expert_m)

  error_distributions <- map(seq_along(expert_m), \(m_i) {
    force(m_i)
    F_e <- function(e) {
      distribution$cdf(error_metric$f_inverse_q(m_i,e))
    }
    f_e <- function(e) {
      distribution$pdf(error_metric$f_inverse_q(m_i,e)) * abs(error_metric$f_prime_inverse_q(m_i,e))
    }
    # calculate new support. Under the assumption that the error metric is
    # monotonic we can calculate it like below
    support1 <- error_metric$f(m_i, distribution$support[1])[1,1]
    support2 <- error_metric$f(m_i, distribution$support[2])[1,1]
    error_support = c(min(support1, support2), max(support1, support2))

    list(cdf = F_e, pdf = f_e, support = error_support)
  })
  error_distributions
}

test_construct_error_distribution_JC_assumption <- function() {
  # test the find median function
  expert_m <- c(-1, 2, 3)
  distribution <- list(cdf = function(x) {pnorm(x, mean = 0, sd = 1)}, pdf = function(x) {dnorm(x, mean = 1, sd = 1)}, support = c(-Inf, Inf))
  error_metric <- get_linear_error_metric()
  error_dists <- construct_error_distribution_JC_assumption(distribution, error_metric, expert_m)
  # plot error dist
  x = seq(-5, 5, by = 0.01)
  y = map(error_dists, \(dist) {
    dist$pdf(x)
  })
  walk(y, \(y_i) {
    plot(x, y_i)
    par(new=TRUE)
  })

}

#' Title
#'
#' @param distributions list of n number of distributions
#' @param error_metric Error metric of S3 class "error_metric"
#' @param assessments nxd matrix or n long vector frame containing the assessments
#'
#' @returns List of nxd error distribution. If assessmenets is a nxd matrix then
#' the output is a nested list with outer length n and inner length d. If
#' assessments is a n long vector then output is a n long list. Each error
#' distribution is a list with fields: cdf, pdf and support.
#' @export
#'
construct_error_distributions_JC_assumption <- function(distributions, error_metric, assessments) {
  is_increasing_vals <- error_metric$f_increasing_q(assessments)

  error_distributions <- imap(distributions, \(distribution, i) {

    is_increasing_vals[i,]

  })

  #
}

find_median_of_fun <- function(fun, support) {
  N = 1000
  param_mesh <- seq(support[1], support[2], length.out = N)
  # pick the guess that maximizes the posterior
  y <- fun(param_mesh)
  # calculate the cum
  cum_area <- cumsum(y) * ((support[2] - support[1])/ (N-1))
  # find the median
  median_guess <- param_mesh[which.min(abs(cum_area - 0.5))]
  median_guess
}


test_find_median_of_fun <- function() {
  # test the find median function
  fun <- function(x) {
    dnorm(x, mean = 1, sd = 1)
  }
  find_median_of_fun(fun, c(-5, 5))
}


k_percentiles_to_colname <- function(k_percentiles) {
  paste0(k_percentiles, "th percentile")
}

single_expert_predict <- function(test_set, expert_id) {
  # For a single expert predict the median
  prediction <- test_set |> filter(expert_id == expert_id) |> dplyr::pull(k_percentiles_to_colname(50))
  prediction
}

estimate_margin_distributions <- function(samples) {
  # Estimate the margins of the samples
  margins <- lapply(samples, function(s) {
    density(s)
  })
  margins
}

copula_fit_and_predict_fit_marginals <- function(training_set, test_set, copula_model="joe", error_metric=NULL, summarizing_function=NULL, k_percentiles=c(5,50,95)) {
  if (is.null(error_metric)) {
    error_metric <- get_ratio_error_metric()
  }
  if (is.null(summarizing_function)) {
    summarizing_function <- get_three_quantiles_summarizing_function()
  }
  nr_experts <- length(unique(training_set$expert_id))
  if (nr_experts == 0) {
    stop("No training data")
  }
  if (nr_experts == 1) {
    # If we only have one expert we predict using the median
    prediction <- single_expert_predict(test_set, training_set$expert_id[1])
    return(list(prediction=prediction, warning=c("ONE_EXPERT")))
  }

  flattened_errors <- split_dataframe_to_error_observations(training_set, error_metric, summarizing_function$f, k_percentiles)
  error_copula <- find_copula(flattened_errors, copula_model)

  # find margins

}


copula_fit_and_predict_JC_assumption <- function(training_set, test_set, copula_model="joe", interpolation="linear", error_metric=NULL, summarizing_function=NULL, k_percentiles=c(5,50,95)){
  if (is.null(error_metric)) {
    error_metric <- get_ratio_error_metric()
  }
  if (is.null(summarizing_function)) {
    summarizing_function <- get_three_quantiles_summarizing_function()
  }
  nr_experts <- length(unique(training_set$expert_id))
  if (nr_experts == 0) {
    stop("No training data")
  }
  if (nr_experts == 1) {
    # If we only have one expert we predict using the median
    prediction <- single_expert_predict(test_set, training_set$expert_id[1])
    return(list(prediction=prediction, warning=c("ONE_EXPERT")))
  }

  flattened_errors <- split_dataframe_to_error_observations(training_set, error_metric, summarizing_function$f, k_percentiles)
  error_copula <- find_copula(flattened_errors, copula_model)

  m_realizations_test <- summarizing_function$f(test_set[k_percentiles_to_colname(k_percentiles)])
  test_set <- add_0_and_100_percentiles(test_set, k_percentiles)
  distributions <- interpolate_distributions(test_set, interpolation=interpolation)

  unnorm_log_posterior <- create_log_unnormalized_posterior_JC(error_copula, distributions, error_metric, m_realizations_test)


  prediction <- mean(sample_log_unnormalized_density(unnorm_log_posterior$logDM, unnorm_log_posterior$support, 1000))
  list(prediction=prediction, warning=NULL)
}


copula_calibration_rejection_fit_and_predict <- function(training_set, test_set, alpha=0.05, copula_model="joe", interpolation="spline", dep_error_metric="rel_error", summarizing_function=NULL){
  # First reject experts with low calibration score then perform the copula method on the remaining experts
  quantile_cols <- k_percentiles_to_colname(c(5, 50, 95))
  assesment_data <- training_set[c(quantile_cols, "realization")] # the first three columns has to be the quantiles
  cal_scores <- training_set |> group_by(expert_id) |> summarise(calibration_score = {
    group_data = assesment_data[cur_group_rows(), ]
    realisations = group_data$realization
    calculateCalibrationScoreForExpert(group_data, realisations)
  })
  # Reject experts with low calibration score
  accepted_experts <- cal_scores |> filter(calibration_score >= alpha) |> pull(expert_id)
  # assert at least one expert is calirated
  below_alpha_warning <- NULL
  if (length(accepted_experts) == 0) {
    below_alpha_warning <- "ALL_BELOW_ALPHA"
    # select the expert with highest cal_score
    accepted_experts <- cal_scores |> arrange(desc(calibration_score)) |> head(1) |> pull(expert_id)
  }

  filtered_training <- training_set |> filter(expert_id %in% accepted_experts)
  filtered_test <- test_set |> filter(expert_id %in% accepted_experts)

  # Then perform the copula method
  result <- copula_fit_and_predict_JC_assumption(filtered_training, filtered_test, copula_model, interpolation, dep_error_metric, summarizing_function)

  result$warning <- append(result$warning, below_alpha_warning)
  result
}


median_average_predict <- function(test_set) {
  median_col_name = k_percentiles_to_colname(50)
  estimates <- test_set[[median_col_name]]
  prediction <- mean(estimates)
  prediction
}

create_cross_validation_sets <- function(study_data) {
  # Nests by question id and then creates all possible combinations of training and test sets
  fold_combinations <- study_data |> tidyr::nest(.by="question_id", .key="question_data") |>
    ExhaustiveFolds(1) |> dplyr::mutate(training = purrr::map(training, unnest, "question_data"), test = purrr::map(test, unnest, "question_data"))
  # returns a tibble with columns question_id, training, test. training and test are themselves also tibbles.
  fold_combinations
}


#' Extracts the percentile columns from the dataframe and returns a matrix
#'
#' @param study_data - a data frame with percentile columns. n rows
#' @param percentiles - the percentiles present in the data frame. Length d
#'
#' @returns - a matrix (nxd) with the assessment data
#' @export
#'
data_frame_to_assessment_matrix <- function(df, percentiles) {
  col_names <- k_percentiles_to_colname(percentiles)
  stopifnot(all(col_names %in% colnames(df)))

  as.matrix(df[col_names])
}

study_data_to_assessment_matrices <- function(study_data, percentiles) {
  col_names <- k_percentiles_to_colname(percentiles)
  stopifnot(all(col_names %in% colnames(study_data)))

  ordered_set <- study_data |> arrange(expert_id, question_id)
  assessments <- as.matrix(ordered_set[col_names])
  split.data.frame(assessments, study_data$expert_id)
}

global_opt_weight_predict <- function(training_set, test_set) {
  # Only works for data with quantiles at 5%, 50% and 90%
  percentiles <- c(5, 50, 95)
  training_assessments <- study_data_to_assessment_matrices(training_set, percentiles)
  realisations <- training_set |> dplyr::distinct(question_id, realization) |> arrange(question_id) |> pull(realization)
  optimal_weights <- perfWeights_opt(training_assessments, realisations, percentiles/100)

  support <- calculate_assessment_support(test_set) |> dplyr::select(L_star, U_star)
  test_assessments <- study_data_to_assessment_matrices(test_set, percentiles)
  DM <- constructDM(test_assessments, optimal_weights, NULL, support$L_star, support$U_star, percentiles/100)
  median <- DM[1,2] # single row, second entry is the median
  median
}

equal_weight_predict <- function(test_set) {
  nr_experts <- nrow(test_set)
  stopifnot(length(unique(test_set$expert_id)) == nr_experts)

  percentiles <- c(5, 50, 95)
  support_df <- calculate_assessment_support(test_set, percentiles, 0.01)
  stopifnot(nrow(support_df) == 1)
  L_star = support_df$L_star[1]
  U_star = support_df$U_star[1]

  assessment_list <- study_data_to_assessment_matrices(test_set, percentiles)

  equal_weight = 1/nr_experts
  DM <- constructDM(assessment_list, rep(equal_weight, nr_experts), NULL, L_star, U_star, percentiles/100)
  median <- DM[1,2] # single row, second entry is the median
  median
}


#' Title
#'
#' @param study_data - a data frame with columns expert_id, question_id, and the
#'  k_percentiles, realization
#' @param k_percentiles
#'
#' @export
#'
study_test_performance <- function(study_data, sim_params=NULL) {
  if (is.null(sim_params)) {
    sim_params <- default_simulation_params()
  }
  p <- sim_params

  fold_combinations <- create_cross_validation_sets(study_data)

  stats <- list(prediction = vector(mode="numeric"), realization = vector(mode="numeric"), rel_error = vector(mode="numeric"), error = vector(mode="numeric"), test_question=vector(mode="numeric"))
  warnings <- c()
  single_expert_warning <- FALSE
  below_alpha_warning <- FALSE

  for (i in seq(nrow(fold_combinations))) {
    test_set = fold_combinations$test[[i]]
    training_set = fold_combinations$training[[i]]

    if (p$prediction_method == "copula") {
      stop("Normal copula prediction method not yet implemented.")
    } else if (p$prediction_method == "copula_assumption") {
      prediction <- copula_fit_and_predict_JC_assumption(training_set, test_set, p$copula_model, p$interpolation, p$error_metric, p$summarizing_function)$prediction
    } else if (p$prediction_method == "median_average") {
      prediction <- median_average_predict(test_set)
    } else if (p$prediction_method == "copula_calibration") {
      prediction_obj <- copula_calibration_rejection_fit_and_predict(training_set, test_set, p$rejection_threshold, p$copula_model, p$interpolation, p$dep_error_metric, p$summarizing_function)
      warnings <- append(warnings, prediction_obj$warning)
      prediction <- prediction_obj$prediction
    } else if (p$prediction_method == "equal_weights") {
      prediction <- equal_weight_predict(test_set)
    } else if (p$prediction_method == "global_opt") {
      prediction <- global_opt_weight_predict(training_set, test_set)
    } else {
      stop("Unknown prediction method")
    }

    realization <- unique(test_set$realization)
    test_question <- unique(test_set$question_id)

    stats$prediction[[i]] <- prediction
    stats$realization[[i]] <- realization
    stats$rel_error[[i]] <- (prediction - realization) / realization
    stats$error[[i]] <- prediction - realization
    stats$test_question[[i]] <- test_question
  }

  # convert to tibble
  stats <- as_tibble(stats)
  # return stats and warnings
  list(stats = stats, warnings = unique(warnings))
}

default_simulation_params <- function(copula_model="joe",
                                      k_percentiles = c(5, 50, 95),
                                      interpolation="linear",
                                      prediction_method="copula_assumption",
                                      rejection_threshold=0.05,
                                      summarizing_function=NULL,
                                      error_metric=NULL) {
  if (is.null(summarizing_function)) {
    summarizing_function <- get_three_quantiles_summarizing_function()
  }
  if (is.null(error_metric)) {
    error_metric <- get_ratio_error_metric()
  }
  params <- list(copula_model=copula_model,
                 interpolation=interpolation,
                 prediction_method=prediction_method,
                 rejection_threshold=rejection_threshold,
                 summarizing_function = summarizing_function,
                 error_metric = error_metric)
  class(params) <- "simulation_parameters"
  params
}

run_analysis_per_study <- function(study_list, simulation_params=NULL) {
  if (is.null(simulation_params)) {
    simulation_params <- default_simulation_params()
  }
  warnings <- list()
  results <- list()
  for (i in seq_along(study_list)) {
    study_id <- unique(study_list[[i]]$study_id) |> head(1)
    print(paste("Running study:", study_id))
    study_data <- study_list[[i]]
    study_result <- study_test_performance(study_data, simulation_params)
    warnings[[i]] <- study_result$warnings

    study_result$stats["study_id"] <- study_id
    results[[i]] <- study_result$stats
  }
  list(results=dplyr::bind_rows(results), warnings=warnings)
}



## make function to plot list of distributions
plot_distributions <- function(distributions, which_plots="pdf&cdf") {
  df <- tibble(x_data = numeric(), y_cdf = numeric(), y_pdf = numeric(), exp_id = numeric())
  for (i in seq_along(distributions)) {
    support_width = distributions[[i]]$support[2] - distributions[[i]]$support[1]
    x_left = distributions[[i]]$support[1] - 0.1 * support_width
    x_right = distributions[[i]]$support[2] + 0.1 * support_width
    x_data <- seq(x_left, x_right, length.out = 1000)
    y_cdf <- distributions[[i]]$cdf(x_data)
    y_pdf <- distributions[[i]]$pdf(x_data)
    exp_id <- i
    df <- bind_rows(df, tibble(x_data, y_cdf, y_pdf, exp_id))
  }
  df$expert <- as.factor(df$exp_id)
  cdfplot <- df |> ggplot(aes(x = x_data, y = y_cdf, color = expert)) + geom_line()
  pdfplot <- df |> ggplot(aes(x = x_data, y = y_pdf, color = expert)) + geom_line()
  if (which_plots == "pdf&cdf") {
    return(cdfplot / pdfplot)
  } else if (which_plots == "pdf") {
    return(pdfplot)
  } else if (which_plots == "cdf")  {
    return(cdfplot)
  }
}


plot_copula_posterior <- function(posterior) {
  support = posterior$support
  x = seq(support[1], support[2], by = 0.1)
  # make ggplot
  df <- tibble(x = x, y = posterior$DM(x))
  df |> ggplot(aes(x = x, y = y)) + geom_line() + labs(y = "DM Density", x = "x")
}
