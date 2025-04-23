


#' Title
#'
#' @param training_estimates QxExD array with the training data.
#' @param training_realizations Q long vector with the realizations of the training data.
#' @param rejection_level numeric value between 0 and 1. The level of rejection.
#' @param test string indicating the test to use.
#' @param decoupler decoupler object. Only used if test is "kruskal".
#' @param min_nr_experts minimum number of experts to accept. If the number of
#' accepted experts is less than this value, the best experts are selected based on p-values.
#'
#' @returns
#' @export
#'
#' @examples
reject_experts <- function(training_estimates, training_realizations,
                           rejection_level=0.1, test="kruskal", decoupler=NULL, min_nr_experts=NULL) {
  checkmate::assert_array(training_estimates, "numeric", d=3)
  checkmate::assert_numeric(training_realizations)
  checkmate::assert_subset(test, c("kruskal", "classical_calibration", "distance_correlation"))
  checkmate::assert_number(rejection_level, lower=0, upper=1)
  checkmate::assert_count(min_nr_experts, positive=TRUE, null.ok=TRUE)
  p_vals <- p_values_test(training_estimates, training_realizations, test, decoupler)

  if (!is.null(min_nr_experts)) {
    accepted_experts <- adaptive_p_value_rejection(p_vals, rejection_level, min_nr_experts)
  } else {
    accepted_experts <- which(p_vals > rejection_level)
  }

  accepted_estimates <- training_estimates[,accepted_experts, ,drop=FALSE]

  return(list(
    accepted_estimates=accepted_estimates,
    p_values=p_vals,
    accepted_experts=accepted_experts,
    rejected_experts=setdiff(1:dim(training_estimates)[2], accepted_experts)
  ))
}

adaptive_p_value_rejection <- function(p_values, rejection_level, min_nr_experts) {
  accepted_experts=which(p_values > rejection_level)

  if (length(accepted_experts) < min_nr_experts) {
    stopifnot(length(p_values) >= min_nr_experts) # check if there are enough experts
    # Sort the p-values and select the best min_nr_experts experts
    sorted_indices <- order(p_values, decreasing=TRUE)
    accepted_experts <- sorted_indices[1:min_nr_experts]
  }

  return(accepted_experts)
}

p_values_test <- function(estimates, realizations, test="kruskal", decoupler=NULL) {
  if (test=="classical_calibration") {
    if (dim(estimates)[3] != 3) {
      stop("For classical calibration we need three quantiles for each expert.")
    }
    return(p_values_classical_calibration(estimates, realizations))
  }
  else if (test=="kruskal") {
    checkmate::assert_class(decoupler, "decoupler")
    return(p_values_kruskal(estimates, realizations, decoupler))
  }
  else if (test == "distance_correlation") {
    checkmate::assert_class(decoupler, "decoupler")
    decoupled_vals <- assessment_array_to_indep_obs(estimates, realizations, decoupler$f)
    return(p_values_distance_correlation(estimates, decoupled_vals))
  }
  else {
    stop("Unknown test.")
  }
}


# Tests the assumption that M and Z are independent
run_distance_correlation_tests  <- function(training_summaries, decoupled_values) {
  # training_summaries (QxExD array)
  # decoupled_values (QxExtilde{D} array)
  # for each expert we have a QxD and a Qx\tilde{D} matrix
  training_summaries_per_E <- purrr::array_branch(training_summaries, c(2))
  decoupled_values_per_E <- purrr::array_branch(decoupled_values, c(2))

  purrr::map2(training_summaries_per_E, decoupled_values_per_E, fixed_dcorT_test)
}

fixed_dcorT_test <- function(x, y) {
  caught_warning <- NULL
  dcor_test <- withCallingHandlers(
    {
      energy::dcorT.test(x, y)
    },
    warning = function(w) {
      if ("NaNs produced" == conditionMessage(w) && "sqrt(XX * YY)" == conditionCall(w)) {
        caught_warning <<- "neg_var_prod"
        invokeRestart("muffleWarning")  # Prevent the warning from printing
      } else if ("NaNs produced" == conditionMessage(w) && "sqrt(1 - bcR^2)" == conditionCall(w)) {
        # This is a different warning, we want to propagate it
        caught_warning <<- "bcR^2_rounding"
        invokeRestart("muffleWarning") # Prevent the warning from printing
      }
    }
  )
  if (is.null(caught_warning)) {
    return(dcor_test)
  }
  if (caught_warning == "neg_var_prod") { # this should have been implemented in the package but when there is a negative term in the square root the adjusted dcor should be zero
    dcor_test$estimate[1] <- 0
    dcor_test$statistic[1] <- 0 # T value is 0 if R adjusted is 0
    dcor_test$p.value[1] <- 0.5 # With T being asymptotically Student T distributed, then half of the distribution is above 0. We are doing a one sided test.
  } else if (caught_warning == "bcR^2_rounding") {
    # This is a different warning, caused by bcR^2 being rounded to larger than 1 when it should be capped to 1.
    # this gives us an infinite test statistic with sign determined from bcR.
    dcor_test$statistic[1] <- Inf * sign(dcor_test$estimate) # T value is Inf if R adjusted is 1)
    dcor_test$p.value[1] <- if (dcor_test$statistic[1] == Inf) 0 else 1
  }
  dcor_test
}

# Tests the assumption that M and Z are independent
p_values_distance_correlation  <- function(training_summaries, decoupled_values) {
  run_distance_correlation_tests(training_summaries, decoupled_values) |>
    purrr::map_dbl(\(test_result) {
      test_result$p.value
    })
}


# estimates is a QxEx3 array
p_values_classical_calibration <- function(estimates, realizations) {
  # Calculate the calibration score for each expert
  stopifnot(length(dim(estimates)) == 3)
  per_expert <- purrr::array_branch(estimates, c(2))

  per_expert |> purrr::map_dbl(\(expert_estimates) {
    # Calculate the error for each expert
    calculateCalibrationScoreForExpert(expert_estimates, realizations)
  })
}

# Here we test against the null hypothesis that the distribution for the decoupled
# random variable is the same for all questions.
p_values_kruskal <- function(estimates, realizations, decoupler) {
  decoupled_vals <- assessment_array_to_indep_obs(estimates, realizations, decoupler$f)
  if (dim(decoupled_vals)[3] != 1) {
    stop("The decoupler must return a single value for each expert for the kruskal test. It is a univarate test.")
  }
  decoupled_vals = abind::adrop(decoupled_vals, drop=3)
  nr_groups = dim(decoupled_vals)[1] # number of questions is the number of potential groups

  per_expert <- purrr::array_branch(decoupled_vals, c(2))

  per_expert |> purrr::map_dbl(\(expert_decoupled_values) { # Q long vector

    stats::kruskal.test(expert_decoupled_values, 1:nr_groups)$p.value
  })
}
