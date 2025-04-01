calculate_error <- function(data, err_metric = "error") {
  # Calculate the error of the median assessment relative to the realization.
  if (err_metric == "error") {
    data$`50th percentile` - data$realization
  } else if (err_metric == "rel_error") {
    (data$`50th percentile` - data$realization) / data$realization
  } else {
    stop("Unknown error metric")
  }
}

create_error_observations <- function(data, err_metric = "error") {
  # if we have E experts we get a E dimensional return
  # if we have N questions then we have N observations
  data$error <- calculate_error(data, err_metric)
  t(data |> split(data$question_id) |> map( ~ .$error) |> simplify2array())
}

widen_error_per_expert <- function(study_data, err_name = "error") {
  # Flatten the data per expert_id
  df <- study_data |> dplyr::select(dplyr::all_of(err_name), expert_id, question_id)
  df |> dplyr::pivot_wider(
    names_from = expert_id,
    names_prefix = "Expert ",
    values_from = dplyr::all_of(err_name)
  )
}

plot_error_between_experts <- function(study_data, err_name = "rel_error") {
  # So if we have E experts we have E * (E-1) / 2 pairs of experts
  # Flatten the data per expert_id
  wide_df <- widen_error_per_expert(study_data, err_name)
  pairs(wide_df,
        main = paste("Error between experts using column name:", err_name))
}

# error_obs is nx(d*E) matrix (n questions, d * E errors)
fit_copula <- function(error_obs, copula_model = "joe",
                       family_set = c("gaussian","indep"),
                       selcrit="mbicv", psi0=0.5) {
  res <- checkmate::check_number(nrow(error_obs), lower=10)
  if (is.character(res)) {
    checkmate::makeAssertion(error_obs, "Number of training samples must be at least 10. This comes from a hard coded limit in the rvinecopula (and vineCopula) package when fitting a copula.",
                             "nrow(error_obs)",
                             NULL)
  }
  error_length = ncol(error_obs)
  lower = NULL
  upper = NULL

  pseudo_obs <- rvinecopulib::pseudo_obs(error_obs)

  # fit a copula to the data
  if (copula_model == "joe") {
    copula_model <- copula::joeCopula(dim = error_length)
    lower = 1
  } else if (copula_model == "indep") {
    return(copula::indepCopula(dim = error_length))
  } else if (copula_model == "frank") {
    copula_model <- copula::frankCopula(dim = error_length)
    lower = 0
  } else if (copula_model == "clayton") {
    copula_model <- copula::claytonCopula(dim = error_length)
  } else if (copula_model == "gumbel") {
    copula_model <- copula::gumbelCopula(dim = error_length)
  } else if (copula_model == "t") {
    copula_model <- copula::tCopula(dim = error_length)
  } else if (copula_model == "normal") {
    copula_model <- copula::normalCopula(dim = error_length, dispstr = "un")
  } else if (copula_model == "vine") {


    cop_fit <- rvinecopulib::vinecop(pseudo_obs, family_set = family_set, cores =2, selcrit=selcrit, psi0=psi0)
    return(cop_fit)

    #return(rvinecopulib::vinecop(pseudo_obs, family_set = "onepar", cores =2))
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
    fitted_params <- copula::fitCopula(
      copula_model,
      pseudo_obs,
      method = "mpl",
      start = alpha_start,
      estimate.variance = FALSE,
      lower = lower,
      upper = upper,
      optim.control = list(factr = 1e8),
      optim.method = "L-BFGS-B"
    )
  } else{
    fitted_vine_cop <- rvinecopulib::vinecop(pseudo_obs)
    fiited_params <- copula::fitCopula(
      copula_model,
      pseudo_obs,
      method = "mpl",
      optim.control = list(factr = 1e8),
      optim.method = "L-BFGS-B"
    )
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
  data |> dplyr::pivot_longer(
    cols = dplyr::contains("percentile"),
    names_to = "k",
    values_to = "k_percentile"
  ) |>
    dplyr::mutate(k = readr::parse_number(k))
}

# Returns a list frame with the L, U, L* and U* for the assessmenet given. Must be usedd for a single study and question.
calculate_assessment_support <- function(assessments,
                                         k_percentiles = c(5, 50, 95),
                                         overshoot = 0.01) {
  # assert a single study and question
  # if expert id is a column name
  if ("question_id" %in% colnames(assessments)) {
    assessments$question_id |> unique() |> checkmate::assert_number()
  }
  if ("study_id" %in% colnames(assessments)) {
    assessments$study_id |> unique() |> checkmate::assert_number()
  }

  # Get the support of the assessments
  min_k = min(k_percentiles)
  min_k_colname = k_percentiles_to_colname(min_k)
  max_k = max(k_percentiles)
  max_k_colname = k_percentiles_to_colname(max_k)
  support_df <- assessments |> dplyr::summarise(across(all_of(min_k_colname), min, .names =
                                                  "L"),
                                         across(all_of(max_k_colname), max, .names = "U"),
                                         .groups = "drop") |> mutate(L_star = L - (U - L) * overshoot,
                                                                     U_star = U + (U - L) * overshoot)
  # single row df to vector
  support_df |> as.list()
}


add_0_and_100_percentiles <- function(data,
                                      k_percentiles = c(5, 50, 95),
                                      overshoot = 0.1) {
  group_by_vars = c("study_id", "question_id")
  support_df <- calculate_assessment_support(data, k_percentiles, overshoot)

  data |>
    dplyr::mutate(`0th percentile` = support_df$L_star,
                  `100th percentile` = support_df$U_star)
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

  error_distributions <- purrr::map(seq_along(expert_m), \(m_i) {
    force(m_i)
    F_e <- function(e) {
      distribution$cdf(error_metric$f_inverse_q(m_i, e))
    }
    f_e <- function(e) {
      distribution$pdf(error_metric$f_inverse_q(m_i, e)) * abs(error_metric$f_prime_inverse_q(m_i, e))
    }
    # calculate new support. Under the assumption that the error metric is
    # monotonic we can calculate it like below
    support1 <- error_metric$f(m_i, distribution$support[1])[1, 1]
    support2 <- error_metric$f(m_i, distribution$support[2])[1, 1]
    error_support = c(min(support1, support2), max(support1, support2))

    list(cdf = F_e,
         pdf = f_e,
         support = error_support)
  })
  error_distributions
}

test_construct_error_distribution_JC_assumption <- function() {
  # test the find median function
  expert_m <- c(-1, 2, 3)
  distribution <- list(
    cdf = function(x) {
      pnorm(x, mean = 0, sd = 1)
    },
    pdf = function(x) {
      dnorm(x, mean = 1, sd = 1)
    },
    support = c(-Inf, Inf)
  )
  error_metric <- get_linear_error_metric()
  error_dists <- construct_error_distribution_JC_assumption(distribution, error_metric, expert_m)
  # plot error dist
  x = seq(-5, 5, by = 0.01)
  y = map(error_dists, \(dist) {
    dist$pdf(x)
  })
  purrr::walk(y, \(y_i) {
    plot(x, y_i)
    par(new = TRUE)
  })

}


find_median_of_fun <- function(fun, support) {
  N = 1000
  param_mesh <- seq(support[1], support[2], length.out = N)
  # pick the guess that maximizes the posterior
  y <- fun(param_mesh)
  # calculate the cum
  cum_area <- cumsum(y) * ((support[2] - support[1]) / (N - 1))
  # find the median
  median_guess <- param_mesh[which.min(abs(cum_area - 0.5))]
  median_guess
}

k_percentiles_to_colname <- function(k_percentiles) {
  # assert that we have a numeric vector with integers
  checkmate::assert_integerish(k_percentiles)
  paste0(k_percentiles, "th percentile")
}

single_expert_predict <- function(test_set, expert_id) {
  # For a single expert predict the median
  prediction <- test_set |> filter(expert_id == expert_id) |> dplyr::pull(k_percentiles_to_colname(50))
  prediction
}

#' Calculates the unnormalized posteriors of $Q$ when conditioning on a single
#' $M_i^e$ random variable and assuming a flat prior.
#'
#' \deqn{
#' f_{Q\mid M_i^e = \dot{m_i^e}}(q) \propto f_{M_i^e \mid Q=q}(\dot{m_i^e} \mid q)
#' }
#'
#' @param error_metric Error metric of S3 class "error_metric"
#' @param error_margins d*E long list of error margins
#' @param m_matrix_observed Exd matrix of summary properties
#'
#' @returns
#'
calculate_marginal_posteriors <- function(error_margins, m_matrix_observed, error_metric) {
  # Compose the error function with the margin function
  E <- nrow(m_matrix_observed)
  d <- ncol(m_matrix_observed)
  create_pdf_closure <- function(observed_m, error_margin, error_metric, support, d_val) {
    force(observed_m)
    force(error_margin)
    force(d_val)
    force(support)
    function(q_vec) {
      ret_vec <- q_vec
      outside_mask <- q_vec <= support[1] | q_vec >= support[2]
      ret_vec[outside_mask] <- 0
      q_inside <- q_vec[!outside_mask]

      m_rep <- rep.int(observed_m, length(q_inside))
      e_values <- error_metric$f_single(m_rep, q_inside, d_val)
      error_pdf_values <- error_margin$pdf(e_values)
      error_deriv <- error_metric$f_single_prime_m(m_rep, q_inside, d_val)

      ret_vec[!outside_mask] <- abs(error_deriv) * error_pdf_values
      return(ret_vec)
    }
  }
  create_cdf_closure <- function(observed_m, error_margin, error_metric, support, d_val) {
    force(observed_m)
    force(error_margin)
    force(d_val)
    force(support)
    function(q_vec) {
      ret_vec <- q_vec
      left_mask <- q_vec <= support[1]
      right_mask <- q_vec >= support[2]
      inside_mask <- left_mask == FALSE & right_mask == FALSE
      ret_vec[left_mask] <- 0
      ret_vec[right_mask] <- 1
      q_inside <- q_vec[inside_mask]

      m_rep <- rep.int(observed_m, length(q_inside))
      e_values <- error_metric$f_single(m_rep, q_inside, d_val)
      is_inc <- error_metric$f_single_increasing_q(observed_m, d_val)
      error_cdf_values <- error_margin$cdf(e_values)
      if (is_inc == FALSE){
        error_cdf_values <- 1 - error_cdf_values
      }
      ret_vec[inside_mask] <- error_cdf_values
      return(ret_vec)
    }
  }

  new_distributions <- list()
  for (e_val in seq_len(E)) {
    for (d_val in seq_len(d)) {
      linear_index = d_E_to_linear_index(d_val, e_val, d)
      observed_m <- m_matrix_observed[e_val, d_val]
      error_margin <- error_margins[[linear_index]]

      new_support <- error_support_to_q_support(error_margin$support, observed_m, d_val, error_metric)

      pdf <- create_pdf_closure(observed_m, error_margin, error_metric, new_support, d_val)
      cdf <- create_cdf_closure(observed_m, error_margin, error_metric, new_support, d_val)


      new_distributions[[linear_index]] <- list(pdf = pdf,
                                                cdf = cdf,
                                                support = new_support,
                                                expert_id = e_val,
                                                d=d_val)
    }
  }
  new_distributions
}

get_error_margins_min_max <- function(flattened_errors) {
  checkmate::assert_matrix(flattened_errors, "numeric")
  n = nrow(flattened_errors)
  # calculate smallest and largest value per col
  purrr::array_branch(flattened_errors, 2) |> purrr::map(\(col) {
    c(min(col), max(col))
  })
}


#' Title
#'
#' @param obs_vec
#' @param support
#' @param overshoot
#' @param out_of_boundary
#' @param clamp_epsilon
#'
#' @returns list of pdf, cdf, support, approx_middle. approx middle will be
#' the mode if the shape parameters are greater than 1. Otherwise it will be
#' the mean.
#'
#' @examples
estimate_margin_beta <- function(obs_vec, support=NULL, overshoot=0.1, out_of_boundary="clamp", clamp_epsilon=1e-3) {
  checkmate::assert_numeric(support, null.ok=TRUE, len=2, any.missing=FALSE, finite=TRUE)
  min_obs <- min(obs_vec)
  max_obs <- max(obs_vec)

  if (is.null(support)) {
    obs_width <- max_obs - min_obs
    support <- c(min_obs - overshoot * obs_width, max_obs + overshoot * obs_width)
  }
  if (out_of_boundary=="clamp") {
    inside_support <- c(support[1] + clamp_epsilon, support[2] - clamp_epsilon)
    obs_vec <- pmin(pmax(obs_vec, inside_support[1]), inside_support[2])
  } else if (out_of_boundary=="discard") {
    obs_vec <- obs_vec |> purrr::discard(\(x) {
      x <= support[1] | x >= support[2]
    })
  } else {
    stop("Unknown out_of_boundary parameter")
  }

  if (any(obs_vec == support[1]) | any(obs_vec == support[2])) {
    warning("Observations on the support boundary. They will be removed to enable fitting of beta distribution.")
  }
  if (support[1] > support[2]) {
    stop("Support is not valid. The lower bound is greater than the upper bound.")
  }
  # check if there are points outside the support
  if (any(obs_vec < support[1]) | any(obs_vec > support[2])) {
    warning("Observations outside the support. They will be removed to enable fitting of beta distribution.")
  }

  if (length(obs_vec) == 0) {
    rlang::abort("No observations left after removing observations on and outside the support boundary.",
                 class="observations_outside_support")
  }
  # sum comes from product of log. The joint if error observations are independent
  neg_log_lik <- \(shape1, shape2) -sum(extraDistr::dnsbeta(obs_vec, shape1, shape2, min=support[1], max=support[2], log=TRUE))
  fit <- stats4::mle(neg_log_lik, start = list(shape1=1.5, shape2=1.5), method="L-BFGS-B", lower=list(shape1=1.001, shape2=1.001))
  params <- fit@coef
  shape1 = params[1]
  shape2 = params[2]

  pdf <- function(e_vec) {
    extraDistr::dnsbeta(e_vec, shape1, shape2, min=support[1], max=support[2])
  }

  cdf <- function(e_vec) {
    extraDistr::pnsbeta(e_vec, shape1, shape2, min=support[1], max=support[2])
  }


  if (shape1 >1 && shape2 > 1) {
    mode_standard <- (shape1 - 1)/(shape1 + shape2 - 2) # formula for mode of standard beta distribution
    scaled_mode <- support[1] + mode_standard * (support[2] - support[1])
    approx_middle <- scaled_mode
  } else {
    mean_standard <- (shape1/(shape1 + shape2)) # formula for mean of standard beta distribution
    scaled_mean <- support[1] + mean_standard * (support[2] - support[1])
    approx_middle <- scaled_mean
  }

  list(pdf = pdf, cdf = cdf, support = support, approx_middle=approx_middle)
}

estimate_margin_kde <- function(obs, support) {
  if (!is.null(support)) {
    kde <- stats::density.default(obs, bw = "SJ", kernel="gaussian", from=support[1], to=support[2])
  } else {
    kde <- stats::density.default(obs, bw = "SJ", kernel="gaussian")
  }

  x <- kde$x
  y <- kde$y
  # emperical cdf y values
  y_cdf = cumsum(y) / cumsum(y)[length(y)]

  #kde <- kde1d::kde1d(obs, xmin=xmin, xmax=xmax)
  pdf<- function (e_vec) {
    #kde1d::dkde1d(e_vec, kde)
    # linear interpolation from x and y
    approx(x, y, e_vec, method="linear", yleft=0, yright=0)$y
  }
  cdf <- function (e_vec) {
    #kde1d::pkde1d(e_vec, kde)
    approx(x, y_cdf, e_vec, method="linear", yleft=0, yright=1)$y
  }
  support <- range(x)
  approx_middle <- support[1] + (support[2]-support[1])/2
  list(pdf = pdf, cdf = cdf, support = support, approx_middle=approx_middle)
}

#' Estimate marginal distributions
#'
#' @param observations - a matrix of observations. n rows and d columns
#'
#' @returns - a list of d fitted marginal distributions
#' @export
#'
#' @examples
estimate_margins <- function(observations, supports=NULL, method="beta", overshoot=0.1, out_of_boundary="clamp") {
  checkmate::assert_matrix(observations, "numeric")
  checkmate::assert_list(supports, null.ok=TRUE, len = ncol(observations), any.missing = FALSE, types = c("numeric"))

  # if null then fill with NULL
  if (is.null(supports)) {
    supports <- rep(NULL, ncol(observations))
  }

  margins <- purrr::map(seq_len(ncol(observations)), function(i) {
    obs <- observations[, i]
    min_obs <- min(obs)
    max_obs <- max(obs)
    if (method == "beta"){
      return(estimate_margin_beta(obs, supports[[i]], overshoot, out_of_boundary=out_of_boundary))

    }
    if (method =="kde") {

      if (!is.null(supports[[i]])) {
        kde <- stats::density.default(obs, bw = "SJ", kernel="gaussian", from=supports[[i]][1], to=supports[[i]][2])
      } else {
        kde <- stats::density.default(obs, bw = "SJ", kernel="gaussian")
      }

      x <- kde$x
      y <- kde$y
      # emperical cdf y values
      y_cdf = cumsum(y) / cumsum(y)[length(y)]

      #kde <- kde1d::kde1d(obs, xmin=xmin, xmax=xmax)
      pdf<- function (e_vec) {
        #kde1d::dkde1d(e_vec, kde)
        # linear interpolation from x and y
        approx(x, y, e_vec, method="linear", yleft=0, yright=0)$y
      }
      cdf <- function (e_vec) {
        #kde1d::pkde1d(e_vec, kde)
        approx(x, y_cdf, e_vec, method="linear", yleft=0, yright=1)$y
      }
      support <- range(x)
      approx_middle <- support[1] + (support[2]-support[1])/2
      list(pdf = pdf, cdf = cdf, support = support, approx_middle=approx_middle)

    } else {stop(glue::glue("Unknown method: {method}"))}

  })
  margins
}

#' Title
#'
#' @param error_intervals list of c(a,b) intervals
#' @param m_test_matrix
#' @param error_metric  Error metric of S3 class "error_metric"
#' @param combination
#'
#' @returns
#' @export
#'
#' @examples
find_target_q_support <- function(error_intervals, m_test_matrix, error_metric, combination="closure") {
  q_supports <- error_supports_to_q_supports(error_intervals, m_test_matrix, error_metric)

  if (combination == "closure") {
    return(q_supports |> support_interval_closure())
  } else if (combination == "intersection") {
    return(q_supports |> support_intersection())
  } else {
    stop("Unknown combination method")
  }
}

target_q_support_from_test_assessments <- function(m_test_matrix) {
  support <- c(min(m_test_matrix), max(m_test_matrix))
}

target_q_support_to_error_supports <- function(target_q_support, m_test_matrix, error_metric) {
  # Calculate the error support from the target q support
  purrr::array_tree(m_test_matrix, margin=c(1,2)) |> purrr::map(\(m_row) {
    purrr::imap(m_row, \(m_observed, d_i) {
      q_support_to_error_support(target_q_support, m_observed, d_i, error_metric)
    })
  }) |> purrr::list_flatten()
}

widen_support <- function(support, overshoot, not_cross_zero=TRUE) {
  # Expand the support by overshoot percent
  checkmate::assert_numeric(support, len=2)
  checkmate::assert_numeric(overshoot, lower=0)
  support_width <- support[2] - support[1]
  new_support <- c(support[1] - overshoot * support_width, support[2] + overshoot * support_width)
  if (not_cross_zero) {
    if (support[1] > 0 && new_support[1] < 0) {
      new_support[1] <- support[1]
    }
    if (support[2] < 0 && new_support[2] > 0) {
      new_support[2] <- support[2]
    }
  }
  new_support
}


#' Title
#'
#' @param training_estimates nxExd array with the training data.
#' @param training_realizations n long vector with the realisations of the training data.
#' @param test_matrix Exd data frame with the test data. Exd order must match the training set.
#' @param copula_model string indicating copula model
#' @param error_metric Error metric of S3 class "error_metric"
#' @param summarizing_function
#' @param k_percentiles
#'
#' @returns
#' @export
#'
#' @examples
fit_and_construct_posterior <- function(training_estimates, training_realizations,
                                                 test_matrix,
                                                 copula_model = "joe",
                                                 error_metric = NULL,
                                        vine_fit_settings = list(),
                                        error_estimation_settings=list(),
                                        q_not_cross_zero=TRUE) {
  if (is.null(error_metric)) {
    error_metric <- get_ratio_error_metric()
  }
  checkmate::assert_array(training_estimates, "numeric", d=3)
  checkmate::assert_numeric(training_realizations)
  checkmate::assert_matrix(test_matrix, "numeric")
  checkmate::assert_class(error_metric, "error_metric")
  checkmate::assert_string(copula_model)
  checkmate::assert_list(vine_fit_settings)
  checkmate::assert_list(error_estimation_settings)

  domain_check <- error_metric$check_domain(abind::abind(training_estimates,
                                         test_matrix, along=1), training_realizations)
  checkmate::makeAssertion(NULL, domain_check, var.name="m and q values", collection = NULL)

  E = nrow(test_matrix)
  D = ncol(test_matrix)
  training_dim <- dim(training_estimates)
  checkmate::assert_set_equal(training_dim, c(length(training_realizations), E, D), ordered=TRUE)
  n <- training_dim[1]
  res <- checkmate::check_number(n, lower=10)
  if (is.character(res)) {
    checkmate::makeAssertion(n, "Number of training samples must be at least 10. This comes from a hard coded limit in the rvinecopula (and vineCopula) package.",
                             "dim(training_estimates)[1]",
                             NULL)
  }

  warning=NULL

  if (E == 1) {
    # If we only have one expert we predict using the median
    #prediction <- test_matrix[1,2]
    #return(list(prediction = prediction, warning = c("ONE_EXPERT")))
    warining=c("ONE_EXPERT")
  }

  errors <- assessment_array_to_errors(training_estimates, training_realizations, error_metric$f)
  flattened_errors <- flatten_3d_array_to_matrix(errors)

  target_q_support <- target_q_support_from_test_assessments(test_matrix)

  #errors_min_max <- get_error_margins_min_max(flattened_errors) # per dimension
  #target_q_support <- find_target_q_support(errors_min_max, test_matrix, error_metric)

  target_q_support <- widen_support(target_q_support, 0.1, not_cross_zero = q_not_cross_zero)

  target_error_supports <- target_q_support_to_error_supports(target_q_support, test_matrix, error_metric)

  # inject here allows us to pass the arg=value pairs of the list as arguments
  margin_distributions <- rlang::inject(estimate_margins(flattened_errors, target_error_supports, !!!error_estimation_settings))
  margin_distributions <- add_d_e_to_list(margin_distributions, D)

  error_copula <- rlang::inject(fit_copula(flattened_errors, copula_model, !!!vine_fit_settings))

  posterior<-create_log_unnormalized_posterior(error_copula,
                                    margin_distributions,
                                    error_metric,
                                    test_matrix)
  list(
    posterior = posterior,
    warning = warning,
    flattened_errors=flattened_errors,
    errors=errors,
    error_copula=error_copula,
    error_margins=margin_distributions
  )
}


copula_fit_and_predict_JC_assumption <- function(training_set,
                                                 test_set,
                                                 copula_model = "joe",
                                                 interpolation = "linear",
                                                 error_metric = NULL,
                                                 summarizing_function = NULL,
                                                 k_percentiles = c(5, 50, 95)) {
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
    return(list(prediction = prediction, warning = c("ONE_EXPERT")))
  }


  flattened_errors <- split_dataframe_to_error_observations(training_set,
                                                            error_metric,
                                                            summarizing_function$f,
                                                            k_percentiles)
  error_copula <- fit_copula(flattened_errors, copula_model)

  m_realizations_test <- summarizing_function$f(test_set[k_percentiles_to_colname(k_percentiles)])
  test_set <- add_0_and_100_percentiles(test_set, k_percentiles)
  distributions <- interpolate_distributions(test_set, interpolation = interpolation)

  unnorm_log_posterior <- create_log_unnormalized_posterior_JC(error_copula,
                                                               distributions,
                                                               error_metric,
                                                               m_realizations_test)

  list(
    posterior = unnorm_log_posterior,
    warning = NULL
  )
}


copula_calibration_rejection_fit_and_predict <- function(training_set,
                                                         test_set,
                                                         alpha = 0.05,
                                                         copula_model = "joe",
                                                         interpolation = "spline",
                                                         dep_error_metric = "rel_error",
                                                         summarizing_function = NULL) {
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
  result <- copula_fit_and_predict_JC_assumption(
    filtered_training,
    filtered_test,
    copula_model,
    interpolation,
    dep_error_metric,
    summarizing_function
  )

  result$warning <- append(result$warning, below_alpha_warning)
  result
}


median_average_predict <- function(test_set) {
  median_col_name = k_percentiles_to_colname(50)
  estimates <- test_set[[median_col_name]]
  prediction <- mean(estimates)
  prediction
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
  optimal_weights <- perfWeights_opt(training_assessments, realisations, percentiles /
                                       100)

  support <- calculate_assessment_support(test_set) |> dplyr::select(L_star, U_star)
  test_assessments <- study_data_to_assessment_matrices(test_set, percentiles)
  DM <- constructDM(
    test_assessments,
    optimal_weights,
    NULL,
    support$L_star,
    support$U_star,
    percentiles / 100
  )
  median <- DM[1, 2] # single row, second entry is the median
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

  equal_weight = 1 / nr_experts
  DM <- constructDM(
    assessment_list,
    rep(equal_weight, nr_experts),
    NULL,
    L_star,
    U_star,
    percentiles / 100
  )
  median <- DM[1, 2] # single row, second entry is the median
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
study_test_performance <- function(study_data, sim_params = NULL) {
  if (is.null(sim_params)) {
    sim_params <- default_simulation_params()
  }
  p <- sim_params

  fold_combinations <- create_cross_validation_sets(study_data)

  stats <- list(
    prediction = vector(mode = "numeric"),
    realization = vector(mode = "numeric"),
    test_question = vector(mode = "numeric"),
    posterior = list()
  )
  warnings <- c()
  single_expert_warning <- FALSE
  below_alpha_warning <- FALSE

  for (i in seq(nrow(fold_combinations))) {
    test_set = fold_combinations$test[[i]]
    training_set = fold_combinations$training[[i]]

    posterior = NULL
    prediction = NA
    if (p$prediction_method == "copula") {
      res <- fit_and_construct_posterior(
        training_set |> study_df_to_assessment_array(p$summarizing_function$f, p$k_percentiles),
        training_set |> study_df_to_realizations(),
        test_set |> study_df_single_question_to_assessment_matrix(p$summarizing_function$f, p$k_percentiles),
        p$copula_model,
        p$error_metric,
        vine_fit_settings = p$vine_fit_settings
      )
      posterior <- res$posterior
      warnings <- append(warnings, res$warning)
    } else if (p$prediction_method == "copula_assumption") {
      res <- copula_fit_and_predict_JC_assumption(
        training_set,
        test_set,
        p$copula_model,
        p$interpolation,
        p$error_metric,
        p$summarizing_function
      )
      posterior <- res$posterior
    } else if (p$prediction_method == "median_average") {
      prediction <- median_average_predict(test_set)
    } else if (p$prediction_method == "copula_calibration") {
      prediction_obj <- copula_calibration_rejection_fit_and_predict(
        training_set,
        test_set,
        p$rejection_threshold,
        p$copula_model,
        p$interpolation,
        p$dep_error_metric,
        p$summarizing_function
      )
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
    stats$test_question[[i]] <- test_question
    stats$posterior[[i]] <- posterior
  }

  # convert to tibble
  stats <- tibble::as_tibble(stats)
  # return stats and warnings
  list(stats = stats, warnings = unique(warnings))
}

default_simulation_params <- function(copula_model = "joe",
                                      k_percentiles = c(5, 50, 95),
                                      interpolation = "linear",
                                      prediction_method = "copula_assumption",
                                      rejection_threshold = 0.05,
                                      summarizing_function = NULL,
                                      error_metric = NULL,
                                      error_estimation_settings = list(),
                                      vine_fit_settings=list()) {
  if (is.null(summarizing_function)) {
    summarizing_function <- get_three_quantiles_summarizing_function()
  }
  if (is.null(error_metric)) {
    error_metric <- get_ratio_error_metric()
  }
  params <- list(
    copula_model = copula_model,
    k_percentiles = k_percentiles,
    interpolation = interpolation,
    prediction_method = prediction_method,
    rejection_threshold = rejection_threshold,
    summarizing_function = summarizing_function,
    error_metric = error_metric,
    vine_fit_settings = vine_fit_settings,
    error_estimation_settings=error_estimation_settings
  )
  class(params) <- "simulation_parameters"
  params
}

run_analysis_per_study <- function(study_list, simulation_params = NULL) {
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
  list(results = dplyr::bind_rows(results),
       warnings = warnings)
}


plot_supports <- function(list_supports) {
  p <- ggplot2::ggplot() + ggplot2::labs(x="Support", y="Expert ID")
  plot_add_supports(p, list_supports)
}

plot_add_supports <- function(p, list_supports) {
  support_df <- list_supports |> purrr::map(\(element) {
    list(left_support=element$support[1], right_support=element$support[2], d=element$d, expert_id=element$e)
  }) |> purrr::list_transpose() |> dplyr::as_tibble()


  support_df$d <- factor(support_df$d, support_df$d |> unique())
  support_df$expert_id <- factor(support_df$expert_id, support_df$expert_id |> unique())

  p + ggplot2::geom_text(data=support_df, aes(x=left_support, y=as.factor(expert_id), label="[", color=d)) +
    ggplot2::geom_text(data=support_df, aes(x=right_support, y=as.factor(expert_id), label="]", color=d))
}


#' Title
#'
#' @param errors NxExD array with the errors
#'
#' @returns
#' @export
#'
#' @examples
plot_3d_errors <- function(errors) {
  q_names <- dimnames(errors)[[1]]
  e_names <- dimnames(errors)[[2]]
  d_names <- dimnames(errors)[[3]]
  errors_flat <- aperm(errors, c(3,2,1)) |> as.vector()
  #errors_flat <- as.vector(errors)
  df <- tibble::tibble(
    error = errors_flat,
    q = rep(q_names, each = length(e_names) * length(d_names)),
    e = rep(rep(e_names, each=length(d_names)), length(q_names)),
    d = rep(d_names, length(q_names) * length(e_names))
  )
  df$q <- factor(df$q, levels = q_names)
  df$e <- factor(df$e, levels = e_names)
  df$d <- factor(df$d, levels = d_names)

  print(glue::glue("Plotting for questions {paste0(q_names, collapse=', ')}"))
  p<-df |> ggplot2::ggplot(ggplot2::aes(x=error, y=e, color=d)) +
    ggplot2::geom_point(alpha=0.9) + ggplot2::labs(x="Error", y="Expert ID")
  p
}

## make function to plot list of distributions
plot_distributions <- function(distributions, which_plots = "pdf&cdf", not_cross_zero=FALSE, D=NULL, overshoot=0.1, fix_range=NULL) {
  df <- tibble::tibble(
    x_data = numeric(),
    y_cdf = numeric(),
    y_pdf = numeric(),
    expert = numeric(),
    d = numeric(),
    row = numeric()
  )
  for (i in seq_along(distributions)) {
    if (is.null(fix_range)) {
      new_support <- widen_support(distributions[[i]]$support, overshoot = overshoot, not_cross_zero = not_cross_zero)
    } else{
      new_support <- fix_range
    }
    x_data <- seq(new_support[1], new_support[2], length.out = 1000)
    y_cdf <- distributions[[i]]$cdf(x_data)
    y_pdf <- distributions[[i]]$pdf(x_data)
    row <- i
    expert <- NA
    d <- NA
    if (is.numeric(distributions[[i]]$expert_id)) {
      expert <- distributions[[i]]$expert_id
    } else if(!is.null(D)) {
      dE <- linear_index_to_d_E(i, D)
      d <- dE$d
      expert <- dE$e
    }
    if (is.numeric(distributions[[i]]$d)) {
      d <- distributions[[i]]$d
    }
    df <- dplyr::bind_rows(df, tibble::tibble(x_data, y_cdf, y_pdf, expert, d, row))
  }
  if (anyNA(df$expert)) {
    df$row <- as.factor(df$row)
    cdfplot <- df |> ggplot2::ggplot(aes(x = x_data, y = y_cdf, color = row))
    pdfplot <- df |> ggplot2::ggplot(aes(x = x_data, y = y_pdf, color = row))
  } else if (anyNA(df$d)) {
    df$expert <- as.factor(df$expert)
    cdfplot <- df |> ggplot2::ggplot(aes(x = x_data, y = y_cdf, color = expert))
    pdfplot <- df |> ggplot2::ggplot(aes(x = x_data, y = y_pdf, color = expert))
  } else {
    df$expert <- as.factor(df$expert)
    df$d <- as.factor(df$d)
    cdfplot <- df |> ggplot2::ggplot(aes(x = x_data, y = y_cdf, color = expert, linetype=d))
    pdfplot <- df |> ggplot2::ggplot(aes(x = x_data, y = y_pdf, color = expert, linetype=d))
  }
  pdfplot <- pdfplot + ggplot2::geom_line()
  cdfplot <- cdfplot + ggplot2::geom_line()
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
