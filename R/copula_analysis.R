
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

bicopula_smoothing_limits <- function(family) {

  lower = -Inf
  upper = Inf
  if (family == "gaussian"){
    lower = -0.2
    upper = 0.2
  } else if (family == "frank") {
    lower = 0
    upper = 1
  } else if (family == "t") {
    # only put a limit on the first correlation parameter
    lower = matrix(c(-0.5 , -Inf), nrow=2, ncol=1)
    upper = matrix(c(0.5 , 2), nrow=2, ncol=1) # 5 is of thick the tails should be
  }

  list(lower=lower, upper=upper)
}

distance_correlation_p_value_adjacency_matrix <- function(decoupler_matrix) {
  # We want to evaluate all pairs of columns in the decoupler matrix
  # and calculate the distance correlation between them

  # get all pairs of columns
  N = nrow(decoupler_matrix)
  D = ncol(decoupler_matrix)
  adj_matrix = matrix(0, nrow=D, ncol=D)
  for (i in seq_len(D)) {
    # loop through triangle of the matrix
    for (j in seq_len(D)) {
      if (i <= j) { # skip when row is smaller or equal than column.
        next
      }
      # calculate the distance correlation between the two columns
      # and store it in the adjacency matrix
      adj_matrix[i,j] <- energy::dcorT.test(decoupler_matrix[,i], decoupler_matrix[,j])$p.value
      adj_matrix[j,i] <- adj_matrix[i,j]
    }
  }

  adj_matrix
}

calculate_adjacency_matrix <- function(decoupler_matrix, method="kendall", threshold=0.7) {
  # we want to find groups of errors that are disjoint
  # determine the correlation between experts
  if (method == "kendall") {
    adj_matrix <- copula::corKendall(decoupler_matrix)
  } else if (method == "pearson") {
    adj_matrix <- stats::cor(decoupler_matrix, method="pearson")
  } else if (method == "distance_correlation") {
    adj_matrix <- distance_correlation_p_value_adjacency_matrix(decoupler_matrix)
    # we want to build edges between dependent experts. That is, experts that are not independent and thus have a low p value
    # Thus we invert the p value
    adj_matrix <- 1 - adj_matrix
  }
  else {
    stop(paste("Unknown method", method))
  }
  adjacency_matrix <- adj_matrix
  diag(adjacency_matrix) <- 0
  # set values with abs value less than threshold to 0 and other to 1
  adjacency_matrix[abs(adj_matrix) < threshold] <- 0
  adjacency_matrix[abs(adj_matrix) >= threshold] <- 1
  adjacency_matrix
}

partition_errors_disjoint <- function(cdf_matrix, method="kendall", threshold=0.7) {
  stopifnot(length(dim(cdf_matrix)) == 2)
  # we want to find groups of errors that are disjoint
  # determine the correlation between experts
  adjacency_matrix <- calculate_adjacency_matrix(cdf_matrix, method=method, threshold=threshold)
  graph <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", diag = FALSE)
  igraph::components(graph)
}

combine_singletons <- function(partition) {
  # Original membership vector
  # Original membership vector
  membership <- partition$membership

  # Compute sizes
  group_sizes <- table(membership)

  # Identify which groups are singletons
  singleton_groups <- as.numeric(names(group_sizes[group_sizes == 1]))

  # Create a new membership vector
  # All singleton groups are mapped to a temporary group ID (e.g. 999)
  temp_membership <- membership
  temp_membership[membership %in% singleton_groups] <- -1

  # Now remap group IDs to consecutive integers starting from 1
  # unique groups, in order of appearance
  #new_ids <- match(temp_membership, unique(temp_membership))

  # Compute updated partition structure
  new_partition <- list(
    membership = temp_membership, #, new_ids,
    csize = as.vector(table(temp_membership)),
    no = length(unique(temp_membership))
  )

  new_partition
}

partition_and_fit_copulas <- function(cdf_matrix, copula_model, threshold, fit_settings, method="kendall") {
  partionings <- partition_errors_disjoint(cdf_matrix, method=method, threshold=threshold)
  partionings <- combine_singletons(partionings)
  groups <- igraph::groups(partionings)
  copulas <- purrr::map2(groups, names(groups), \(group, name) {
    group_obs <- cdf_matrix[,group,drop=FALSE]
    if (name=="-1") {
      fit_copula(group_obs, "indep")
    } else {
      rlang::inject(fit_copula(group_obs, copula_model, !!!fit_settings))
    }
  })

  density <- function(u_mat, log=FALSE) { # u_mat is a vec that is dE long or a nx(d*E) matrix
    if (is.vector(u_mat)) {
      u_mat <- matrix(u_mat, nrow=1)
    }
    stopifnot(is.matrix(u_mat))
    # use membership info to know what copula to feed what subsections of the u vec
    densities <- purrr::map2(groups, seq_len(length(groups)), \(group, i) {
      u_subset <- u_mat[,group, drop=FALSE]
      copulas[[i]]$density(u_subset, log=log)
    }) # list of length nr of groups. Each element is a vector of same length as number of rows or u_mat
    # product of densities
    purrr::reduce(densities, \(x,y) x * y)
  }

  distribution <- function(u_mat) { # u_mat is a vec that is dE long or a nx(d*E) matrix
    if (is.vector(u_mat)) {
      u_mat <- matrix(u_mat, nrow=1)
    }
    stopifnot(is.matrix(u_mat))
    # use membership info to know what copula to feed what subsections of the u vec
    distributions <- purrr::imap(groups, \(group, i) {
      u_subset <- u_mat[,group, drop=FALSE]
      copulas[[i]]$distribution(u_subset)
    }) # list of length nr of groups. Each element is a vector of same length as number of rows or u_mat
    # product of densities
    purrr::reduce(distributions, \(x,y) x * y)
  }
  list(
    density = density,
    distribution = distribution,
    copulas = copulas
  )
}

wrap_copula <- function(copula, method=NULL) {
  if (inherits(copula, "vinecop")) {
    density_function <- function(u_vec, log=FALSE) {
      if (log) {
        log(rvinecopulib::dvinecop(u_vec, copula))
      } else {
        rvinecopulib::dvinecop(u_vec, copula)
      }
    }
    distribution_function <- function(u_vec) {
      rvinecopulib::pvinecop(u_vec, copula)
    }

    ret = list(
      density = density_function,
      distribution = distribution_function,
      copula=copula
    )
  } else if (inherits(copula, "copula") || inherits(copula, "indepCopula")) {
    ret = list(
      density = \(u, log=FALSE) copula::dCopula(u, copula, log=log),
      distribution = \(u) copula::pCopula(u, copula),
      copula=copula
    )
  } else {
    stop(paste("Unknown copula type", class(copula)))
  }
  if (!is.null(method)) {
    ret$method = method
  }
  ret
}

get_wrapped_independence_copula <- function(D, method="indep") {
  # create a independence copula
  indep_copula <- copula::indepCopula(dim = D)
  wrap_copula(indep_copula, method=method)
}

fit_hiarchical_copula <- function(pseudo_obs, eta=1, recover_numerical_failure=FALSE, eps=1e-10) {
  # fit a hierarchical copula
  pseudo_obs <- pmin(pmax(pseudo_obs, eps), 1 - eps) # prevent Inf and -Inf.
  normal_space <- qnorm(pseudo_obs) # stan model expects this transformation already to be made
  N = nrow(pseudo_obs)
  D = ncol(pseudo_obs)


  data_list <- list(
    N = N,
    D = D,
    eta = eta,
    normal_latent=normal_space
  )

  # load copula_hiarch_model.stan
  mod <- load_stan_copula_model()
  fit_map <- tryCatch({
     rstan::optimizing(
      mod,
      data = data_list #, verbose = TRUE
    )
    }, error = function(e) {
      if (recover_numerical_failure) {
        rlang::inform(paste0("Stan model did not converge. Error message is: ", e$message, ". Recovering to independence copula."),
                     class="stan_optimization_warning")
        return("recover")
      } else {
        rlang::abort(paste0("Stan model did not converge. Error message is: ", e$message),
                     class="stan_optimization_error", parent=e)
      }
    }, warning = function(w) {
      if (recover_numerical_failure) {
        rlang::inform(paste0("Stan model did not converge. Warning message is: ", w$message, ". Recovering to independence copula."),
                      class="stan_optimization_warning")
        return("recover")
      } else {
        rlang::warn(paste0("Stan model did not converge. Warning message is: ", w$message),
                    class="stan_optimization_warning")
    }
  })

  if (is.character(fit_map) && fit_map == "recover") {
    # Return an independence copula
    return(get_wrapped_independence_copula(D, method="recovered hierarchical"))
  }

  if (fit_map$return_code != 0) {
    rlang::abort(paste0("Stan model did not converge. Warning code is: ", fit_map$return_code),
          "stan_optimization_error")
  }
  # extract the parameters with names Sigma[i,j]
  nr_params = D*D*2
  # L and Sigma are returned and the first half of the parameters are the L values
  sigma_values <- fit_map$par[(nr_params/2+1):nr_params]
  sigma_matrix <- matrix(sigma_values, nrow = D, ncol = D)

  fitted_norm_cop <- copula::normalCopula(copula::P2p(sigma_matrix), dim = D, dispstr="un")

  wrapped_copula <- wrap_copula(fitted_norm_cop)
  wrapped_copula$eta <- eta
  wrapped_copula$method <- "hierarchical"
  wrapped_copula
}

# pseudo_obse rror_obs is nx(d*E) matrix (n questions, d * E errors)
fit_copula <- function(pseudo_obs, copula_model = "joe",
                       family_set = c("gaussian","indep"),
                       selcrit="mbicv", psi0=0.5,
                       copula_fit_method="itau",
                       threshold=0.5, eta=1, recover_numerical_failure=FALSE) {
  if (ncol(pseudo_obs) == 1) {
    indep1DCopula = copula::indepCopula(dim=1)

    return(wrap_copula(indep1DCopula))
  }
  if (copula_model == "indep") {
    return(wrap_copula(copula::indepCopula(dim = ncol(pseudo_obs))))
  }

  if (copula_model == "hierarchical") {
    # fit a hierarchical copula
    return(fit_hiarchical_copula(pseudo_obs, eta=eta, recover_numerical_failure=recover_numerical_failure))
  }

  res <- checkmate::check_number(nrow(pseudo_obs), lower=10)
  if (is.character(res)) {
    checkmate::makeAssertion(pseudo_obs, "Number of training samples must be at least 10. This comes from a hard coded limit in the rvinecopula (and vineCopula) package when fitting a copula.",
                             "nrow(pseudo_obs)",
                             NULL)
  }
  error_length = ncol(pseudo_obs)
  lower = NULL
  upper = NULL
  #pseudo_obs <- rvinecopulib::pseudo_obs(error_obs)

  # fit a copula to the data
  if (copula_model == "joe") {
    copula_model <- copula::joeCopula(dim = error_length)
    lower = 1
    upper = 18.7 # 18.7 Corresponds to an Kendall tau of 0.9. For numerical reasons we limit it to this.
  } else if (copula_model == "frank") {
    copula_model <- copula::frankCopula(dim = error_length)
    lower = 0
    upper = 38.28 # 38.28 Corresponds to an Kendall tau of 0.9. For numerical reasons we limit it to this.
  } else if (copula_model == "clayton") {
    copula_model <- copula::claytonCopula(dim = error_length)
    lower = 0
    upper = 18
  } else if (copula_model == "gumbel") {
    copula_model <- copula::gumbelCopula(dim = error_length)
    lower = 1
    upper = 10
  } else if (copula_model == "t") {
    copula_model <- copula::tCopula(dim = error_length)
    lower = -9.8
    upper = 9.8
  } else if (copula_model == "normal") {
    copula_model <- copula::normalCopula(dim = error_length, dispstr = "un")
    lower = -9.8
    upper = 9.8
    copula_fit_method="itau"
  } else if (copula_model == "vine") {


    cop_fit <- rvinecopulib::vinecop(pseudo_obs, family_set = family_set,
                                     cores =2, selcrit=selcrit, psi0=psi0,
                                     threshold = threshold)
    cop_fit$pair_copulas <- purrr::modify_tree(cop_fit$pair_copulas,
                                       leaf = \(bicop) {
                                         limits <- bicopula_smoothing_limits(bicop$family)
                                         bicop$parameters = pmax(bicop$parameters, limits$lower)
                                         bicop$parameters = pmin(bicop$parameters, limits$upper)
                                         bicop
                                       })

    return(wrap_copula(cop_fit))
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
    fitted_params <- copula::fitCopula(
      copula_model,
      pseudo_obs,
      method = copula_fit_method,
      estimate.variance = FALSE,
      optim.control = list(factr = 1e8),
      optim.method = "L-BFGS-B",
      lower=lower,
      upper=upper
    )
  }

  wrap_copula(fitted_params@copula)
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

add_0_and_100_percentiles_matrix <- function(quantiles, overshoot=0.1, support_restriction=NULL) {
  checkmate::assert_matrix(quantiles, "numeric")
  checkmate::assert_number(overshoot, lower=0)
  checkmate::assert_string(support_restriction, null.ok=TRUE)
  support <- range(quantiles)
  support <- widen_support(support, overshoot=overshoot, support_restriction = support_restriction)
  # add the 0th and 100th percentiles
  cbind(support[1], quantiles, support[2])
}

add_0_and_100_percentiles_matrix_per_row <- function(quantiles, overshoot=0.1, support_restriction=NULL) {
  checkmate::assert_matrix(quantiles, "numeric")
  checkmate::assert_number(overshoot, lower=0)
  checkmate::assert_string(support_restriction, null.ok=TRUE)
  supports <- widen_supports(quantiles[,c(1, ncol(quantiles))], overshoot=overshoot, support_restriction = support_restriction)
  # add the 0th and 100th percentiles
  quantiles <- cbind(supports[,1], quantiles, supports[,2])
  quantiles
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

mean_variance_to_beta_params <- function(mean, variance, L, U, epsilon=1e-6) {
  stopifnot(variance < 0.25) # variance must be less than 0.25 to be a beta distribution
  # mean = alpha / (alpha + beta)
  # variance = alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
  # solve for alpha and beta
  # (mean-U)*(mean-L)
  #tmp_term <- (L * U - L * mean - U * mean + mean^2 + variance) / (variance * (L - U))
  #alpha <- (L-mean) *tmp_term
  #beta <- (mean - U) * tmp_term
  if (mean < L || mean > U) {
    stop(glue::glue("Mean {round(mean, digits=3)} is outside of the support ({L},{U})."))
  }
  max_variance <- (mean-L) * (U-mean)
  if (variance >= max_variance) {
    warning(glue::glue("Variance {round(variance, digits=3)} is too large with mean {round(mean, digits=3)} and support ({L},{U}). Setting variance to (mean-L) * (U-mean)-epsilon={round(max_variance-epsilon, digits=3)}"))
    variance <- max_variance - epsilon
  }

  width = U-L
  mean_scaled <- (mean - L) / width
  variance_scaled <- variance / (width^2)

  alpha = mean_scaled*(mean_scaled * (1-mean_scaled)/variance_scaled - 1)
  beta = alpha/mean_scaled - alpha

  alpha <- pmax(alpha, epsilon) # rounding errors
  beta <- pmax(beta, epsilon) # rounding errors
  c(alpha, beta)

  #alpha <- mean * (mean * (1 - mean) / variance - 1)
  #beta <- (1 - mean) * (mean * (1 - mean) / variance - 1)
  #alpha <- pmax(alpha, min_value) # rounding errors
  #beta <- pmax(beta, min_value) # rounding errors
  #c(alpha, beta) # shape, rate
}

mean_variance_to_gamma_params <- function(mean, variance) {
  # mean = alpha / beta
  # variance = alpha / (beta^2)
  # solve for alpha and beta
  alpha <- mean^2 / variance
  beta <- mean / variance
  c(alpha, beta) # shape, rate
}

estimate_margin_beta_prior <- function(support, beta_mean, beta_var) {
  beta_vars <- mean_variance_to_beta_params(beta_mean, beta_var, support[1], support[2])

  pdf <- function(e_vec, log=FALSE) {
    extraDistr::dnsbeta(e_vec, beta_vars[1], beta_vars[2], min=support[1], max=support[2], log=log)
  }
  cdf <- function(e_vec, log=FALSE) {
    extraDistr::pnsbeta(e_vec, beta_vars[1], beta_vars[2], min=support[1], max=support[2], log.p=log)
  }
  list(pdf = pdf, cdf = cdf, support = support, approx_middle=beta_approx_middel(beta_vars[1], beta_vars[2], support[1], support[2]))
}

estimate_margin_beta_hiarch <- function(obs_vec, support, beta_mean=0.5, beta_var=1/12, prior_var=0.1, clamp_epsilon=0.001, out_of_boundary="clamp") {
  beta_param_centers <- mean_variance_to_beta_params(beta_mean, beta_var, support[1], support[2]) # Wait I am doing this a bit wrong... Needs to be adjusted for the scaling
  A_prior_params <- mean_variance_to_gamma_params(beta_param_centers[1], prior_var)
  B_prior_params <- mean_variance_to_gamma_params(beta_param_centers[2], prior_var)

  support_width <- support[2] - support[1]
  # if (support_width != 1) {
  #   stop(glue::glue("Support width must be 1. Not {support_width}."))
  # }
  if (out_of_boundary=="clamp") {
    boundry_offset <- support_width * clamp_epsilon
    inside_support <- c(support[1] + boundry_offset, support[2] - boundry_offset)
    obs_vec <- pmin(pmax(obs_vec, inside_support[1]), inside_support[2])
  }
  #zero_to_one_obs <- (obs_vec - support[1])/support_width

  model <- load_stan_beta_model()
  N = length(obs_vec)
  data_list <- list(
    N = N,
    Z = obs_vec,
    L=support[1],
    U=support[2],
    alpha_a = A_prior_params[1],
    beta_a = A_prior_params[2],
    alpha_b = B_prior_params[1],
    beta_b = B_prior_params[2]
  )

  fit_map <- rstan::optimizing(
    model,
    #verbose=TRUE,
    data = data_list # , check_data = FALSE - not sure if I should have this or not.
  )
  fit_map

  map_a <- fit_map$par[[1]]
  map_b <- fit_map$par[[2]]

  pdf <- function(e_vec, log=FALSE) {
    extraDistr::dnsbeta(e_vec, map_a, map_b, min=support[1], max=support[2], log=log)
  }

  cdf <- function(e_vec) {
    extraDistr::pnsbeta(e_vec, map_a, map_b, min=support[1], max=support[2])
  }

  list(pdf = pdf, cdf = cdf, support = support, approx_middle=beta_approx_middel(map_a, map_b, support[1], support[2]))
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
  support_width = support[2] - support[1]
  if (out_of_boundary=="clamp") {
    boundry_offset <- support_width * clamp_epsilon
    inside_support <- c(support[1] + boundry_offset, support[2] - boundry_offset)
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
  fit <- stats4::mle(neg_log_lik, start = list(shape1=1, shape2=1), method="L-BFGS-B", lower=list(shape1=0.001, shape2=0.001))
  params <- fit@coef
  shape1 = params[1]
  shape2 = params[2]

  pdf <- function(e_vec) {
    extraDistr::dnsbeta(e_vec, shape1, shape2, min=support[1], max=support[2])
  }

  cdf <- function(e_vec) {
    extraDistr::pnsbeta(e_vec, shape1, shape2, min=support[1], max=support[2])
  }

  list(pdf = pdf, cdf = cdf, support = support, approx_middle=beta_approx_middel(shape1, shape2, support[1], support[2]))
}

beta_approx_middel <- function(shape1, shape2, L, U) {
  if (shape1 > 1 && shape2 > 1) {
    mode_standard <- (shape1 - 1)/(shape1 + shape2 - 2) # formula for mode of standard beta distribution
    mode_general <- mode_standard * (U - L) + L
    return(mode_general)
  } else {
    mean_standard <- (shape1/(shape1 + shape2)) # formula for mean of standard beta distribution
    mean_general <- mean_standard * (U - L) + L
    return(mean_general)
  }

}

estimate_margin_kde <- function(obs, support, bw="SJ", package="stats") {
  if (!is.null(support)) {
    if (package=="stats") {
      kde <- stats::density.default(obs, bw = bw, kernel="gaussian", from=support[1], to=support[2])
    } else if (package=="ks") {
      kde <- ks::kde.boundary(obs, xmin=support[1], xmax=support[2], h=bw, compute.cont=FALSE)
    } else {
      stop("Unknown package for kde")
    }
  } else {
    kde <- stats::density.default(obs, bw = bw, kernel="gaussian")
  }

  if (package=="ks") {
    x <- kde$eval.points
    y <- kde$estimate
  } else {
    x <- kde$x
    y <- kde$y
  }
  # emperical cdf y values
  y_cdf = cumsum(y) / cumsum(y)[length(y)]

  #kde <- kde1d::kde1d(obs, xmin=xmin, xmax=xmax)
  pdf<- function (e_vec) {
    # linear interpolation from x and y
    approx(x, y, e_vec, method="linear", yleft=0, yright=0)$y
  }
  # because of the cut off with the support, it is not guaranteed the area is 1 so we normalize it here.
  if (package=="ks") {
    area = 1 # advantage of ks package is that it is normalized properly even when boundary is set
  } else {
    area <-  integrate(pdf, min(x), max(x))$value
  }
  pdf_normalized <- function (e_vec) {
    pdf(e_vec) / area
  }

  cdf <- function (e_vec) {
    #kde1d::pkde1d(e_vec, kde)
    approx(x, y_cdf, e_vec, method="linear", yleft=0, yright=1)$y
  }
  support <- range(x)
  approx_middle <- support[1] + (support[2]-support[1])/2
  list(pdf = pdf_normalized, cdf = cdf, support = support, approx_middle=approx_middle)
}

#' Estimate marginal distributions
#'
#' @param observations - a matrix of observations. n rows and d columns
#'
#' @returns - a list of d fitted marginal distributions
#' @export
#'
#' @examples
estimate_margins <- function(observations, supports=NULL, method="beta", overshoot=0.1, out_of_boundary="clamp", bw="SJ", beta_mean=0.5, beta_var=1/12, prior_var=0.1) {
  checkmate::assert_matrix(observations, "numeric")
  nr_dims = ncol(observations)
  checkmate::assert(
    checkmate::check_list(supports, null.ok=TRUE, len = ncol(observations), any.missing = FALSE, types = c("numeric")),
    checkmate::check_numeric(supports, len = 2, any.missing = FALSE, finite=TRUE),
  )
  checkmate::assert(
    checkmate::check_numeric(beta_mean, len = nr_dims, lower=0, upper=1),
    checkmate::check_number(beta_mean, lower=0, upper=1)
  )
  checkmate::assert(
    checkmate::check_numeric(beta_var, len = nr_dims, lower=0, upper=0.25), # variance of beta distributions are always less than 0.25
    checkmate::check_number(beta_var, lower=0, upper=0.25)
  )
  if (is.numeric(supports)) { # if a single support, replicate it for all columns
    supports <- rep(list(supports), ncol(observations))
  }
  else if (is.null(supports)) {
    # if null then fill with NULL
    supports <- rep(list(NULL), ncol(observations))
  }
  if (length(beta_mean) == 1) {
    beta_mean <- rep(beta_mean, nr_dims)
  }
  if (length(beta_var) == 1) {
    beta_var <- rep(beta_var, nr_dims)
  }

  margins <- purrr::map(seq_len(ncol(observations)), function(i) {
    obs <- observations[, i]
    if (method == "beta_MLE" || method=="beta"){
      return(estimate_margin_beta(obs, supports[[i]], overshoot, out_of_boundary=out_of_boundary))
    } else if (method =="kde") {
      return(estimate_margin_kde(obs, supports[[i]], bw=bw))
    } else if (method == "beta_MAP") {
      return(estimate_margin_beta_hiarch(obs, supports[[i]], beta_mean=beta_mean[i], beta_var=beta_var[i], prior_var=prior_var))
    } else if (method == "uniform") {
      if (is.null(supports[[i]])) {
        supports[[i]] <- range(obs) |> widen_support(overshoot=overshoot)
      }
      return(estimate_margin_uniform(supports[[i]]))
    } else if (method == "beta_prior") {
      return(estimate_margin_beta_prior(supports[[i]], beta_mean[i], beta_var[i]))
    }
    else {
      stop(glue::glue("Unknown method: {method}"))
    }

  })
  margins
}

uniform_prediction_posterior <- function(test_matrix, overshoot=0.1, support_restriction=NULL) {
  support <- range(test_matrix) |> widen_support(overshoot=0.1, support_restriction=support_restriction)
  uniform_dist <- estimate_margin_uniform(support)

  logDM <- function(q) {
    uniform_dist$pdf(q, log=TRUE)
  }
  uniform_dist$logDM <- logDM
  return(
    uniform_dist
  )

}

estimate_margin_uniform <- function(support) {
  pdf <- function(z, log=FALSE) {
    dunif(z, min=support[1], max=support[2], log=log)
  }
  cdf <- function(z, log=FALSE) {
    punif(z, min=support[1], max=support[2], log.p=log)
  }
  list(
    pdf=pdf,
    cdf=cdf,
    mean=mean(support),
    median=mean(support),
    support=support,
    approx_middle=median
  )
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

target_q_support_to_independence_supports <- function(target_q_support, m_test_matrix, decoupler) {
  decoup_values <- decoupler$f(target_q_support, m_test_matrix) # 2xExD array
  purrr::array_branch(decoup_values, c(2,3)) |> purrr::map(\(decoup_values) {
    c(min(decoup_values), max(decoup_values))
  })
}

widen_supports <- function(supports_matrix, overshoot, support_restriction=NULL) {
  # Expand the support by overshoot percent
  checkmate::assert_subset(support_restriction, c("non_negative", "strict_positive", NULL))
  checkmate::assert_matrix(supports_matrix, "numeric", ncols=2)
  checkmate::assert_numeric(overshoot, lower=0)
  support_widths <- supports_matrix[,2] - supports_matrix[,1]
  new_left <- supports_matrix[,1] - overshoot * support_widths
  new_right <- supports_matrix[,2] + overshoot * support_widths
  if (is.null(support_restriction)) {
    return(cbind(new_left, new_right))
  }
  if (support_restriction == "strict_positive") {
    if (any(supports_matrix[,1] <= 0)) {
      stop("All supports are not strictly positive. Cannot widen support with 'strict_positive' restriction.")
    }
    which_non_positive <- which(new_left <= 0)
    new_left[which_non_positive] <- new_left[which_non_positive] * (1-overshoot)
  } else if (support_restriction == "non_negative") {
    if (any(supports_matrix[,1] < 0)) {
      stop("All support are not non negative. Cannot widen support with 'non_negative' restriction.")
    }
    which_non_positive <- which(new_left < 0)
    new_left[which_non_positive] <- 0
  }
  cbind(new_left, new_right)
}

widen_support <- function(support, overshoot, support_restriction=NULL) {
  # Expand the support by overshoot percent
  checkmate::assert_subset(support_restriction, c("non_negative", "strict_positive", NULL))
  checkmate::assert_numeric(support, len=2, sorted=TRUE)
  checkmate::assert_numeric(overshoot, lower=0)
  support_width <- support[2] - support[1]
  new_support <- c(support[1] - overshoot * support_width, support[2] + overshoot * support_width)
  if (is.null(support_restriction)) {
    return(new_support)
  }
  if (support_restriction == "strict_positive") {
    if (support[1] <= 0) {
      stop("Support is not positive. Cannot widen support with 'strict_positive' restriction.")
    }
    if (new_support[1] <= 0) {
      new_support[1] <- support[1] * (1-overshoot)
    }
  } else if (support_restriction == "non_negative") {
    if (support[1] < 0) {
      stop("Support is negative. Cannot widen support with 'non_negative' restriction.")
    }
    if (new_support[1] < 0) {
      new_support[1] <- 0
    }
  }
  new_support
}

ideal_means_vars_non_finite_to_uniform <- function(means_vars, supports) {
  # if the means and vars are not finite, set them to uniform
  target_mean <- supports |> purrr::map_dbl(mean) |> matrix(nrow=nrow(means_vars$means), ncol=ncol(means_vars$means), byrow=TRUE)
  target_var <- supports |> purrr::map_dbl(\(support) (support[2] - support[1])^2/12) |> matrix(nrow=nrow(means_vars$vars), ncol=ncol(means_vars$vars), byrow=TRUE)
  means_vars$means[!is.finite(means_vars$means)] <- target_mean[!is.finite(means_vars$means)]
  means_vars$vars[!is.finite(means_vars$vars)] <- target_var[!is.finite(means_vars$vars)]
  means_vars
}

#' Title
#'
#' @param training_estimates
#' @param training_realizations
#' @param test_matrix
#' @param decoupler
#' @param copula_model
#' @param vine_fit_settings List of arguments to pass to the `fit_copula` function together with the copula_model argument.
#' @param error_estimation_settings List of arguments to pass to the `estimate_margins` function
#' @param q_support_restriction
#' @param rejection_threshold
#' @param rejection_test
#'
#' @returns
#' @export
#'
#' @examples
fit_and_construct_posterior_indep <- function(training_estimates, training_realizations,
                                           test_matrix,
                                           decoupler,
                                           copula_model = "joe",
                                           vine_fit_settings = list(),
                                           error_estimation_settings=list(),
                                           q_support_restriction=NULL,
                                           q_support_overshoot=0.1,
                                           rejection_threshold=NULL,
                                           rejection_min_experts=3,
                                           rejection_test="distance_correlation",
                                           connection_threshold=NULL,
                                           connection_metric="kendall",
                                           recover_numerical_failure=FALSE) {
  checkmate::assert_array(training_estimates, "numeric", d=3)
  checkmate::assert_numeric(training_realizations)
  checkmate::assert_matrix(test_matrix, "numeric")
  checkmate::assert_class(decoupler, "decoupler")
  checkmate::assert_string(copula_model)
  checkmate::assert_list(vine_fit_settings)
  checkmate::assert_list(error_estimation_settings)
  checkmate::assert_subset(q_support_restriction, c("non_negative", "strict_positive", NULL))
  checkmate::assert_number(q_support_overshoot, lower=0, finite=TRUE)
  checkmate::assert_number(rejection_threshold, null.ok=TRUE, lower=0, upper=1)
  checkmate::assert_count(rejection_min_experts, positive=TRUE)
  checkmate::assert_subset(rejection_test, c("kruskal", "classical_calibration", "distance_correlation"))
  checkmate::assert_number(connection_threshold, null.ok=TRUE, lower=0, upper=1)

  # TODO implement a domain check prior?
  # domain_check <- decoupler$check_domain(abind::abind(training_estimates,
  #                                                        test_matrix, along=1), training_realizations)
  #checkmate::makeAssertion(NULL, domain_check, var.name="m and q values", collection = NULL)

  E = nrow(test_matrix)
  D = ncol(test_matrix)
  training_dim <- dim(training_estimates)
  checkmate::assert_set_equal(training_dim, c(length(training_realizations), E, D), ordered=TRUE)
  res <- checkmate::check_number(training_dim[1], lower=10)
  if (is.character(res) && !(copula_model %in% c("indep", "hierarchical"))) {
    stop("Number of training samples ('dim(training_estimates)[1]') must be at least 10. This comes from a hard coded limit in the rvinecopula (and vineCopula) package.")
  }

  if (!is.null(rejection_threshold)) {
    rejection_results <- reject_experts(training_estimates, training_realizations, rejection_threshold, test=rejection_test,
                   decoupler=decoupler, min_nr_experts=rejection_min_experts)
    if (length(rejection_results$accepted_experts) == 0) {
      stop("No experts accepted. Please check the rejection threshold or set rejection_min_experts > 0.")
    }
    training_estimates <- rejection_results$accepted_assessments
    test_matrix <- test_matrix[rejection_results$accepted_experts,, drop=FALSE]
  }

  warning=NULL

  if (E == 1) {
    warining=c("ONE_EXPERT")
  }

  decouple_values <- assessments_to_decoupler_observations(training_estimates, training_realizations, decoupler$f)
  decouple_flattened <- flatten_3d_array_to_matrix(decouple_values)

  target_q_support <- target_q_support_from_test_assessments(test_matrix)

  widend_support <- widen_support(target_q_support, q_support_overshoot,  support_restriction=q_support_restriction)

  decoupler_support <- target_q_support_to_independence_supports(widend_support, test_matrix, decoupler)

  if (is.function(decoupler$ideal_mean_var)) {
    ideal_vars <- decoupler$ideal_mean_var(test_matrix) |> ideal_means_vars_non_finite_to_uniform(supports=decoupler_support)
    ideal_means <- ideal_vars$means |> flatten_matrix_row_by_row()
    ideal_vars <- ideal_vars$vars |> flatten_matrix_row_by_row()
    error_estimation_settings$beta_mean <- ideal_means
    error_estimation_settings$beta_var <- ideal_vars
  }
  # inject here allows us to pass the arg=value pairs of the list as arguments
  margin_distributions <- rlang::inject(estimate_margins(decouple_flattened, decoupler_support, !!!error_estimation_settings))
  D_tilde = dim(decouple_values)[3]
  margin_distributions <- add_d_e_to_list(margin_distributions, D_tilde)
  # convert the decoupled values to cdf values using the estimated margins
  cdf_values <- purrr::imap(margin_distributions, \(margin, i) margin$cdf(decouple_flattened[,i])) |> do.call(what=cbind)

  if (is.null(connection_threshold)){
    decoupler_copula <- rlang::inject(fit_copula(cdf_values, copula_model, !!!vine_fit_settings))
  } else {
    decoupler_copula <- partition_and_fit_copulas(cdf_values, copula_model, connection_threshold, vine_fit_settings)
  }

  posterior<-create_log_unnormalized_posterior_indep(decoupler_copula,
                                               margin_distributions,
                                               decoupler,
                                               test_matrix,
                                               support=target_q_support)
  ret <- list(
    posterior = posterior,
    warning = warning,
    flattened_errors=decouple_flattened,
    errors=decouple_values,
    decoupler_copula=decoupler_copula,
    decoupler_margins=margin_distributions
  )
  if (exists("rejection_results")) {
    ret$accepted_experts <- rejection_results$accepted_experts
  }
  ret
}

posterior_product_predictor <- function(test_matrix,
                                     overshoot = 0.1,
                                     k_percentiles=c(5,50,95),
                                        q_support_restriction=NULL) {

  decoupler <- get_CDF_decoupler(scale="linear", overshoot=overshoot, k_percentiles=k_percentiles, support_restriction = q_support_restriction)
  decoupler_fix_m <- decoupler$fix_m(test_matrix)

  logDM <- function(q) {
    pdf_vals <- decoupler_fix_m$f(q) |> abind::adrop(drop=3) # nxEx1 to nxE
    matrixStats::rowProds(pdf_vals)
  }

  list(
    posterior = list(logDM=logDM, support=decoupler_fix_m$support)
  )
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
                                        q_support_restriction=NULL) {
  if (is.null(error_metric)) {
    error_metric <- get_ratio_error_metric()
  }
  checkmate::assert_subset(q_support_restriction, c("non_negative", "strict_positive", NULL))
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

  target_q_support <- widen_support(target_q_support, 0.1, support_restriction = q_support_restriction)

  target_error_supports <- target_q_support_to_error_supports(target_q_support, test_matrix, error_metric)

  # inject here allows us to pass the arg=value pairs of the list as arguments
  margin_distributions <- rlang::inject(estimate_margins(flattened_errors, target_error_supports, !!!error_estimation_settings))
  margin_distributions <- add_d_e_to_list(margin_distributions, D)


  cdf_values <- purrr::imap(margin_distributions, \(margin, i) margin$cdf(flattened_errors[,i])) |> do.call(what=cbind)
  error_copula <- rlang::inject(fit_copula(cdf_values, copula_model, !!!vine_fit_settings))

  posterior<-create_log_unnormalized_posterior(error_copula,
                                    margin_distributions,
                                    error_metric,
                                    test_matrix,
                                    target_q_support)
  list(
    posterior = posterior,
    warning = warning,
    flattened_errors=flattened_errors,
    errors=errors,
    error_copula=error_copula,
    error_margins=margin_distributions
  )
}

decoupler_support_to_q_support <- function(decouple_support, m, e, d, decoupler) {
  # Calculate the q support from the error support
  q_support <- decoupler$f_inverse(decouple_support, m, e, d)
  if (!decoupler$f_increasing(m)[d,e]) {
    q_support <- c(q_support[2], q_support[1])
  }
  q_support
}

construct_posterior_margins <- function(margin_distributions, decoupler, m_matrix) {
  decoupler <- decoupler$fix_m(m_matrix)

  force(decoupler)
  single_output_f <- function(q, e, d) {
      decoupler$f(q)[e,d]
  }

  single_output_f_prime <- function(q, e, d) {
    decoupler$f_prime(q)[e,d]
  }

  is_increasing <- decoupler$is_increasing()

  margin_distributions |> purrr::map(
    \(dist) {
      e = dist$expert_id
      d = dist$d
      # create a new distribution with the decoupling function
      pdf <- function(q) {
        pdf_val <- dist$pdf(single_output_f(q, e, d))
        pdf_val * abs(single_output_f_prime(q, e, d))
      }
      cdf <- function(q) {
        cdf_val <- dist$cdf(single_output_f(q, e, d))
        cdf_val <- if (is_increasing[e, d]) cdf_val else 1 - cdf_val
        cdf_val
      }
      support <- decoupler_support_to_q_support(dist$support, m_matrix, e, d, decoupler)
      list(pdf = pdf,
           cdf = cdf,
           support = support)
    }
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

  m_realizations_test <- summarizing_function$f(test_set[k_percentiles_to_colname(k_percentiles)])
  test_set <- add_0_and_100_percentiles(test_set, k_percentiles)
  distributions <- interpolate_distributions(test_set, interpolation = interpolation)

  cdf_values <- purrr::imap(distributions, \(margin, i) margin$cdf(flattened_errors[,i])) |> do.call(what=cbind)


  error_copula <- fit_copula(cdf_values, copula_model)
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
  cal_scores <- training_set |> dplyr::group_by(expert_id) |> dplyr::summarise(calibration_score = {
    group_data = assesment_data[cur_group_rows(), ]
    realisations = group_data$realization
    calculateCalibrationScoreForExpert(group_data, realisations)
  })
  # Reject experts with low calibration score
  accepted_experts <- cal_scores |> dplyr::filter(calibration_score >= alpha) |> dplyr::pull(expert_id)
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

split_training_array_by_expert <- function(training_array) {
  checkmate::assert_array(training_array, "numeric", d=3)
  E = dim(training_array)[2]
  training_list <- purrr::map(seq_len(E), function(e) {
    training_array[,e,, drop=TRUE] # returns NxD matrices
  })
}

split_test_matrix_by_expert <- function(test_matrix, row_matrix_out=TRUE) {
  checkmate::assert_matrix(test_matrix, "numeric")
  E = nrow(test_matrix)
  test_list <- purrr::map(seq_len(E), function(e) {
    if (row_matrix_out) {
      test_matrix[e,, drop=FALSE] # returns 1xD matrices
    } else {
      test_matrix[e, drop=TRUE] # returns D vectors

    }
  })
  test_list
}

global_opt_weight_predict <- function(training_array, realizations, test_set_matrix, overshoot=0.1, support_restriction=NULL) {
  # Only works for data with quantiles at 5%, 50% and 90%
  cdf_values <- c(0.05, 0.50, 0.95)
  nr_experts = nrow(test_set_matrix)
  training_list <- split_training_array_by_expert(training_array)
  #training_assessments <- study_data_to_assessment_matrices(training_set, percentiles)
  #realisations <- training_set |> dplyr::distinct(question_id, realization) |> arrange(question_id) |> pull(realization)
  optimal_weights <- perfWeights_opt(training_list, realizations, cdf_values)

  quantiles <- add_0_and_100_percentiles_matrix(test_set_matrix, overshoot=overshoot, support_restriction=support_restriction)
  dists <- linear_distribution_interpolation_matrix(quantiles, c(0,cdf_values,1)) # E long
  support <- range(quantiles)
  # legacy code to get median
  DM <- constructDM(
    split_test_matrix_by_expert(test_set_matrix),
    optimal_weights,
    NULL,
    support[1],
    support[2],
    cdf_values
  )
  median <- DM[1, 2] # single row, second entry is the median

  DM_list <- (linear_combination_predict(dists, support, optimal_weights))
  DM_list$median = median
  DM_list
}

linear_combination_predict <- function(dists, support, weights=NULL) {
  if (is.null(weights)) {
    weights <- rep(1/length(dists), length(dists))
  }

  combined_mean = purrr::map2_dbl(dists, weights, \(dist, w) {
    dist$mean * w
  }) |> sum()

  logDM <- function(q) {
    pdf_vals <- purrr::map2(dists, weights, \(dist, w) {
      dist$pdf(q) * w
    }) |> do.call(what=cbind) # ends up with NxE matrix
    return(log(matrixStats::rowSums2(pdf_vals)))
  }

  return(list(logDM=logDM, support=support, mean=combined_mean))

}

equal_weight_predict <- function(test_matrix, accumulated_probabilities = c(0.05,0.5,0.95), overshoot=0.1, support_restriction=NULL) {
  checkmate::assert_matrix(test_matrix, "numeric")
  nr_experts <- nrow(test_matrix)
  quantiles <- add_0_and_100_percentiles_matrix(test_matrix, overshoot=overshoot, support_restriction=support_restriction)
  dists <- linear_distribution_interpolation_matrix(quantiles, c(0,accumulated_probabilities,1)) # E long
  support = range(quantiles)

  DM_list <- linear_combination_predict(dists, support)
  DM <- constructDM(
    split_test_matrix_by_expert(test_matrix),
    rep(1/nr_experts, nr_experts),
    NULL,
    support[1],
    support[2],
    accumulated_probabilities
  )
  median <- DM[1, 2] # single row, second entry is the median
  DM_list$median = median
  DM_list
}

equal_weight_predict_2 <- function(test_set) {
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
    posterior = list(),
    accepted_experts = vector(mode = "numeric"),
    experts = vector(mode = "numeric")
  )
  warnings <- c()
  single_expert_warning <- FALSE
  below_alpha_warning <- FALSE
  total_nr_experts <- length(unique(study_data$expert_id))

  for (i in seq(nrow(fold_combinations))) {
    test_set = fold_combinations$test[[i]]
    training_set = fold_combinations$training[[i]]
    arr_format <- df_format_to_array_format(training_set, test_set, p$summarizing_function$f, p$k_percentiles)
    test_matrix <- arr_format$test_summaries
    training_array <- arr_format$training_summaries
    training_realizations <- arr_format$training_realizations

    posterior = NA
    prediction = NA
    nr_experts <- total_nr_experts
    if (p$prediction_method == "copula") {
      if (is(p$error_metric, "decoupler")) {
        res <- tryCatch( {
         fit_and_construct_posterior_indep(
           training_array,
           training_realizations,
          test_matrix,
          p$error_metric,
          p$copula_model,
          vine_fit_settings = p$vine_fit_settings,
          error_estimation_settings = p$error_estimation_settings,
          q_support_restriction = p$q_support_restriction,
          q_support_overshoot = p$q_support_overshoot,
          rejection_test = p$rejection_test,
          rejection_threshold = p$rejection_threshold,
          rejection_min_experts = p$rejection_min_experts,
          connection_threshold = p$connection_threshold,
          recover_numerical_failure = p$recover_numerical_failure,
        )
        }, error = \(e) {
          warning(e)
          list(
          posterior = NULL,
          warning = c("ERROR_IN_FIT", e$message)
          )
        }
        )
      } else {
        res <- tryCatch( {
        fit_and_construct_posterior(
          arr_format$training_summaries,
          arr_format$training_realizations,
          test_matrix,
          p$copula_model,
          p$error_metric,
          vine_fit_settings = p$vine_fit_settings,
          error_estimation_settings = p$error_estimation_settings,
          q_support_restriction = p$q_support_restriction
        )
        }, error = \(e) {
          warning(e)
          list(
            posterior = NULL,
            warning = c("ERROR_IN_FIT", e$message)
          )
        }
        )
      }
      posterior <- res$posterior
      warnings <- append(warnings, res$warning)
      nr_experts <- dim(res$errors)[2]
    } else if (p$prediction_method == "density_product") {
      res <- posterior_product_predictor(
        df_format_to_array_format(training_set, test_set, p$summarizing_function$f, p$k_percentiles)$test_summaries,
        overshoot = p$overshoot,
        k_percentiles = p$k_percentiles,
        q_support_restriction = p$q_support_restriction
      )
      posterior <- res$posterior
    } else if (p$prediction_method == "perfect_expert") {
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
    } else if (p$prediction_method == "equal_weights") {
      posterior <- equal_weight_predict(test_matrix)
    } else if (p$prediction_method == "classical_global_opt") {
      posterior <- global_opt_weight_predict(training_array, training_realizations, test_matrix, p$overshoot, p$support_restriction)
    } else if (p$prediction_method == "uniform") {
      posterior <- uniform_prediction_posterior(test_matrix, p$overshoot, p$support_restriction)
    }
    else {
      stop(glue::glue("Unknown prediction method: {p$prediction_method}"))
    }

    realization <- unique(test_set$realization)
    test_question <- unique(test_set$question_id)

    stats$prediction[[i]] <- prediction
    stats$realization[[i]] <- realization
    stats$test_question[[i]] <- test_question
    stats$posterior[[i]] <- posterior
    stats$accepted_experts[[i]] <- nr_experts
    stats$experts[[i]] <- total_nr_experts
  }

  # convert to tibble
  stats <- tibble::as_tibble(stats)
  # return stats and warnings
  list(stats = stats, warnings = unique(warnings))
}

new_simulation_params <- function(copula_model,
                                      k_percentiles,
                                      interpolation,
                                      prediction_method,
                                      rejection_test = "distance_correlation",
                                      rejection_min_experts=1,
                                      rejection_threshold,
                                      summarizing_function,
                                      error_metric = NULL,
                                      error_estimation_settings = list(),
                                      vine_fit_settings=list(),
                                      q_support_restriction=NULL) {
  structure(list(
    copula_model = copula_model,
    k_percentiles = k_percentiles,
    interpolation = interpolation,
    prediction_method = prediction_method,
    rejection_test = rejection_test,
    rejection_threshold = rejection_threshold,
    rejection_min_experts=rejection_min_experts,
    summarizing_function = summarizing_function,
    error_metric = error_metric,
    vine_fit_settings = vine_fit_settings,
    error_estimation_settings=error_estimation_settings,
    q_support_restriction=q_support_restriction
  ), class="simulation_parameters")
}

default_simulation_params <- function(copula_model = "joe",
                                      k_percentiles = c(5, 50, 95),
                                      interpolation = "linear",
                                      prediction_method = "perfect_expert",
                                      rejection_test = "distance_correlation",
                                      rejection_min_experts = 1,
                                      rejection_threshold = 0.05,
                                      summarizing_function = get_three_quantiles_summarizing_function(),
                                      error_metric = get_CDF_decoupler(),
                                      error_estimation_settings = list(),
                                      vine_fit_settings=list(),
                                      q_support_restriction=NULL,
                                      q_support_overshoot=0.1,
                                      connection_threshold=NULL,
                                      connection_metric="kendall",
                                      recover_numerical_failure=FALSE) {
  checkmate::assert_string(copula_model)
  checkmate::assert_numeric(k_percentiles, sorted=TRUE, lower=0, upper=100)
  checkmate::assert_string(interpolation)
  checkmate::assert_string(prediction_method)
  checkmate::assert_subset(rejection_test, c("kruskal", "classical_calibration", "distance_correlation"))
  checkmate::assert_count(rejection_min_experts)
  checkmate::assert_number(rejection_threshold, null.ok=TRUE, lower=0, upper=1) # if NULL then no rejection happens
  checkmate::assert_class(summarizing_function, "summarizing_function", null.ok=TRUE)
  checkmate::assert(
    checkmate::check_class(error_metric, "decoupler"),
    checkmate::check_class(error_metric, "error_metric")
  )
  checkmate::assert_list(error_estimation_settings)
  checkmate::assert_list(vine_fit_settings)
  checkmate::assert_string(q_support_restriction, null.ok=TRUE)
  checkmate::assert_number(q_support_overshoot, lower=0, finite=TRUE)
  checkmate::assert_number(connection_threshold, null.ok=TRUE, lower=0, upper=1)


  params <- list(
    copula_model = copula_model,
    k_percentiles = k_percentiles,
    interpolation = interpolation,
    prediction_method = prediction_method,
    rejection_test = rejection_test,
    rejection_threshold = rejection_threshold,
    rejection_min_experts=rejection_min_experts,
    summarizing_function = summarizing_function,
    error_metric = error_metric,
    vine_fit_settings = vine_fit_settings,
    error_estimation_settings=error_estimation_settings,
    q_support_restriction=q_support_restriction,
    q_support_overshoot=q_support_overshoot,
    connection_threshold=connection_threshold,
    connection_metric=connection_metric,
    recover_numerical_failure=recover_numerical_failure
  )
  class(params) <- "simulation_parameters"
  params
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

  p + ggplot2::geom_text(data=support_df, aes(x=left_support, y=as.factor(expert_id), label="[", color=d, group="Support")) +
    ggplot2::geom_text(data=support_df, aes(x=right_support, y=as.factor(expert_id), label="]", color=d, group="Support")) +
    ggplot2::labs(group="al")
}


#' Title
#'
#' @param errors NxExD array with the errors
#'
#' @returns
#' @export
#'
#' @examples
plot_3d_errors <- function(errors, x_lab="Decoupler Value") {
  dims <- dim(errors)
  q_names <- dimnames(errors)[[1]]
  q_names <- if (is.null(q_names)) seq(dims[1]) else q_names
  e_names <- dimnames(errors)[[2]]
  e_names <- if (is.null(e_names)) seq(dims[2]) else e_names
  d_names <- dimnames(errors)[[3]]
  d_names <- if (is.null(d_names)) seq(dims[3]) else d_names
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

  #print(glue::glue("Plotting for questions {paste0(q_names, collapse=', ')}"))
  p<-df |> ggplot2::ggplot(ggplot2::aes(x=error, y=e, color=d)) +
    ggplot2::geom_point(alpha=0.9) + ggplot2::labs(x=x_lab, y="Expert ID", color="Error dimension")
  p
}

## make function to plot list of distributions
plot_distributions <- function(distributions, which_plots = "pdf&cdf", q_support_restriction=NULL, D=NULL, overshoot=0.1, fix_range=NULL, pdf_element_name="pdf", cdf_element_name="cdf") {
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
      new_support <- widen_support(distributions[[i]]$support, overshoot = overshoot, support_restriction = q_support_restriction)
    } else{
      new_support <- fix_range
    }
    x_data <- seq(new_support[1], new_support[2], length.out = 1000)
    y_cdf <- distributions[[i]][[cdf_element_name]](x_data)
    y_pdf <- distributions[[i]][[pdf_element_name]](x_data)
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
