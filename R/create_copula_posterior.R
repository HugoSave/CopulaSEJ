
## According to jouini and clemen (1996) equation 27.
create_copula_posterior_unnormalized_simple <- function(error_copula, expert_distributions) {
  DM <- function(q_vec) {
    survival_func_values <- sapply(expert_distributions, function(p) {
      # 1 - p$cdf(q_vec)
       p$cdf(q_vec)
    })

    pdf_values <- lapply(expert_distributions, function(p) {
      p$pdf(q_vec)
    })
    multiplied_pdf_values <- Reduce(`*`, pdf_values)

    return(copula::dCopula(survival_func_values, error_copula) * multiplied_pdf_values)
  }

  support = Reduce(function(support1, support2) {
    c(max(support1[1], support2[1]), min(support1[2], support2[2]))
  }, lapply(expert_distributions, function(p) {
    p$support
  }))

  if (support[1] >= support[2]) {
    warning("Support is 0")
  }

  return(list(DM=DM, support=support))
}

calc_sum_log_e_prime_m <- function(error_metric, m, q_vec) {
  error_to_3d_array(error_metric$f_prime_m, m, q_vec) |>
    abs()|> log() |> rowSums()
}

calc_sum_log_error_terms_JC <- function(error_metric, m, q_vec) {
  E = nrow(m)
  added_log_error_metric_values <- purrr::map_dbl(q_vec, function(q_i) {

    e_prime_m_values <- error_metric$f_prime_m(m, rep(q_i, E))
    e_prime_inverse_values <- error_metric$f_prime_inverse_q(m, e_prime_m_values)
    sum(log(abs(e_prime_m_values * e_prime_inverse_values)))
  })  # n long vector
  added_log_error_metric_values
}

calc_sum_log_pdf_values_JC <- function(expert_distributions, q_vec, d) {
  log_pdf_values <- purrr::map(expert_distributions, function(p) {
    log(p$pdf(q_vec))
  }) |> do.call(what = cbind) # is now a nxE matrix

  # normally sum over d but they are constant wrt d so we can multiply instead
  added_log_pdf_values = rowSums(log_pdf_values) * d
  added_log_pdf_values
}

flip_cdf_values_if_decreasing <- function(cdf_matrix, is_increasing_flat) {
  checkmate::assert_matrix(cdf_matrix)
  if (length(cdf_matrix) == 0)
    return(cdf_matrix)
  n = nrow(cdf_matrix)
  width = ncol(cdf_matrix)
  checkmate::assert_logical(is_increasing_flat, any.missing=FALSE, len = width)
  is_increasing_matrix <- matrix(is_increasing_flat, nrow=n, ncol=width, byrow=TRUE)
  cdf_matrix[!is_increasing_matrix] <- 1 - cdf_matrix[!is_increasing_matrix]
  cdf_matrix
}

calc_cdf_values_JC <- function(expert_distributions, q_vec, is_increasing_flat, d) {
  n = length(q_vec)
  E = length(expert_distributions)
  stopifnot(length(is_increasing_flat) == d*E)

  cdf_values <- purrr::map(expert_distributions, function(p) {
    p$cdf(q_vec)
  }) |> do.call(what = cbind)
  # is now a nxE data.frame

  # duplicates each expert's cdf values d times.
  cdf_values <- cdf_values[,rep(seq_len(E), each=d), drop=FALSE] # nx(d*E) matrix
  return(flip_cdf_values_if_decreasing(cdf_values, is_increasing_flat))
}

#' Title
#'
#' @param eval_points nxExd matrix
#' @param distributions list of E*d list elements with a cdf and pdf function
#'
#' @returns list with two nx(E*d) matrix cdf and pdf values
#'
#' @examples
calc_pdf_and_cdf_values_flat <- function(eval_points, distributions) {
  stopifnot(is.array(eval_points))
  dims <- dim(eval_points)
  n=dims[1]
  E = dims[2]
  D = dims[3]
  stopifnot(D*E==length(distributions))
  cdf_values <- array(NA, dim=c(n, D, E)) # not conventional dimension order (nxExd) to make flattening easier
  pdf_values <- array(NA, dim=c(n, D, E))
  for (i in seq_len(D*E)) {
    indices <- linear_index_to_d_E(i, D)
    e = indices$e
    d = indices$d
    evaluation <- eval_points[,e,d] # n vector
    cdf_values[,d,e] <- distributions[[i]]$cdf(evaluation)
    pdf_values[,d,e] <- distributions[[i]]$pdf(evaluation)
  }
  # n x (d * E) matrix
  flat_cdf_values <- array(cdf_values, dim=c(n, d*E))
  flat_pdf_values <- array(pdf_values, dim=c(n, d*E))

  # flat_cdf_values <- flip_cdf_values_if_decreasing(flat_cdf_values, is_increasing_flat)

  list(cdf_values=flat_cdf_values, pdf_values=flat_pdf_values)
}


calc_is_increasing_flat <- function(error_metric, m) {
  is_increasing = error_metric$f_increasing_q(m) # E x d matrix
  is_increasing_flat = flatten_matrix_row_by_row(is_increasing) # E*d long vector
  is_increasing_flat
}

#' Returns the log of the unnormalized posterior distribution of the posterior
#' where the JC assumption of perfect experts is in place.
#'
#' @param error_copula class "Copula" with dimension d*E. The copula should have its
#' arguments ordered such that the the first d arguments belongs to expert 1 and
#' so on.
#' @param expert_distributions list of expert distributions. Length E.
#' @param error_metric
#' @param m Exd matrix of the summaries properties
#'
#' @export
#'
create_log_unnormalized_posterior_JC <- function(error_copula, expert_distributions, error_metric, m) {
  d = ncol(m)
  E = nrow(m)
  if (is(error_copula, "vinecop")) {
    c_dim = length(error_copula$var_types)
    is_vine = T
  } else if (is(error_copula, "Copula")) {
    c_dim = error_copula@dimension
    is_vine = F
  } else {
    stop("Unknown copula type")
  }
  stopifnot(c_dim==d*E)
  stopifnot(E==length(expert_distributions))
  is_increasing_flat = calc_is_increasing_flat(error_metric, m) # E*d long vector
  stopifnot(!any(is.na(is_increasing_flat)))

  support = support_intersection(expert_distributions, nested_list=TRUE)

  logDM <- function(q_vec) { # q_vec is n long
    stopifnot(is.vector(q_vec))

    # For each expert we get n cdf values
    cdf_values <- calc_cdf_values_JC(expert_distributions, q_vec, is_increasing_flat, d)

    added_log_pdf_values <- calc_sum_log_pdf_values_JC(expert_distributions, q_vec, d)

    added_log_error_metric_values <- calc_sum_log_error_terms_JC(error_metric, m, q_vec)

    log_copula_density_values <- error_copula$density(cdf_values, log=TRUE)

    return(log_copula_density_values + added_log_pdf_values + added_log_error_metric_values )
  }

  return(list(logDM=logDM, support=support))
}

create_log_unnormalized_posterior_indep <- function(indep_copula, indep_margins, indep_class, m_matrix, support=NULL) {
  checkmate::assert_matrix(m_matrix, "numeric", any.missing = FALSE)
  E = nrow(m_matrix)
  checkmate::assert_list(indep_copula)
  checkmate::assert_subset("density",names(indep_copula))
  nr_dims <- dim(indep_copula)[[1]]
  checkmate::assert_list(indep_margins, len=nr_dims)

  indep_fix_m <- indep_class$fix_m(m_matrix)

  #is_increasing_flat <- calc_is_increasing_flat(indep_copula, m)

  stopifnot(!is.null(support))
  # TODO fix support. Have to implement inverse of the indep_class function.
  #support <- indep_margins |> purrr::map(\(x) x$support) |> error_supports_to_q_supports(m, error_metric) |> support_intersection(nested_list=FALSE)
  #support <- calc_support_error_marginals(error_margins, error_metric, m) |> support_intersection(nested_list=FALSE)


  logDM <- function(q_input) {
    ret_density <- vector("double", length=length(q_input))
    out_of_support <- q_input <= support[1] | q_input >= support[2]
    ret_density[out_of_support] <- -Inf

    q_vec <-q_input[!out_of_support]
    if (length(q_vec) == 0) {
      return(ret_density)
    }

    n = length(q_vec)
    indep_values <- indep_fix_m$f(q_vec) # nxExD array
     #error_to_3d_array(error_metric$f, m, q_vec) |> # nxExd array
    pdf_cdf_vals <- calc_pdf_and_cdf_values_flat(indep_values, indep_margins)

    pdf_values <- pdf_cdf_vals$pdf_values
    cdf_values <- pdf_cdf_vals$cdf_values

    log_copula_density_values <- indep_copula$density(cdf_values, log=TRUE)
    # if (is(indep_copula, "vinecop_dist")) {
    #   copula_density_values <- rvinecopulib::dvinecop(cdf_values, indep_copula, cores=1)
    #   log_copula_density_values <- log(copula_density_values)
    # } else {
    #   log_copula_density_values <- copula::dCopula(cdf_values, indep_copula, log=TRUE)
    # }
    added_log_pdf_values <- rowSums(log(pdf_values))

    added_log_e_prime_m <- indep_fix_m$f_prime_q(q_vec) |> abs() |> log() |> rowSums()

    #added_log_e_prime_m <- calc_sum_log_e_prime_m(error_metric, m, q_vec)
    ret_density[!out_of_support] <- log_copula_density_values + added_log_pdf_values + added_log_e_prime_m
    ret_density
  }
  return(list(logDM=logDM, support=support))
}

#' Title
#'
#' @param error_copula Already fitted error copula
#' @param error_margins
#' @param error_metric
#' @param m Exd matrix of the summaries properties
#'
#' @returns
#' @export
#'
#' @examples
create_log_unnormalized_posterior <- function(error_copula, error_margins, error_metric, m, support=NULL) {
  checkmate::assert_matrix(m, "numeric", any.missing = FALSE)
  d = ncol(m)
  E = nrow(m)
  checkmate::assert_list(error_margins, len=d*E)

  checkmate::assert_set_equal(dim(error_copula)[[1]], d*E)
  is_increasing_flat <- calc_is_increasing_flat(error_metric, m)


  if (is.null(support)) {
    support <- error_margins |> purrr::map(\(x) x$support) |> error_supports_to_q_supports(m, error_metric) |> support_intersection(nested_list=FALSE)
  }

  logDM <- function(q_input) {
    ret_density <- vector("double", length=length(q_input))
    out_of_support <- q_input <= support[1] | q_input >= support[2]
    ret_density[out_of_support] <- -Inf

    q_vec <-q_input[!out_of_support]
    if (length(q_vec) == 0) {
      return(ret_density)
    }

    n = length(q_vec)
    pdf_cdf_vals <- error_to_3d_array(error_metric$f, m, q_vec) |> # nxExd array
      calc_pdf_and_cdf_values_flat(error_margins) # evaluate error pdfs and cdfs at the error values

    pdf_values <- pdf_cdf_vals$pdf_values
    cdf_values <- flip_cdf_values_if_decreasing(pdf_cdf_vals$cdf_values, is_increasing_flat)

    log_copula_density_values <- error_copula$density(cdf_values, log=TRUE)
    added_log_pdf_values <- rowSums(log(pdf_values))

    added_log_e_prime_m <- calc_sum_log_e_prime_m(error_metric, m, q_vec)
    ret_density[!out_of_support] <- log_copula_density_values + added_log_pdf_values + added_log_e_prime_m
    ret_density

  }
  return(list(logDM=logDM, support=support))
}

calc_normalization_constant <- function(unnorm_posterior, support, method="numerical") {
  if (method == "numerical") {
    area = integrate(unnorm_posterior, support[1], support[2], subdivisions = 10000)
    return(area$value)
  } else if (method == "uniform importance") {
    samples <- runif(10000, support[1], support[2])
    area <- mean(unnorm_posterior(samples))
    return(area)
  } else if (method=="bridge sampler") {
    stop("Not implemented yet ")
  }
  else {
    stop(glue::glue("Unknown normalization method: {method}"))
  }
}

# numerical integration of log posterior
calc_normalization_constant_log_posterior <- function(unnorm_log_posterior, support, method="numerical") {

  non_log_post <- function(q_vec) {
    exp(unnorm_log_posterior(q_vec))
  }

  calc_normalization_constant(non_log_post, support, method)
}



## According to jouini and clemen (1996) equation 27.
create_copula_posterior_numerical_integration <- function(copula, expert_distributions) {
   res <- create_copula_posterior_unnormalized_simple(copula, expert_distributions)

   area = integrate(res$DM, res$support[1], res$support[2], subdivisions = 10000)

   DM_normalized <- function(q_vec) {
     res$DM(q_vec) / area$value
   }

  return(list(DM=DM_normalized, support=res$support))
}
