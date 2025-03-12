
create_copula_posterior_unnormalized <- function(error_copula, expert_distributions) {
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

create_log_copula_posterior_unnormalized <- function(error_copula, error_distributions, error_metric) {
  stop("Not implemented yet")
  logDM <- function(q_vec) {
    cdf_values <- sapply(error_distributions, function(p) {
      p$cdf(q_vec)
    })

    log_pdf_values <- lapply(error_distributions, function(p) {
      log(p$pdf(q_vec))
    })
    added_log_pdf_values <- Reduce(`+`, log_pdf_values )

    return(log(copula::dCopula(survival_func_values, error_copula)) + added_log_pdf_values )
  }

  support = Reduce(function(support1, support2) {
    c(max(support1[1], support2[1]), min(support1[2], support2[2]))
  }, lapply(error_distributions, function(p) {
    p$support
  }))

  if (support[1] >= support[2]) {
    warning("Support is 0 wide")
  }

  return(list(logDM=logDM, support=support))

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

calc_sum_log_pdf_values <- function(expert_distributions, q_vec, d) {
  log_pdf_values <- purrr::map(expert_distributions, function(p) {
    log(p$pdf(q_vec))
  }) |> do.call(what = cbind) # is now a nxE matrix

  # normally sum over d but they are constant wrt d so we can multiply instead
  added_log_pdf_values = rowSums(log_pdf_values) * d
  added_log_pdf_values
}

calc_cdf_values <- function(expert_distributions, q_vec, is_increasing_flat, d) {
  n = length(q_vec)
  E = length(expert_distributions)
  stopifnot(length(is_increasing_flat) == d*E)

  cdf_values <- purrr::map(expert_distributions, function(p) {
    p$cdf(q_vec)
  }) |> do.call(what = cbind)
  # is now a nxE data.frame

  # duplicates each expert's cdf values d times.
  cdf_values <- cdf_values[,rep(seq_len(E), each=d)] # nx(d*E) matrix
  is_increasing_matrix <- matrix(is_increasing_flat, nrow=n, ncol=d*E, byrow=TRUE)
  cdf_values[!is_increasing_matrix] <- 1 - cdf_values[!is_increasing_matrix]
  cdf_values
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

  support = Reduce(function(support1, support2) {
    c(max(support1[1], support2[1]), min(support1[2], support2[2]))
  }, lapply(expert_distributions, function(p) {
    p$support
  }))

  if (support[1] >= support[2]) {
    warning("Support is 0 wide")
  }

  logDM <- function(q_vec) { # q_vec is n long
    stopifnot(is.vector(q_vec))

    # For each expert we get n cdf values
    cdf_values <- calc_cdf_values(expert_distributions, q_vec, is_increasing_flat, d)

    added_log_pdf_values <- calc_sum_log_pdf_values(expert_distributions, q_vec, d)

    added_log_error_metric_values <- calc_sum_log_error_terms_JC(error_metric, m, q_vec)

    if (is_vine) {
      copula_density_values <- rvinecopulib::dvinecop(cdf_values, error_copula, cores=1)
      # Not an obvious performance bootst on my system using more cores here.
      #copula_density_values <- rvinecopulib::dvinecop(cdf_values, error_copula, cores=4)
    } else {
      copula_density_values <- copula::dCopula(cdf_values, error_copula)
    }

    return(log(copula_density_values) + added_log_pdf_values + added_log_error_metric_values )
  }



  return(list(logDM=logDM, support=support))
}





## According to jouini and clemen (1996) equation 27.
create_copula_posterior_numerical_integration <- function(copula, expert_distributions) {
   res <- create_copula_posterior_unnormalized(copula, expert_distributions)

   area = integrate(res$DM, res$support[1], res$support[2], subdivisions = 10000)

   DM_normalized <- function(q_vec) {
     res$DM(q_vec) / area$value
   }

  return(list(DM=DM_normalized, support=res$support))
}
