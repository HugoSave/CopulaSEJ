


#' Constructor for error metric class.
#'
#' @param name name of the function
#' @param f_single,f_single_prime_m The error function,
#' its derivative wrt m. Both functions have three
#' arguments (m,q,d) where m and q are equal length vectors and d is a an
#' integer specifying the kind of
#' property m is. Value of d can thus change the error function used.
#' @param f,f_prime_m The error metric and the derivative of f with respect to
#' m. Both functions have two arguments (m,q) where m is a nxd matrix and
#' q is a vector of length n. Returns a matrix of the same shape as m.
#' @param f_single_inverse_q function of (m,e,d) where m and e are vectors
#' of same length and d is an integer. Outputs a vector of the same length as m.
#' corresponding to the q values that would have generated the e vector.
#' @param f_inverse_q,f_prime_inverse_q The inverse of f with respect to
#' q and the derivative of the inverse respectively. Have two arguments (m,e)
#' where m and e are matrices both with size nxd. Outputs a matrix of nxd.
#' @param f_single_increasing_q Function of (m, d) where m is a vector and d is
#' a scalar
#' @param f_increasing_q Function of a single argument (m) where m is a matrix
#' with d columns.
#' @param check_domain function of (m, q) that checks if the error metric is
#' well defined on the given values. The function returns a string if there is
#' an error or TRUE if the values are OK. If this function is not provided, it
#' defaults to function always returning TRUE.
#' @param d The dimension of m this error metric is applicable to. Defaults to
#' Inf which means it is applicable to any dimension.
#'
#' @returns A new error metric object
#' @export
#'
new_error_metric <- function(name,
                             f_single,
                             f,
                             f_single_prime_m,
                             f_prime_m,
                             f_prime_inverse_q,
                             f_single_increasing_q,
                             f_increasing_q,
                             f_single_inverse_q,
                             f_inverse_q = NULL,
                             check_domain = NULL,
                             d = Inf) {
  checkmate::assert_string(name)
  checkmate::assert_function(f_single)
  checkmate::assert_function(f_single_prime_m)
  checkmate::assert_function(f)
  checkmate::assert_function(f_prime_m)
  checkmate::assert_function(f_prime_inverse_q)
  checkmate::assert_function(f_single_increasing_q)
  checkmate::assert_function(f_increasing_q)
  checkmate::assert_function(f_single_inverse_q)
  checkmate::assert_function(f_inverse_q, null.ok = TRUE)
  checkmate::assert_function(check_domain, null.ok = TRUE)
  checkmate::assert(checkmate::check_count(d, positive = TRUE),
                    checkmate::check_number(d, lower = Inf))

  if (is.null(check_domain)) {
    check_domain = function(m, q) {TRUE}
  }

  return(structure(
    list(
      name = name,
      f_single = f_single,
      f = f,
      f_single_prime_m=f_single_prime_m,
      f_prime_m = f_prime_m,
      f_prime_inverse_q = f_prime_inverse_q,
      f_single_increasing_q = f_single_increasing_q,
      f_increasing_q = f_increasing_q,
      f_single_inverse_q=f_single_inverse_q,
      f_inverse_q = f_inverse_q,
      check_domain=check_domain,
      d = d
    ),
    class = "error_metric"
  ))
}

new_decoupler <- function(name, f, f_prime_q, f_increasing, fix_m, f_inverse=NULL) {
  checkmate::assert_string(name)
  checkmate::assert_function(f, args=c("q", "m"), ordered=TRUE)
  checkmate::assert_function(f_prime_q, args=c("q", "m"), ordered=TRUE)
  checkmate::assert_function(f_increasing, args=c("m"), ordered=TRUE)
  checkmate::assert_function(fix_m, args=c("m"), ordered=TRUE)
  checkmate::assert_function(f_inverse, args=c("z", "m"), null.ok=TRUE)

  return(structure(
    list(
      name = name,
      f = f,
      f_prime_q = f_prime_q,
      f_increasing=f_increasing,
      fix_m = fix_m,
      f_inverse=f_inverse
    ),
    class = "decoupler"
  ))
}

print.decoupler <- function(object, ...) {
  cat(object$name, "decoupler\n")
  print(glue::glue("Includes functions with output dimension nxExD:"))
  cat("f(q,m)\n")
  cat("f_prime_q(q,m)\n")
  cat("f_increasing(m)\n")
  cat("fix_m(m)\n")
}

print.fixed_decopuler <- function(object, ...) {
  cat(object$name, "fixed decoupler\n")
  print(glue::glue("Includes functions with output dimension nxEx{object$D}:"))
  cat("f(q)\n")
  cat("f_prime_q(q)\n")
  cat("f_increasing()\n")
}

new_fixed_decoupler <- function(name, f, f_prime_q, f_increasing, D, f_inverse=NULL) {
  checkmate::assert_string(name)
  checkmate::assert_function(f, args=c("q"), ordered=TRUE)
  checkmate::assert_function(f_prime_q, args=c("q"), ordered = TRUE)
  checkmate::assert_function(f_increasing, args=character(0))
  checkmate::assert_count(D, positive=TRUE)
  checkmate::assert_function(f_inverse, args=c("z"), ordered=TRUE, null.ok=TRUE)

  return(structure(
    list(
      name = name,
      f = f,
      f_prime_q = f_prime_q,
      f_increasing=f_increasing,
      D=D,
      f_inverse=f_inverse
    ),
    class = "fixed_decoupler"
  ))
}


validate_error_metric_on_observations <- function(error_metric, m_training, q_training, m_test) {

  is_inc <- error_metric$f_increasing_q(m_training)
  if (any(is.nan(is_inc))) {
    stop("Error: some of the assessments make the error metric are not strictly monotinc wrt q")
  }

  is_inc <- error_metric$f_increasing_q(m_test)
  if (any(is.nan(is_inc))) {
    stop("Error: some of the assessments make the error metric are not strictly monotinc wrt q")
  }
}

q_over_m_check_domain <- function(m, q, domain_limit_epsilon = 1e-3) {
  collection <- checkmate::makeAssertCollection()

  checkmate::assert(
    checkmate::check_array(m, "numeric"),
    checkmate::check_numeric(m),
    add=collection
  )
  checkmate::assert_numeric(q, add=collection)

  if (min(abs(m)) < domain_limit_epsilon) {
    collection$push(glue::glue("some of the assessments, m, are too close to zero for the ratio error metric. Must be further away than {domain_limit_epsilon}."))
  }
  if (min(abs(q)) < domain_limit_epsilon) {
    collection$push(glue::glue("some of the realizations, q, are too close to zero for the ratio error metric. Must be further away than {domain_limit_epsilon}."))
  }

  m <- as.matrix(m, nrow=length(m), ncol=1) # flatten but it should be a matrix to comply with function definitions

  is_inc <- as.vector(q_over_m_increasing_q(m))
  if (any(is.na(is_inc))) {
    specific_m_value <- m[which(is.na(is_inc))] |> head(1)
    collection$push(glue::glue("Error: some of the assessments, m, make the error metric undefined. For example m = {specific_m_value}"))
  }

  if (collection$isEmpty()) {
    return(TRUE)
  } else {
    return(collection$getMessages() |> paste0(collapse="\n"))
  }
}

ratio_error_check_domain <- function(m, q, domain_limit_epsilon = 1e-3) {
  collection <- checkmate::makeAssertCollection()

  checkmate::assert(
    checkmate::check_array(m, "numeric"),
    checkmate::check_numeric(m),
    add=collection
  )
  checkmate::assert_numeric(q, add=collection)

  if (min(abs(m)) < domain_limit_epsilon) {
    collection$push(glue::glue("some of the assessments, m, are too close to zero for the ratio error metric. Must be further away than {domain_limit_epsilon}."))
  }
  if (min(abs(q)) < domain_limit_epsilon) {
    collection$push(glue::glue("some of the realizations, q, are too close to zero for the ratio error metric. Must be further away than {domain_limit_epsilon}."))
  }

  m <- as.matrix(m, nrow=length(m), ncol=1) # flatten but it should be a matrix to comply with function definitions

  # all q must be either all positive or all negative
  checkmate::assert(
    checkmate::check_true(all(q > 0)),
    checkmate::check_true(all(q < 0)),
    add=collection
  )

  is_inc <- as.vector(ratio_error_increasing_q(m))

  if (any(is.na(is_inc))) {
    specific_m_value <- m[which(is.na(is_inc))] |> head(1)
    collection$push(glue::glue("Error: some of the assessments, m, make the error metric undefined. For example m = {specific_m_value}"))
  }

  if (collection$isEmpty()) {
    return(TRUE)
  } else {
    return(collection$getMessages() |> paste0(collapse="\n"))
  }
}

get_valid_rows_error_metric <- function(m_values, q, error_metric) {
  # exclude questions with no assessments
  e <- error_metric(m_values, q)
  # return row indicies with finite values
  e_finite = is.finite(e)
  return(which(apply(e_finite, 1, all)))
}

linear_error_single <- function(m, q, d) {
  stopifnot(is.vector(q))
  stopifnot(is.vector(m))
  return(q - m)
}

linear_error <- function(m, q) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  return(q - m)
}

linear_error_prime_m <- function(m, q) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  m[] = -1
  return(m)
}

linear_error_single_prime_m <- function(m, q, d) {
  stopifnot(is.vector(q))
  stopifnot(is.vector(m))
  m[] = -1
  return(m)
}

linear_error_single_inverse_q <- function(m, epsilon, d) {
  # convert m and epsilon to matrices in case they are provided as vectors
  stopifnot(is.vector(epsilon), is.vector(m))
  return(m + epsilon)
}


linear_error_inverse_q <- function(m, epsilon) {
  # convert m and epsilon to matrices in case they are provided as vectors
  argcheck <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- argcheck[[1]]
  epsilon <- argcheck[[2]]

  return(m + epsilon)
}


linear_error_prime_inverse_q <- function(m, epsilon) {
  # TODO e is not a vector but could be a matrix.
  argcheck <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- argcheck[[1]]
  epsilon <- argcheck[[2]]

  return(matrix(1, nrow(m), ncol(m)))

  stopifnot(is.vector(e), is.vector(m))
  return(matrix(1, length(e), length(m)))
}

linear_error_prime_inverse_q_2 <- function(m, epsilon) {
  # TODO e is not a vector but could be a matrix.
  m <- typecheck_and_convert_matrix_vector(m, vector())
  epsilon <- typecheck_and_convert_matrix_vector(epsilon, vector())
  return(matrix(1, ncol(m), length(m)))
}


linear_error_single_increasing_q <- function(m, d) {
  stopifnot(is.vector(m))
  return(rep.int(TRUE, length(m)))
}

linear_error_increasing_q <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(matrix(TRUE, nrow(m), ncol(m)))
}


# linear error metric
get_linear_error_metric <- function() {
  return(
    new_error_metric(
      "linear error q-m",
      linear_error_single,
      linear_error,
      linear_error_single_prime_m,
      linear_error_prime_m,
      linear_error_prime_inverse_q,
      linear_error_single_increasing_q,
      linear_error_increasing_q,
      linear_error_single_inverse_q,
      linear_error_inverse_q
    )
  )
}

is_positive_no_zero <- function(x) {
  if (is.vector(x)) {
    x_ret <- vector("logical", length(x))
  } else if (is.array(x)) {
    x_ret <- array(NA, dim(x), dimnames = dimnames(x))
  }


  na_values = which(x == 0)
  true_values = which(x > 0)
  false_values = which(x < 0)
  x_ret[true_values] = TRUE
  x_ret[false_values] = FALSE
  x_ret[na_values] = NA
  x_ret
}

# checks if first argument is a vector or matrix and if second argument is a
# vector. Converts first argument to a matrix.
typecheck_and_convert_matrix_vector <- function(m, q) {
  stopifnot(is.vector(q))
  if (is.data.frame(m)) {
    m <- as.matrix(m)
  }
  stopifnot(is.vector(m) || is.matrix(m))
  if (is.vector(m)) {
    # make it a single column matrix
    m <- matrix(m, ncol = 1)
  }
  assertthat::are_equal(nrow(m), length(q))
  m
}

typecheck_and_convert_matrix_matrix <- function(m, e) {
  if (is.data.frame(m)) {
    m <- as.matrix(m)
  }
  if (is.data.frame(e)) {
    e <- as.matrix(e)
  }

  stopifnot(is.vector(e) || is.matrix(e))
  stopifnot(is.vector(m) || is.matrix(m))
  if (is.vector(m)) {
    # make it a single column matrix
    m <- matrix(m, ncol = 1)
  }
  if (is.vector(e)) {
    e <- matrix(e, ncol = 1)
  }
  stopifnot(all(dim(m) == dim(e)))
  list(m, e)
}

q_over_m_single <- function(m, q, d) {
  stopifnot(is.vector(q))
  stopifnot(is.vector(m))
  return(q/m)
}

q_over_m <-function(m, q) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  # repeats q over each column of m to end up with the same size.
  q_rep_matrix = outer(q, rep.int(1L, ncol(m)))
  return(q_rep_matrix / m)
}

q_over_m_prime_m <- function(m, q) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  q_rep_matrix = outer(q, rep.int(1L, ncol(m)))
  return(-q_rep_matrix / (m^2))
}

q_over_m_single_prime_m <- function(m, q, d) {
  stopifnot(is.vector(q))
  stopifnot(is.vector(m))
  return(- q / m^2)
}

q_over_m_single_inverse_q <- function(m, epsilon, d) {
  stopifnot(is.vector(epsilon), is.vector(m))
  return(m * epsilon)
}

q_over_m_inverse_q <- function(m, epsilon) {
  args_checked <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- args_checked[[1]]
  epsilon <- args_checked[[2]]
  return(m * epsilon)

}

q_over_m_prime_inverse_q <- function(m, epsilon) {
  args_checked <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- args_checked[[1]]
  epsilon <- args_checked[[2]]
  return(m)
}

q_over_m_single_increasing_q <- function(m, d) {
  stopifnot(is.vector(m))
  return(is_positive_no_zero(m))
}

q_over_m_increasing_q <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(is_positive_no_zero(m))
}

get_q_over_m_error_metric <- function() {
  return(
    new_error_metric(
      "ratio error q/m",
      q_over_m_single,
      q_over_m,
      q_over_m_single_prime_m,
      q_over_m_prime_m,
      q_over_m_prime_inverse_q,
      q_over_m_single_increasing_q,
      q_over_m_increasing_q,
      q_over_m_single_inverse_q,
      q_over_m_inverse_q,
      q_over_m_check_domain
    )
  )
}



############################
# Implementation of m/q

ratio_error_single <- function(m, q, d) {
  stopifnot(is.vector(q))
  stopifnot(is.vector(m))
  return(m / q)
}

ratio_error <- function(m, q) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  # repeats q over each column of m to end up with the same size.
  q_rep_matrix = outer(q, rep.int(1L, ncol(m)))
  return(m / q_rep_matrix)
}

ratio_error_prime_m <- function(m, q) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  q_rep_matrix = outer(q, rep.int(1L, ncol(m)))
  return(1 / q_rep_matrix)
}

ratio_error_single_prime_m <- function(m, q, d) {
  stopifnot(is.vector(q))
  stopifnot(is.vector(m))
  return(1 / q)
}

ratio_error_single_inverse_q <- function(m, epsilon, d) {
  stopifnot(is.vector(epsilon), is.vector(m))
  return(m / epsilon)
}

# epsilon is a vector with length L
# m is L long
# output is L
ratio_error_inverse_q <- function(m, epsilon) {
  args_checked <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- args_checked[[1]]
  epsilon <- args_checked[[2]]
  return(m / epsilon)

}

ratio_error_prime_inverse_q <- function(m, epsilon) {
  args_checked <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- args_checked[[1]]
  epsilon <- args_checked[[2]]

  return(-m / (epsilon^2))

  # - m/e^2
  # repeat the epsilon over multiple columns
  epsilon_rep_matrix = outer(epsilon, rep.int(1L, length(m)))
  # repeat m for each row
  m_rep_matrix = outer(rep.int(1L, length(epsilon)), m)
  return(-m_rep_matrix / (epsilon_rep_matrix^2))
}

ratio_error_single_increasing_q <- function(m, d) {
  stopifnot(is.vector(m))
  return(!is_positive_no_zero(m))
}

ratio_error_increasing_q <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(!is_positive_no_zero(m))
}

get_ratio_error_metric <- function() {
  return(
    new_error_metric(
      "ratio error m/q",
      ratio_error_single,
      ratio_error,
      ratio_error_single_prime_m,
      ratio_error_prime_m,
      ratio_error_prime_inverse_q,
      ratio_error_single_increasing_q,
      ratio_error_increasing_q,
      ratio_error_single_inverse_q,
      ratio_error_inverse_q,
      ratio_error_check_domain
    )
  )
}

############################
# Implementation of relative error (q-m)/q=1-m/q
relative_error_single <- function(m, q, d) {
  return(1 - ratio_error_single(m,q,d))
}
relative_error <- function(m, q) {
  return(1 - ratio_error(m,q))
}

relative_error_prime_m <- function(m, q) {
  return(-ratio_error_prime_m(m,q))
}

relative_error_single_prime_m <- function(m, q, d) {
  return(-ratio_error_single_prime_m(m,q,d))
}

relative_error_single_inverse_q <- function(m, epsilon, d) {
  stopifnot(is.vector(epsilon), is.vector(m))
  return(ratio_error_single_inverse_q(m, 1-epsilon, d))
}

relative_error_inverse_q <- function(m, epsilon) {
  return(ratio_error_inverse_q(m, 1-epsilon))
}

relative_error_prime_inverse_q <- function(m, epsilon) {
  return(-ratio_error_prime_inverse_q(m, 1-epsilon))
}

relative_error_single_increasing_q <- function(m, d) {
  stopifnot(is.vector(m))
  return(is_positive_no_zero(m))
}

relative_error_increasing_q <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(is_positive_no_zero(m))
}

get_relative_error_metric <- function() {
  return(
    new_error_metric(
      "relative error (1-m/q)",
      f_single=relative_error_single,
      f=relative_error,
      f_single_prime_m=relative_error_single_prime_m,
      f_prime_m=relative_error_prime_m,
      f_prime_inverse_q=relative_error_prime_inverse_q,
      f_single_increasing_q=relative_error_single_increasing_q,
      f_increasing_q=relative_error_increasing_q,
      f_single_inverse_q=relative_error_single_inverse_q,
      f_inverse_q=relative_error_inverse_q,
      check_domain = ratio_error_check_domain,
      d=Inf
    )
  )
}

not_yet_implemented <- function(...) {
  stop("Not yet implemented")
}

#### Linear decoupler
decoupler_linear <- function(q, m) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  q_rep_matrix = array(q, dim=(c(N, E, D)))
  q - m_rep_matrix
}

decoupler_linear_prime_q <- function(q, m) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  array(1, dim=(c(N, E, D)))
}

decoupler_linear_increasing <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(matrix(TRUE, nrow=nrow(m), ncol=ncol(m)))
}

decoupler_linear_inverse <- function(z, m) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  z_rep_matrix = array(z, dim=(c(N,E,D)))
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  (z_rep_matrix + m_rep_matrix)
}

decoupler_linear_fix_m <- function(m) {
  increasing_matrix = decoupler_linear_increasing(m)
  new_fixed_decoupler(
    name="linear",
    f=\(q) decoupler_linear(q,m),
    f_prime_q=\(q) decoupler_linear_prime_q(q,m),
    f_increasing=\() increasing_matrix,
    D=ncol(m),
    f_inverse=\(z) decoupler_linear_inverse(z,m)
  )
}

get_linear_decoupler <- function() {
  return(
    new_decoupler(
      name="linear decoupler",
      f=decoupler_linear,
      f_prime_q=decoupler_linear_prime_q,
      f_increasing=decoupler_linear_increasing,
      fix_m = decoupler_linear_fix_m,
      f_inverse=decoupler_linear_inverse
    )
  )
}


##### CDF error
indep_CDF <- function(q, m, overshoot=0.1, k_percentiles=c(5,50,95), scale="linear"){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles, scale)
  ind_cond_m$f(q)
}

indep_CDF_prime_q <- function(q, m, overshoot=0.1, k_percentiles=c(5,50,95), scale="linear"){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles, scale)
  ind_cond_m$f_prime_q(q)
}

indep_CDF_increasing <- function(m, overshoot=0.1, k_percentiles=c(5,50,95)){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles)
  ind_cond_m$f_increasing()
}

decouple_CDF_inverse <- function(z, m, overshoot=0.1, k_percentiles=c(5,50,95), scale="linear"){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles, scale)
  ind_cond_m$f_inverse(z)
}


get_CDF_decoupler <- function(scale="linear") {
  return(
    new_decoupler(
      name="CDF",
      f=\(q,m) indep_CDF(q, m, scale=scale),
      f_prime_q=\(q,m) indep_CDF_prime_q(q, m, scale=scale),
      f_increasing=indep_CDF_increasing,
      fix_m = \(m) get_cdf_indep_function_fix_m(m, scale=scale),
      f_inverse=\(z, m) decouple_CDF_inverse(z, m, scale=scale)
    )
  )
}

#' Title
#'
#' @param m A matrix ExD
#' @param overshoot A scalar
#' @param k_percentiles A vector of percentiles
#'
#' @returns
#' @export
#'
#' @examples
get_cdf_indep_function_fix_m <- function(m, overshoot=0.1, k_percentiles = c(5,50,95), scale="linear") {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  checkmate::assert_numeric(k_percentiles, len=ncol(m))
  N = nrow(m)
  is_increasing_matrix <- matrix(TRUE, nrow=N, ncol=1) # single CDF value per expert
  L = min(m)
  U = max(m)
  if (scale=="linear"){
    new_support <- widen_support(c(L, U), overshoot, support_restriction = NULL)
  } else if (scale=="log") {
    new_support <- widen_support(c(L, U), overshoot, support_restriction = "strict_positive")
  }
  L_star = new_support[1]
  U_star = new_support[2]
  cdf_values <- c(0, k_percentiles/100, 1)
  extended_m_matrix <- abind::abind(matrix(L_star, nrow=N, ncol=1), m, matrix(U_star, nrow=N, ncol=1), along=2)
  if (scale=="linear") {
    distributions <- purrr::array_branch(extended_m_matrix, 1) |>
      purrr::map(\(fractiles) {
        linear_distribution_interpolation(fractiles, cdf_values)
      })
  } else if(scale=="log") {
    distributions <- purrr::array_branch(extended_m_matrix, 1) |>
      purrr::map(\(fractiles) {
        log_linear_distribution_interpolation(fractiles, cdf_values)
      })
  } else {
    stop("Unknown scale")
  }

  # returns ZxNx1 if N is rows of m and Z is length of z
  f_inverse <- function(z) {
    values <- distributions |> purrr::map(\(p) {
      p$cdf_inv(z)
    }) |> do.call(what = cbind)
    dim(values) <- c(nrow(values), ncol(values), 1)
    values
  }

  f <- function(q) {
    values <- distributions |> purrr::map(\(p) {
      p$cdf(q)
    }) |> do.call(what = cbind)
    dim(values) <- c(nrow(values), ncol(values), 1)
    values
  }


  f_prime_q_fix_m <- function(q) {
    values <- distributions |> purrr::map(\(p) {
      p$pdf(q)
    }) |> do.call(what = cbind)
    # add a dimension 1 to values
    dim(values) <- c(nrow(values), ncol(values), 1)
    values

  }

  f_increasing <- function() {
    return(is_increasing_matrix)
  }

  new_fixed_decoupler(
    name = "CDF",
    f=f,
    f_prime_q=f_prime_q_fix_m,
    f_increasing = f_increasing,
    D=1,
    f_inverse=f_inverse
  )
}


###### Sigmoid functions to compose above error functions with

sigmoid_centered_image <- function(x, L, x_0, k) {
  return(L / (1 + exp(-k * (x - x_0))) - L / 2)
}

sigmoid_centered_prime <- function(x, L, x_0, k) {
  exp_val = exp(-k * (x - x_0))
  return(k * L * exp_val / (1 + exp_val)^2)
}

sigmoid_centered_image_inverse <- function(y, L, x_0, k) {
  return(x_0 - log(L / (y + L / 2) - 1) / k)
}

sigmoid_centered_image_inverse_prime <- function(y, L, x_0, k) {
  return(4*L/(k*(L - 2*y)*(L+2*y)))
}

get_sigmoid_q_over_m_error_metric <- function(k=1) {
  error_metric <- get_q_over_m_error_metric()
  L=2
  x_0 = 1
  compose_sigmoid_with_error(error_metric,
                             sigmoid=purrr::partial(sigmoid_centered_image, L=L, x_0=x_0, k=k),
                             sigmoid_prime=purrr::partial(sigmoid_centered_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse = purrr::partial(sigmoid_centered_image_inverse, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_prime = purrr::partial(sigmoid_centered_image_inverse_prime, L=L, x_0=x_0, k=k)
  )
}

get_sigmoid_linear_error_metric <- function(k=0.05) {
  error_metric <- get_linear_error_metric()
  L=2
  x_0 = 0
  compose_sigmoid_with_error(error_metric,
                             sigmoid=purrr::partial(sigmoid_centered_image, L=L, x_0=x_0, k=k),
                             sigmoid_prime=purrr::partial(sigmoid_centered_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse = purrr::partial(sigmoid_centered_image_inverse, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_prime = purrr::partial(sigmoid_centered_image_inverse_prime, L=L, x_0=x_0, k=k)
  )

}

get_sigmoid_relative_error_metric <- function(k=0.05) {
  error_metric <- get_relative_error_metric()
  L=2
  x_0 = 0
  compose_sigmoid_with_error(error_metric,
                             sigmoid=purrr::partial(sigmoid_centered_image, L=L, x_0=x_0, k=k),
                             sigmoid_prime=purrr::partial(sigmoid_centered_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse = purrr::partial(sigmoid_centered_image_inverse, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_prime = purrr::partial(sigmoid_centered_image_inverse_prime, L=L, x_0=x_0, k=k)
  )

}

chain_rule <- function(outer_prime, inner, inner_prime) {
  chained_fun <- function(...) {
    args_list = list(...)
    inner_vals <- rlang::inject(inner(!!!args_list))
    inner_prime_vals <- rlang::inject(inner_prime(!!!args_list))
    return(outer_prime(inner_vals) *inner_prime_vals)
  }
  return(chained_fun)
}

chain_rule_mqd <- function(outer_prime, inner, inner_prime) {
  return(function(m, q, d) {
    return(outer_prime(inner(m,q,d)) * inner_prime(m,q,d))
  })
}

chain_rule_mq <- function(outer_prime, inner, inner_prime) {
  return(function(m, q) {
    return(outer_prime(inner(m,q)) * inner_prime(m,q))
  })
}

#' Title
#'
#' @param error_metric
#' @param sigmoid Should be an increasing function
#' @param sigmoid_prime
#' @param sigmoid_inverse
#' @param sigmoid_inverse_prime
#'
#' @returns
#' @export
#'
#' @examples
compose_sigmoid_with_error <- function(error_metric, sigmoid, sigmoid_prime, sigmoid_inverse, sigmoid_inverse_prime) {

  f_inverse_q <- if (is.null(error_metric$f_inverse_q))  NULL else \(m,e) {
    error_metric$f_inverse_q(m, sigmoid_inverse(e))
    }

  f_prime_inverse_q=\(m,e) {
    error_metric$f_prime_inverse_q(m, sigmoid_inverse(e)) * sigmoid_inverse_prime(e)
  }

  f_single_inverse_q=\(m,e, d) {
    error_metric$f_single_inverse_q(m, sigmoid_inverse(e), d)
  }

  new_error_metric(glue::glue("{error_metric$name} with sigmoid"),
                   f_single=purrr::compose(sigmoid, error_metric$f_single),
                   f=purrr::compose(sigmoid, error_metric$f),
                   f_single_prime_m=chain_rule(sigmoid_prime, error_metric$f_single, error_metric$f_single_prime_m),
                   f_prime_m=chain_rule(sigmoid_prime, error_metric$f, error_metric$f_prime_m),
                   f_prime_inverse_q=f_prime_inverse_q,
                   f_single_increasing_q=error_metric$f_single_increasing_q,
                   f_increasing_q=error_metric$f_increasing_q,
                   f_single_inverse_q=f_single_inverse_q,
                   f_inverse_q=f_inverse_q,
                   check_domain=error_metric$check_domain,
                   d=error_metric$d
                   )

}




##### Error uitls

#' Useful when you want to calculate an error for a fix m matrix but for
#' multiple different q values.
#'
#' @param e_fun Takes a axb matrix and a long vector and returns a axb matrix
#' @param m Exd long matrix
#' @param q N long vector
#'
#' @returns a NxExd matrix
#' @export
#'
#' @examples
error_to_3d_array <- function(e_fun, m, q) {
  stopifnot(is.function(e_fun))
  stopifnot(is.matrix(m))
  stopifnot(is.vector(q))
  E = nrow(m)
  d = ncol(m)
  n = length(q)
  # first extend q and m to N*E and (N*E)xd respectively
  q_ext <- rep_len(q, n * E)
  row_indices <- rep(1:E, each = n)
  m_ext <- m[row_indices, ]

  e_values = e_fun(m_ext, q_ext)

  dim(e_values) = c(n, E, d)
  e_values

}
