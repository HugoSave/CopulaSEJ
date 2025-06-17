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

#' Title
#'
#' @param name
#' @param f,f_prime_q (q,m) returns (NxExD)
#' @param f_increasing
#' @param D_tilde
#' @param fix_m
#' @param f_inverse
#' @param short_name
#' @param ideal_mean_var
#'
#' @returns
#' @export
#'
#' @examples
new_decoupler <- function(name, f, f_prime_q, f_increasing, D_tilde, fix_m=NULL, f_inverse=NULL, short_name=NULL, ideal_mean_var=NULL) {
  checkmate::assert_string(name)
  checkmate::assert_function(f, args=c("q", "m"), ordered=TRUE)
  checkmate::assert_function(f_prime_q, args=c("q", "m"), ordered=TRUE)
  checkmate::assert_function(f_increasing, args=c("m"), ordered=TRUE)
  checkmate::assert_function(fix_m, args=c("m"), ordered=TRUE, null.ok=TRUE)
  checkmate::assert_function(f_inverse, args=c("z", "m"), null.ok=TRUE)
  checkmate::assert_string(short_name, null.ok=TRUE)
  checkmate::assert_function(ideal_mean_var, null.ok=TRUE) # function of m and quantiles (with 0 100 support)
  checkmate::assert(
    checkmate::check_function(D_tilde),
    checkmate::check_count(D_tilde, positive=TRUE),
    checkmate::check_choice(D_tilde, c("D")),
  )
  short_name <- if (is.null(short_name)) name else short_name
  if (is.null(fix_m)) {
    fix_m = \(m) default_fix_m(m, f, f_prime_q, f_increasing, name, f_inverse)
  }
  if (is.character(D_tilde) && D_tilde == "D") {
    D_tilde_f = \(D) D
  } else if (is.numeric(D_tilde)) {
    D_tilde_f = \(D) D_tilde
  } else {
    D_tilde_f = D_tilde
  }

  return(structure(
    list(
      name = name,
      short_name=short_name,
      f = f,
      f_prime_q = f_prime_q,
      f_increasing=f_increasing,
      fix_m = fix_m,
      f_inverse=f_inverse,
      D_tilde=D_tilde_f,
      ideal_mean_var=ideal_mean_var
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
  m <- typecheck_and_convert_matrix_vector(m, z)
  D = ncol(m)
  E = nrow(m)
  N = length(z)
  z_rep_matrix = array(z, dim=(c(N,E,D)))
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  (z_rep_matrix + m_rep_matrix)
}

decoupler_linear_fix_m <- function(m, name="linear") {
  increasing_matrix = decoupler_linear_increasing(m)
  new_fixed_decoupler(
    name=name,
    f=\(q) decoupler_linear(q,m),
    f_prime_q=\(q) decoupler_linear_prime_q(q,m),
    f_increasing=\() increasing_matrix,
    D=ncol(m),
    f_inverse=\(z) decoupler_linear_inverse(z,m)
  )
}

preprocess_m_to_median <- function(m) {
  stopifnot(ncol(m)==3)
  m[,2,drop=FALSE]
}

preprocess_m_from_string <- function(m_string) {
  if (m_string == "median") {
    return(preprocess_m_to_median)
  } else if (m_string == "mean_G") {
    return(\(m) m_mean_estimate(m, quantiles_probs=c(0.05,0.50,0.95), overshoot=0.1, support_restriction = NULL, global_support=TRUE))
  } else if (m_string == "mean_E") {
    return(\(m) m_mean_estimate(m, quantiles_probs=c(0.05,0.50,0.95), overshoot=0.1, support_restriction = NULL, global_support=FALSE))
  } else {
    stop("Unknown m preprocessing function")
  }
}

m_mean_estimate <- function(assessments, quantiles_probs = c(0.05,0.50,0.95), overshoot=0.1, support_restriction=NULL, global_support=FALSE) {
  stopifnot(overshoot > 0) # some overshoot is needed for the interpolation to make sense
  if (is.data.frame(assessments)) {
    assessments = as.matrix(assessments)
  }

  # check that assessments has 3 cols
  if (ncol(assessments) != 3) {
    stop("Assessments must have 3 columns")
  }

  if (global_support) {
    quantiles <- add_0_and_100_percentiles_matrix(assessments, overshoot=overshoot, support_restriction=support_restriction)
  } else {
    quantiles <- add_0_and_100_percentiles_matrix_per_row(assessments, overshoot=overshoot, support_restriction=support_restriction)
  }

  means <- estimate_mean_from_quantiles(quantiles, cdf_values=quantiles_probs)
  matrix(means, nrow=nrow(assessments), ncol=1, dimnames = list(names(means), "mean"))
}

estimate_mean_from_quantiles <- function(quantiles, cdf_values) {
  stopifnot(!(0 %in% cdf_values))
  stopifnot(!(1 %in% cdf_values))
  D = ncol(quantiles)
  stopifnot(length(cdf_values) == (D-2))
  cdf_values <- c(0,cdf_values, 1)
  cdf_values_diff <- cdf_values[2:D] - cdf_values[1:(D-1)]
  cdf_values_matrix <- matrix(cdf_values_diff, nrow=nrow(quantiles), ncol=D-1, byrow=TRUE)
  # check that quantiles has 3 cols
  Matrix::rowSums(cdf_values_matrix * (quantiles[, 2:D] + quantiles[, 1:(D-1)])) / 2
}

clamp_vars_for_beta <- function(means_vars, epsilon=1e-6) {
  vars_limit <- means_vars$means * (1-means_vars$means) # upper limit for variance of a beta function
  means_vars$vars <- pmin(means_vars$vars, vars_limit)
  means_vars
}

get_linear_decoupler <- function(D_tilde, compose_sigmoid=TRUE, m_preprocess=NULL, name=NULL, short_name=NULL, k=0.05, overshoot=0.1) { #6e-3) {
  checkmate::assert(
    checkmate::check_function(m_preprocess, null.ok = TRUE),
    checkmate::check_string(m_preprocess, null.ok = TRUE)
  )
  name <- ifelse(is.null(name), "linear decoupler", name)
  short_name <- ifelse(is.null(short_name), "q-m", short_name)

  decoupler <- new_decoupler(
      name=name,
      f=decoupler_linear,
      f_prime_q=decoupler_linear_prime_q,
      f_increasing=decoupler_linear_increasing,
      fix_m = decoupler_linear_fix_m,
      f_inverse=decoupler_linear_inverse,
      short_name=short_name,
      D_tilde=D_tilde
    )

  if (!is.null(m_preprocess)) {
    if (is.character(m_preprocess)) {
      m_preprocess <- preprocess_m_from_string(m_preprocess)
    }
    decoupler <- compose_decoupler_m_preprocessing(m_preprocess, decoupler,
                                      name=name,
                                      short_name = short_name)
  }
  if (compose_sigmoid) {
    ideal_mean_var <- \(quantiles) {
      quantiles_extended <- add_0_and_100_percentiles_matrix(quantiles, overshoot=overshoot)
      sigmoid_composed_affine_decoupler_ideal_mean_var(decoupler, k=k, m=quantiles, quantiles=quantiles_extended, cum_prob = c(0,0.05, 0.5, 0.95, 1)) |>
        clamp_vars_for_beta()
    }
    L=1
    x_0 = 0
    decoupler_composed <- compose_sigmoid_with_decoupler(decoupler,
                                                         sigmoid_f=purrr::partial(sigmoid, k=k),
                                                         sigmoid_prime_f=purrr::partial(sigmoid_prime, k=k),
                                                         sigmoid_inverse_f = purrr::partial(sigmoid_inverse, k=k)
    )
    decoupler_composed$ideal_mean_var = ideal_mean_var
    return(decoupler_composed)
  }

  decoupler
}

decoupler_relative <- function(q, m, epsilon=0.001) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  q_rep_matrix = array(q, dim=(c(N, E, D)))
  (q_rep_matrix - m_rep_matrix) / (abs(m_rep_matrix) + epsilon)
}

decoupler_relative_prime_q <- function(q, m, epsilon=1e-6) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  1 / (abs(m_rep_matrix) + epsilon)
}

decoupler_relative_increasing <- function(m, epsilon=1e-6) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(is_positive_no_zero(abs(m)+epsilon))
}

decoupler_relative_inverse <- function(z, m, epsilon=1e-6) {
  m <- typecheck_and_convert_matrix_vector(m, z)
  D = ncol(m)
  E = nrow(m)
  N = length(z)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  z_rep_matrix = array(z, dim=(c(N,E,D)))
  z_rep_matrix * (abs(m_rep_matrix) + epsilon) + m_rep_matrix
}

get_relative_decoupler <- function(D_tilde, compose_sigmoid=FALSE, m_preprocess=NULL, name="relative decoupler", short_name="Rel.", k=6e-2, overshoot=0.1, epsilon=0.001) {
  checkmate::assert(
    checkmate::check_function(m_preprocess, null.ok = TRUE),
    checkmate::check_string(m_preprocess, null.ok = TRUE)
  )
  decoupler <- new_decoupler(
      name=name,
      f=\(q,m) decoupler_relative(q,m,epsilon=epsilon),
      f_prime_q=\(q,m) decoupler_relative_prime_q(q,m, epsilon=epsilon),
      f_increasing=decoupler_relative_increasing,
      f_inverse=\(z,m) decoupler_relative_inverse(z,m, epsilon=epsilon),
      short_name=short_name,
      D_tilde=D_tilde
    )


  if (!is.null(m_preprocess)) {
    if (is.character(m_preprocess)) {
      m_preprocess <- preprocess_m_from_string(m_preprocess)
    }
    decoupler <- compose_decoupler_m_preprocessing(m_preprocess, decoupler,
                                                   name=name,
                                                   short_name = short_name)
  }
  if (compose_sigmoid) {
    ideal_mean_var <- \(quantiles) {
      quantiles_extended <- add_0_and_100_percentiles_matrix(quantiles, overshoot=overshoot)
      sigmoid_composed_affine_decoupler_ideal_mean_var(decoupler, k=k, m=quantiles, quantiles=quantiles_extended, cum_prob = c(0,0.05, 0.5, 0.95, 1)) |>
        clamp_vars_for_beta()
    }
    L=1
    x_0 = 0
    decoupler_composed <- compose_sigmoid_with_decoupler(decoupler,
                                                sigmoid_f=purrr::partial(sigmoid_clamped, k=k),
                                                sigmoid_prime_f=purrr::partial(sigmoid_prime, k=k),
                                                sigmoid_inverse_f = purrr::partial(sigmoid_inverse, k=k)
    )
    decoupler_composed$ideal_mean_var = ideal_mean_var
    class(decoupler_composed) <- c(class(decoupler_composed), "relative_decoupler")
    decoupler_composed$epsilon=epsilon
    decoupler_composed$k=k
    return(decoupler_composed)
  }
  decoupler$epsilon=epsilon
  decoupler$k=k
  class(decoupler) <- c(class(decoupler), "relative_decoupler")
  decoupler
}

#### Decoupler Ratio
decoupler_ratio <- function(q, m) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  q_rep_matrix = array(q, dim=(c(N, E, D)))
  q_rep_matrix / m_rep_matrix
}

decoupler_ratio_prime_q <- function(q, m) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  D = ncol(m)
  E = nrow(m)
  N = length(q)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  1 / m_rep_matrix
}

decoupler_ratio_increasing <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(is_positive_no_zero(m))
}

decoupler_ratio_inverse <- function(z, m) {
  m <- typecheck_and_convert_matrix_vector(m, z)
  D = ncol(m)
  E = nrow(m)
  N = length(z)
  m_rep_matrix = array(m, dim=(c(E, D, N))) |> aperm(c(3,1,2)) # this flips it to NxExD
  z_rep_matrix = array(z, dim=(c(N,E,D)))
  z_rep_matrix * m_rep_matrix
}

decoupler_ratio_fix_m <- function(m, name="ratio") {
  increasing_matrix = decoupler_ratio_increasing(m)
  new_fixed_decoupler(
    name=name,
    f=\(q) decoupler_ratio(q,m),
    f_prime_q=\(q) decoupler_ratio_prime_q(q,m),
    f_increasing=\() increasing_matrix,
    D=ncol(m),
    f_inverse=\(z) decoupler_ratio_inverse(z,m)
  )
}

get_ratio_decoupler <- function(D_tilde, compose_sigmoid=FALSE,m_preprocess=NULL,name="ratio decoupler", short_name="Ratio", k=6e-2) {
  checkmate::assert(
    checkmate::check_function(m_preprocess, null.ok = TRUE),
    checkmate::check_string(m_preprocess)
  )
  decoupler <- new_decoupler(
      name=name,
      f=decoupler_ratio,
      f_prime_q=decoupler_ratio_prime_q,
      f_increasing=decoupler_ratio_increasing,
      fix_m = decoupler_ratio_fix_m,
      f_inverse=decoupler_ratio_inverse,
      short_name=short_name,
      D_tilde=D_tilde
    )

  if (!is.null(m_preprocess)) {
    if (is.character(m_preprocess)) {
      m_preprocess <- preprocess_m_from_string(m_preprocess)
    }
    decoupler <- compose_decoupler_m_preprocessing(m_preprocess, decoupler,
                                                   name=name,
                                                   short_name = short_name)
  }

  if (compose_sigmoid) {
    L=1
    x_0 = 1
    decoupler <- compose_sigmoid_with_decoupler(decoupler,
                                                sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                                                sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                                                sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k)
    )
  }
  decoupler
}

decoupler_ratio_mean <- function(q, m, quantile_probs=c(0.05,0.5,0.95), overshoot=0.1, support_restriction=NULL) {
  m <- typecheck_and_convert_matrix_vector(m, q)
  stopifnot(ncol(m) == 3)
  m_extended <- add_0_and_100_percentiles_matrix(m, overshoot=overshoot, support_restriction=support_restriction)
  means <- estimate_mean_from_quantiles(m_extended, quantile_probs)
  decoupler_ratio(q, means)
}

get_mean_linear_decoupler <- function(global_support=TRUE, quantile_probs=c(0.05,0.5,0.95), overshoot=0.1, support_restriction=NULL) {
  # here we take 3 quantiles, calculate the mean and then apply the ratio
  get_means <- function(m) {
    m_mean_estimate(m, quantiles_probs = quantile_probs, support_restriction=support_restriction, global_support=global_support, overshoot=overshoot)
  }
  f <- function(q, m) {
    means <- get_means(m)
    decoupler_linear(q, means)
  }
  f_prime_q <- function(q, m) {
    means <- get_means(m)
    decoupler_linear_prime_q(q, means)
  }
  f_increasing <- function(m) {
    means <- get_means(m)
    decoupler_linear_increasing(means)
  }
  f_inverse <- function(z, m) {
    means <- get_means(m)
    decoupler_linear_inverse(z, means)
  }
  f_fix_m <- function(m) {
    means <- get_means(m)
    decoupler_linear_fix_m(means, "mean linear")
  }
  name <- if(global_support) {
    "mean_G linear decoupler"
  } else {
    "mean_E linear decoupler"
  }
  short_name <- if(global_support) {
    "q-Mn_G"
  } else {
    "q-Mn_E"
  }
  return(
    new_decoupler(
      name=name,
      f=f,
      f_prime_q=f_prime_q,
      f_increasing=f_increasing,
      fix_m = f_fix_m,
      f_inverse=f_inverse,
      short_name=short_name,
      D_tilde=1
    )
  )
}

default_fix_m <- function(m, f, f_prime_q, f_increasing, name, f_inverse=NULL) {
  D = ncol(f_increasing(m))
  f_inverse_new <- function(z) {
    if (is.null(f_inverse)) {
      NULL
      }
    else  {
      f_inverse(z,m)
    }
  }
  new_fixed_decoupler(
    name = glue::glue("fixed {name}"),
    f = \(q) f(q, m),
    f_prime_q = \(q) f_prime_q(q, m),
    f_increasing = \() f_increasing(m),
    D=D,
    f_inverse = f_inverse_new
  )
}


compose_decoupler_m_preprocessing <- function(process_m, decoupler, name, short_name) {
  f <- function(q, m) {
    m_new <- process_m(m)
    decoupler$f(q, m_new)
  }
  f_prime_q <- function(q,m) {
    decoupler$f_prime_q(q, process_m(m))
  }
  f_increasing <- function(m) {
    decoupler$f_increasing(process_m(m))
  }
  fix_m <- function(m) {
    decoupler$fix_m(process_m(m))
  }
  f_inverse <- if (is.null(decoupler$f_inverse)) {
    NULL
  } else {
    function(z, m) {
      decoupler$f_inverse(z, process_m(m))
    }
  }

  return(
    new_decoupler(
      name=name,
      f=f,
      f_prime_q = f_prime_q,
      f_increasing=f_increasing,
      fix_m=fix_m,
      f_inverse=f_inverse,
      short_name=short_name,
      D_tilde = decoupler$D_tilde
    )
  )
}


get_scaled_linear_decoupler <- function(point_estimate="mean_G", quantile_probs=c(0.05,0.5,0.95), overshoot=0.1, support_restriction=NULL) {
  checkmate::assert_subset(
    point_estimate,
    c("mean_G", "mean_E", "median", "3Q"),
    empty.ok = FALSE
  )
  D_tilde = if (point_estimate == "3Q") {
    3
  } else {
    1
  }
  linear_decoupler <- get_linear_decoupler(D_tilde=D_tilde)
  ratio_decoupler <- get_ratio_decoupler(D_tilde=D_tilde)

  get_supports_and_means <- function(m) {
    m <- typecheck_and_convert_matrix_vector(m, vector())
    stopifnot(ncol(m) == 3)
    #m_extended <- add_0_and_100_percentiles_matrix(m, overshoot=overshoot, support_restriction=support_restriction)
    m_extended_per_expert <- add_0_and_100_percentiles_matrix_per_row(m, overshoot=overshoot, support_restriction=support_restriction)
    L_bound <- m_extended_per_expert[,1, drop=FALSE] # ensure that the lower bound is not below the mean
    U_bound <- m_extended_per_expert[,ncol(m_extended_per_expert), drop=FALSE] # ensure that the upper bound is not above the mean
    if (point_estimate == "median") {
      estimates <- preprocess_m_to_median(m)
    } else if (point_estimate == "mean_G") {
      estimates <- m_mean_estimate(m, quantiles_probs = quantile_probs, support_restriction=support_restriction, global_support=TRUE, overshoot=overshoot)
    } else if (point_estimate == "mean_E") {
      estimates <- m_mean_estimate(m, quantiles_probs = quantile_probs, support_restriction=support_restriction, global_support=FALSE, overshoot=overshoot)
    } else if (point_estimate == "3Q") {
      estimates <- m
    }
    else {
      stop(glue::glue("Unknown point estimate type: {point_estimate}"))
    }

    mat <- cbind(
      estimates,
      U_bound - L_bound # Nx1 matrix with lower bounds
    ) # Nx(D+1) matrix with estimates and widths
    rownames(mat) <- rownames(estimates)
    if (is.null(colnames(estimates))) {
      colnames(mat)[ncol(mat)] <-  "Width"
    } else {
      colnames(mat) <- c(colnames(estimates), "Width")
    }
    mat
  }

  f <- function(q, m) {
    D = ncol(m) - 1 # D
    E = nrow(m) # number of experts
    estimates = m[,1:D, drop=FALSE] # ExD
    mean_errors <- linear_decoupler$f(q, estimates) # returns NxExD
    N = dim(mean_errors)[1] # number of assessments
    support_widths <- m[,D+1, drop=TRUE] # E
    # extend widths to NxExD
    support_widths <- array(support_widths, dim=c(E, D, N)) |> aperm(c(3,1,2)) # NxExD
    mean_errors/support_widths
  }

  f_prime_q <- function(q, m) {
    D = ncol(m) - 1 # D
    support_widths <- m[,D+1, drop=TRUE] # E
    # extend widths to NxExD
    support_widths <- array(support_widths, dim=c(E, D, N)) |> aperm(c(3,1,2)) # NxExD
    1/ support_widths
  }

  f_increasing <- function(m) {
    D = ncol(m) - 1
    matrix(TRUE, nrow=nrow(m), ncol=D)
  }

  f_inverse <- function(z, m) {
    D = ncol(m) - 1 # D
    E = nrow(m) # number of experts
    N = length(z) # number of assessments
    # extend all components to matching dimensions
    z = array(z, dim=c(length(z), E, D)) # NxExD

    support_widths <- m[,D+1, drop=TRUE] # E
    # extend widths to NxExD
    support_widths <- array(support_widths, dim=c(E, D, N)) |> aperm(c(3,1,2)) # NxExD

    estimates <- m[,1:D, drop=FALSE] # ExD
    estimates <- array(estimates, dim=c(E, D, N)) |> aperm(c(3,1,2)) # NxExD

    z * support_widths + estimates
  }

  tmp_decoupler <- new_decoupler(
    name = "scaled linear tmp",
    f=f,
    f_prime_q = f_prime_q,
    f_increasing = f_increasing,
    f_inverse = f_inverse,
    D_tilde = 1
  )

  # create name as function of the point estimate type
  if (point_estimate == "mean_G") {
    point_estimate_name <- "mu_G"
  } else if (point_estimate == "mean_E") {
    point_estimate_name <- "mu_E"
  } else if (point_estimate == "median") {
    point_estimate_name <- "Md"
  } else if (point_estimate == "3Q") {
    point_estimate_name <- "3Q"
  } else {
    stop("Unknown point estimate type")
  }

  compose_decoupler_m_preprocessing(get_supports_and_means, tmp_decoupler, "scaled linear {point_estimate_name}", glue::glue("Sc.Lin.{point_estimate_name}"))
}


get_mean_ratio_decoupler <- function(global_support=TRUE, quantile_probs=c(0.05,0.5,0.95), overshoot=0.1, support_restriction=NULL) {
  # here we take 3 quantiles, calculate the mean and then apply the ratio
  get_means <- function(m) {
    m_mean_estimate(m, quantiles_probs = quantile_probs, support_restriction=support_restriction, global_support=global_support, overshoot=overshoot)
  }
  f <- function(q, m) {
    means <- get_means(m)
    decoupler_ratio(q, means)
  }
  f_prime_q <- function(q, m) {
    means <- get_means(m)
    decoupler_ratio_prime_q(q, means)
  }
  f_increasing <- function(m) {
    means <- get_means(m)
    decoupler_ratio_increasing(means)
  }
  f_inverse <- function(z, m) {
    means <- get_means(m)
    decoupler_ratio_inverse(z, means)
  }
  f_fix_m <- function(m) {
    means <- get_means(m)
    decoupler_ratio_fix_m(means, "mean ratio")
  }
  name <- if(global_support) {
    "mean_G ratio decoupler"
  } else {
    "mean_E ratio decoupler"
  }
  short_name <- if(global_support) {
    "q/Mn_G"
  } else {
    "q/Mn_E"
  }
  return(
    new_decoupler(
      name=name,
      f=f,
      f_prime_q=f_prime_q,
      f_increasing=f_increasing,
      fix_m = f_fix_m,
      f_inverse=f_inverse,
      short_name=short_name,
      D_tilde=1
    )
  )
}

##### CDF error
indep_CDF <- function(q, m, overshoot=0.1, k_percentiles=c(5,50,95), scale="linear", support_restriction=NULL){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles, scale, support_restriction)
  ind_cond_m$f(q)
}

indep_CDF_prime_q <- function(q, m, overshoot=0.1, k_percentiles=c(5,50,95), scale="linear", support_restriction=NULL){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles, scale, support_restriction)
  ind_cond_m$f_prime_q(q)
}

indep_CDF_increasing <- function(m){
  # fit distributions to m. m should be a matrix with E rows and D columns
  is_increasing_matrix <- matrix(TRUE, nrow=nrow(m), ncol=1)
}

decouple_CDF_inverse <- function(z, m, overshoot=0.1, k_percentiles=c(5,50,95), scale="linear", support_restriction=NULL){
  # fit distributions to m. m should be a matrix with E rows and D columns
  ind_cond_m <- get_cdf_indep_function_fix_m(m, overshoot, k_percentiles, scale, support_restriction)
  ind_cond_m$f_inverse(z)
}

uniform_mean_var <- function(E, D) {
  # corresponds to a uniform distribution
  uniform_mean <- 0.5
  uniform_var <- 1/12
  means <- matrix(uniform_mean, nrow=E, ncol=D)
  vars <- matrix(uniform_var, nrow=E, ncol=D)
  dim_names <- list(default_E_names(seq_len(E)), default_D_names(seq_len(D)))
  dimnames(means) <- dim_names
  dimnames(vars) <- dim_names

  list(
    means = means,
    vars = vars
  )
}


get_CDF_decoupler <- function(scale="linear", overshoot=0.1, k_percentiles=c(5,50,95), support_restriction=NULL) {
  decoupler <-
    new_decoupler(
      name="CDF",
      f=\(q,m) indep_CDF(q, m, overshoot=overshoot, k_percentiles=k_percentiles, scale=scale, support_restriction=support_restriction),
      f_prime_q=\(q,m) indep_CDF_prime_q(q, m, overshoot=overshoot, k_percentiles=k_percentiles, scale=scale, support_restriction=support_restriction),
      f_increasing=indep_CDF_increasing,
      fix_m = \(m) get_cdf_indep_function_fix_m(m, overshoot=overshoot, k_percentiles=k_percentiles, scale=scale, support_restriction=support_restriction),
      f_inverse=\(z, m) decouple_CDF_inverse(z, m, scale=scale),
      D_tilde = 1,
      ideal_mean_var = \(quantiles) uniform_mean_var(nrow(quantiles), 1), # always one dimensional output
      short_name="CDF"
    )
  class(decoupler) <- c(class(decoupler), "CDF_decoupler")
  decoupler$overshoot=overshoot
  decoupler$k_percentiles=k_percentiles
  decoupler$support_restriction=support_restriction
  decoupler

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
get_cdf_indep_function_fix_m <- function(m, overshoot=0.1, k_percentiles = c(5,50,95), scale="linear", support_restriction=NULL) {
  distributions <- distributions_from_percentile_matrix(m, overshoot=overshoot, k_percentiles=k_percentiles, scale=scale, support_restriction=support_restriction)
  support <- distributions$support
  distributions <- distributions$distributions
  is_increasing_matrix <- matrix(TRUE, nrow=nrow(m), ncol=1) # single CDF value per expert

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

  decoupler <- new_fixed_decoupler(
    name = "CDF",
    f=f,
    f_prime_q=f_prime_q_fix_m,
    f_increasing = f_increasing,
    D=1,
    f_inverse=f_inverse
  )
  decoupler$support <- support
  decoupler
}


###### Sigmoid functions to compose above error functions with

sigmoid <- function(x, k) {
  return(1 / (1 + exp(-k * x)))
}

sigmoid_clamped <- function(x, k, epsilon=1e-6) {
  # this is a sigmoid with output clamped to [epsilon, 1-epsilon]
  # this is useful for numerical stability
  stopifnot(epsilon > 0)
  stopifnot(epsilon < 0.5)
  clamped_value <- sigmoid(x, k)
  clamped_value <- pmax(clamped_value, epsilon)
  clamped_value <- pmin(clamped_value, 1 - epsilon)
  return(clamped_value)
}

sigmoid_prime <- function(x, k) {
  exp_val = exp(-k * (x))
  return(k * exp_val / ((1 + exp_val)^2))
}

sigmoid_inverse <- function(y, k) {
  return(log(y / (1-y)) / k)
}

shifted_sigmoid <- function(x, L, x_0, k) {
  # if x is Inf or -Inf return L/2 and -L/2 respectively
  stopifnot(L=1)
  stopifnot(x_0=0)
  return(L / (1 + exp(-k * (x - x_0))) - L / 2)
}

shifted_sigmoid_prime <- function(x, L, x_0, k) {
  exp_val = exp(-k * (x - x_0))
  return(k * L * exp_val / (1 + exp_val)^2)
}

shifted_sigmoid_inverse <- function(y, L, x_0, k) {
  return(x_0 - log(L / (y + L / 2) - 1) / k)
}

sigmoid_centered_image_inverse_prime <- function(y, L, x_0, k) {
  return(4*L/(k*(L - 2*y)*(L+2*y)))
}

sigmoid_primitive <- function(x, k) {
  return(log(1+exp(k * x))/k)
}

sigmoid_squared <- function(x, k) {
  return(sigmoid(x,k)^2)
}

sigmoid_squared_primitive <- function(x, k) {
  term1 <- 1/(k*(1+exp(k * x)))
  #term2 <- log1pexp(k*x)/k # this seemed to give more hidden numerical issues.
  term2 <- log(1+exp(k * x))/k
  return(term1 + term2 - 1) # -1 makes it 0 at x=-Inf
}

sigmoid_composed_primitive <- function(affine_decoupler, k, q, m, primitive) {
  # This corresponds to the inegral of $\phi(q, m)$ wrt q.

  derivs <- affine_decoupler$f_prime_q(q, m) # constant wrt q but plugging it in to get the output dimensions
  decoupler_vals <- affine_decoupler$f(q, m)
  primitive(decoupler_vals, k=k) / derivs
}

sigmoid_composed_affine_decoupler_ideal_mean_var <- function(affine_decoupler, k, m, quantiles, cum_prob=c(0, 0.05, 0.5, 0.95, 1)) {
  checkmate::assert_matrix(quantiles, ncols=length(cum_prob))
  stopifnot(0 %in% cum_prob)
  stopifnot(1 %in% cum_prob)
  # diff works column wise but we want it over the rows
  cum_prob_diff <- diff(cum_prob)

  # loop per expert
  means_and_vars <- purrr::array_branch(quantiles, margin=1) |> purrr::map2(seq_len(nrow(quantiles)), \(expert_quantiles, e) {
    primitive_values_squared <- sigmoid_composed_primitive(affine_decoupler, k, expert_quantiles, m, sigmoid_squared_primitive)[,e,] # DxEx\tilde{D} to D after automatic drop (\tilde{D}=1)
    primitive_values_diffs <- diff(primitive_values_squared)
    primitive_values <-sigmoid_composed_primitive(affine_decoupler, k, expert_quantiles, m, sigmoid_primitive)[,e,] # DxEx\tilde{D} to D after automatic drop (\tilde{D}=1)
    primitive_value_diffs <- diff(primitive_values)
    expert_quantile_diffs <- diff(expert_quantiles)
    expected_value <- sum(cum_prob_diff * primitive_value_diffs / expert_quantile_diffs)
    expected_value_squared <- sum(cum_prob_diff * primitive_values_diffs / expert_quantile_diffs)
    varaiance <- expected_value_squared - expected_value^2
    list(
      mean=expected_value,
      var=varaiance
    )
  }) |> purrr::list_transpose()
  mean_matrix <- matrix(means_and_vars$mean, nrow=nrow(quantiles), ncol=1, dimnames = list(rownames(quantiles), "D1"))
  var_matrix <- matrix(means_and_vars$var, nrow=nrow(quantiles), ncol=1, dimnames = list(rownames(quantiles), "D1"))
  list(
    means=mean_matrix,
    vars=var_matrix
  )
  # We are evaluating the primitive at the quantiles but also providing the quantiles for the decoupler function
  #matrix(expected_values, nrow=nrow(quantiles), ncol=1, dimnames = list(rownames(quantiles), "mean"))
}

get_sigmoid_ratio_decoupler <- function(k=6e-2) {
  decoupler <- get_ratio_decoupler()
  L=1
  x_0 = 1
  compose_sigmoid_with_decoupler(decoupler,
                             sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                             sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k)
  )
}

get_sigmoid_relative_decoupler <- function(k=1) {
  decoupler <- get_relative_decoupler(D_tilde=1, )
  L=1
  x_0 = 0
  decoupler_composed <- compose_sigmoid_with_decoupler(decoupler,
                             sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                             sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k))



  decoupler_composed
}

get_sigmoid_mean_ratio_decoupler <- function(k=6e-2, global_support=TRUE) {
  decoupler <- get_mean_ratio_decoupler(global_support=global_support)
  L=1
  x_0 = 1
  compose_sigmoid_with_decoupler(decoupler,
                                 sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                                 sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                                 sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k)
  )
}

get_sigmoid_support_ratio_decoupler <- function(k=1, global_support=TRUE, quantile_probs=c(0.05,0.5,0.95), overshoot=0.1, support_restriction=NULL) {
  decoupler <- get_support_ratio_decoupler(global_support=global_support, quantile_probs=quantile_probs, overshoot=overshoot, support_restriction=support_restriction)
  L=1
  x_0=0

  compose_sigmoid_with_decoupler(decoupler,
                             sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                             sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k))
}

get_sigmoid_mean_linear_decoupler <- function(k=6e-2, global_support=TRUE) {
  decoupler <- get_mean_linear_decoupler(global_support=global_support)
  L=1
  x_0 = 1
  compose_sigmoid_with_decoupler(decoupler,
                                 sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                                 sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                                 sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k)
  )
}

get_sigmoid_linear_error_metric <- function(k=0.05) {
  error_metric <- get_linear_error_metric()
  L=2
  x_0 = 0
  compose_sigmoid_with_error(error_metric,
                             sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                             sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_prime = purrr::partial(sigmoid_centered_image_inverse_prime, L=L, x_0=x_0, k=k)
  )

}

get_sigmoid_relative_error_metric <- function(k=0.05) {
  error_metric <- get_relative_error_metric()
  L=2
  x_0 = 0
  compose_sigmoid_with_error(error_metric,
                             sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                             sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_prime = purrr::partial(sigmoid_centered_image_inverse_prime, L=L, x_0=x_0, k=k)
  )
}

get_sigmoid_linear_decoupler <- function(k=0.05) {
  decoupler <- get_linear_decoupler(compose_sigmoid=FALSE)
  L=1
  x_0 = 0
  compose_sigmoid_with_decoupler(decoupler,
                             sigmoid_f=purrr::partial(shifted_sigmoid, L=L, x_0=x_0, k=k),
                             sigmoid_prime_f=purrr::partial(shifted_sigmoid_prime, L=L, x_0=x_0, k=k),
                             sigmoid_inverse_f = purrr::partial(shifted_sigmoid_inverse, L=L, x_0=x_0, k=k)
                             )
}

#' Title
#'
#' @param decoupler
#' @param sigmoid Should be an increasing function
#' @param sigmoid_prime
#' @param sigmoid_inverse
#'
#' @returns
#' @export
#'
#' @examples
compose_sigmoid_with_decoupler <- function(decoupler, sigmoid_f, sigmoid_prime_f, sigmoid_inverse_f) {
  checkmate::assert_class(decoupler, "decoupler")

  f <- \(q,m) {
    sigmoid_f(decoupler$f(q, m))
  }

  f_inverse <- if (is.null(decoupler$f_inverse))  NULL else \(z,m) {
    decoupler$f_inverse(sigmoid_inverse_f(z), m)
  }

  f_prime_q <- \(q,m) {
    sigmoid_prime_f(decoupler$f(q, m)) * decoupler$f_prime_q(q, m)
  }

  fix_m <- \(m) {
    fixed_decoupler <- decoupler$fix_m(m)
    new_fixed_decoupler(
      name=glue::glue("{decoupler$name} with sigmoid"),
      f = \(q) sigmoid_f(fixed_decoupler$f(q)),
      f_prime_q= \(q) sigmoid_prime_f(fixed_decoupler$f(q)) * fixed_decoupler$f_prime_q(q),
      f_increasing=fixed_decoupler$f_increasing,
      D=fixed_decoupler$D,
      f_inverse= if (is.null(fixed_decoupler$f_inverse)) NULL else \(z) fixed_decoupler$f_inverse(sigmoid_inverse_f(z))
    )
  }

  new_decoupler(glue::glue("{decoupler$name} with sigmoid"),
                   f=f,
                   f_prime_q = f_prime_q,
                   f_increasing = decoupler$f_increasing,
                   fix_m =fix_m,
                   f_inverse=f_inverse,
                short_name = glue::glue("s({decoupler$short_name})"),
                D_tilde = decoupler$D_tilde,
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
