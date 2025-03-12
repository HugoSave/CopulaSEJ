
#' Constructor for error metric class
#'
#' @param name name of the function
#' @param f,f_prime_m The error metric and the derivative of f with respect to
#' m. Both functions have two arguments (m,q) where m is a nxd matrix and
#' q is a vector of length n. Returns a matrix of the same shape as m.
#' @param f_inverse_q,f_prime_inverse_q The inverse of f with respect to
#' q and the derivative of the inverse respectively. Have two arguments (m,e)
#' where m and e are vectors with length d and L. Outputs a matrix of size Lxd
#' @param f_increasing_q Function of a single argument (m). f_inverse_q is
#' optional.
#' Outputs a Boolean vector or matrix in the same shape as m. True if f is
#' increasing in q, False if not, NA if not defined.
#'
#' @returns A new error metric object
#' @export
#'
new_error_metric <- function(name, f, f_prime_m, f_prime_inverse_q, f_increasing_q, f_inverse_q=NULL) {
  stopifnot(is.character(name))
  if (!is.function(f) || !is.function(f_prime_m) ||
      !is.function(f_prime_inverse_q) || !is.function(f_increasing_q) ||
      (!is.null(f_inverse_q) && !is.function(f_inverse_q))) {
    stop("f,f_prime_m,f_prime_inverse_q and f_increasing_q must all be functions. f_inverse_q must be null or a function optional")
  }
  return(structure(list(name=name, f = f, f_prime_m = f_prime_m, f_prime_inverse_q = f_prime_inverse_q, f_increasing_q=f_increasing_q, f_inverse_q=f_inverse_q), class = "error_metric"))
}

get_valid_rows_error_metric <- function(m_values, q, error_metric) {
  # exclude questions with no assessments
  e <- error_metric(m_values, q)
  # return row indicies with finite values
  e_finite = is.finite(e)
  return(which(apply(e_finite, 1, all)))
}

linear_error <- function(m,q) {
  m <- typecheck_and_convert_matrix_vector(m,q)
  return(q-m)
}

linear_error_prime_q <- function(m,q) {
  m <- typecheck_and_convert_matrix_vector(m,q)
  m[] = -1
  return(m)
}

linear_error_inverse_q_2 <- function(m,epsilon) {
  stopifnot(is.vector(epsilon))
  stopifnot(is.vector(m))
  epsilon_rep_matrix = outer(epsilon, rep.int(1L, length(m)))
  m_rep_matrix = outer(rep.int(1L, length(epsilon)), m)
  return(m_rep_matrix+epsilon_rep_matrix)
}

linear_error_inverse_q <- function(m,epsilon) {
  # convert m and epsilon to matrices in case they are provided as vectors
  argcheck <- typecheck_and_convert_matrix_matrix(m,epsilon)
  m <- argcheck[[1]]
  epsilon <- argcheck[[2]]

  return(m+epsilon)
}


linear_error_prime_inverse_q <- function(m,epsilon) {
  # TODO e is not a vector but could be a matrix.
  argcheck <- typecheck_and_convert_matrix_matrix(m,epsilon)
  m <- argcheck[[1]]
  epsilon <- argcheck[[2]]

  return(matrix(1, nrow(m), ncol(m)))

  stopifnot(is.vector(e), is.vector(m))
  return(matrix(1, length(e), length(m)))
}

linear_error_prime_inverse_q_2 <- function(m,epsilon) {
  # TODO e is not a vector but could be a matrix.
  m <- typecheck_and_convert_matrix_vector(m, vector())
  epsilon <- typecheck_and_convert_matrix_vector(epsilon, vector())
  return(matrix(1, ncol(m), length(m)))
}

linear_error_increasing_q <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(matrix(TRUE, nrow(m), ncol(m)))
}


# linear error metric
get_linear_error_metric <- function() {
  return(new_error_metric("linear error q-m", linear_error,
                          linear_error_prime_q, linear_error_prime_inverse_q,
                          linear_error_increasing_q, linear_error_inverse_q))
}

is_positive_no_zero <- function(x) {
  na_values = which(x == 0)
  true_values = which(x >0)
  false_values = which(x < 0)
  x[true_values] = TRUE
  x[false_values] = FALSE
  x[na_values] = NA
  x
}

# checks if first argument is a vector or matrix and if second argument is a
# vector. Converts first argument to a matrix.
typecheck_and_convert_matrix_vector <- function(m,q) {
  stopifnot(is.vector(q))
  if (is.data.frame(m)) {
    m <- as.matrix(m)
  }
  stopifnot(is.vector(m) || is.matrix(m))
  if (is.vector(m)) {
    # make it a single column matrix
    m <- matrix(m, ncol=1)
  }
  m
}

typecheck_and_convert_matrix_matrix <- function(m,e) {
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
    m <- matrix(m, ncol=1)
  }
  if (is.vector(e)) {
    e <- matrix(e, ncol=1)
  }
  stopifnot(all(dim(m) == dim(e)))
  list(m,e)
}

ratio_error <- function(m,q) {
  m <- typecheck_and_convert_matrix_vector(m,q)
  # repeats q over each column of m to end up with the same size.
  q_rep_matrix = outer(q, rep.int(1L, ncol(m)))
  return(m/q_rep_matrix)
}

ratio_error_prime_q <- function(m,q) {
  m <- typecheck_and_convert_matrix_vector(m,q)
  q_rep_matrix = outer(q, rep.int(1L, ncol(m)))
  return(1/q_rep_matrix)
}

# epsilon is a vector with length L
# m is d
# output is Lxd
ratio_error_inverse_q <- function(m, epsilon) {
  args_checked <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- args_checked[[1]]
  epsilon <- args_checked[[2]]
  return(m / epsilon)
  epsilon_rep_matrix = outer(epsilon, rep.int(1L, length(m)))
  m_rep_matrix = outer(rep.int(1L, length(epsilon)), m)
  return(m_rep_matrix/(epsilon_rep_matrix))
}

ratio_error_prime_inverse_q <- function(m, epsilon) {
  args_checked <- typecheck_and_convert_matrix_matrix(m, epsilon)
  m <- args_checked[[1]]
  epsilon <- args_checked[[2]]

  return(-m/(epsilon^2))

  # - m/e^2
  # repeat the epsilon over multiple columns
  epsilon_rep_matrix = outer(epsilon, rep.int(1L, length(m)))
  # repeat m for each row
  m_rep_matrix = outer(rep.int(1L, length(epsilon)), m)
  return(-m_rep_matrix/(epsilon_rep_matrix^2))
}

ratio_error_increasing_q <- function(m) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  return(!is_positive_no_zero(m))
}

get_ratio_error_metric <- function() {
  return(new_error_metric("ratio error m/q", ratio_error, ratio_error_prime_q, ratio_error_prime_inverse_q, ratio_error_increasing_q, ratio_error_inverse_q))
}



