



flatten_matrix_row_by_row <- function(matrix) {
  as.vector(t(matrix))
}


#' Converts a matrix to a nested list by row
#'
#' @param matrix matrix
#'
#' @returns A nested list such that list[[n]][[m]] returns the element in the n-th
#' row and m-th column of the matrix
#'
matrix_to_nested_list_by_row <- function(matrix) {
  checkmate::assert_matrix(matrix)
  purrr::array_tree(matrix, c(1,2))
}


#' Splits data frame by questions and structures error order
#'
#' Takes a data frame with the columns expert_id, question_id, and the k_percentiles
#' and splits the data frame by question_id (in ascending order). For each sub
#' data frame the errors are calculated and put into a d*E long vector ordered
#' such that the errors from a single expert are kept adjacent. These vectors are then
#' put row by row in a matrix resulting in a nx(d*E) matrix where n is the
#' number of questions.
#'
#' @param dataframe frame data frame with columns expert_id, question_id, and the k_percentiles
#' @param error_metric error metric object
#' @param summarizing_function summarizing function
#' @param k_percentiles vector of percentiles
#'
#' @returns n x (d*E) matrix
#'
split_dataframe_to_error_observations <- function(dataframe,
                                                     error_metric,
                                                     summarizing_function,
                                                     k_percentiles = c(5, 50, 95)) {
  # Split the training set into error observations
  # training set contains n*E rows where n is the number of questions and E is the number of experts
  ordered_training <- dataframe |> dplyr::arrange(question_id, expert_id)
  assessments <- percentiles_from_dataframe(ordered_training, k_percentiles)
  m <- summarizing_function(assessments) # dim = (n*E, d)
  errors <- error_metric$f(m, ordered_training$realization) # dim = (n*E, d)
  # orders the second dimension such that a the errors for a single error type are kept adjacent
  # that is, after splitting by question id, the smaller error matrix is flattened in row order
  flattened_errors <- errors |> split.data.frame(ordered_training$question_id) |>
    purrr::map(flatten_matrix_row_by_row) |> do.call(what = rbind) # dim = (n, d*E)
  flattened_errors
}


#' Assessment 3d array to flattened errors
#'
#' @param assessments NxExd matrix
#' @param realizations N vector
#' @param error_fun function taking Mxd matrix and M vector returning Mxd matrix
#'
#' @returns Nx(d*E) matrix
#' @export
#'
assessment_array_to_flattened_errors <- function(assessments, realizations, error_fun) {
  checkmate::assert_array(assessments, d=3)
  checkmate::assert_numeric(realizations)
  checkmate::assert_function(error_fun, nargs=2)
  N = dim(assessments)[1]
  E = dim(assessments)[2]
  D = dim(assessments)[3]
  checkmate::assert_set_equal(N, length(realizations))

  over_Q_estimates <- purrr::array_branch(assessments, 1) # list with Exd matrices

  flattened_erros <- purrr::map2(over_Q_estimates, realizations, \(estimates, realization) {
    q_repeated_for_E <- rep(realization, E) # same realisation for each expert
    error_values <- error_fun(estimates, q_repeated_for_E) # an Exd matrix
    flatten_matrix_row_by_row(error_values)
  }) |> do.call(what = rbind) # Nx(d*E) matrix

  row_names <- paste0("Q", 1:N)
  colnames(flattened_erros) <- flattened_E_D_names(E, D)
  rownames(flattened_erros) <- row_names
  flattened_erros
}

default_Q_names <- function(Q_indices) {
  paste0("Q", Q_indices)
}

default_E_names <- function(E_indices) {
  paste0("E", E_indices)
}

flattened_E_D_names <- function(E, D) {
  purrr::map(1:E, \(i) paste0("E", i, "D", 1:D)) |> unlist()
}

default_D_names <- function(D_indices) {
  paste0("D", D_indices)
}


flatten_3d_array_to_matrix <- function(arr) {
  # go from QxExD to QX(E*D) where E*D is flattened row by row
  dims <- dim(arr)
  dimnames_1 <- dimnames(arr)[[1]]
  dimnames_2 <- dimnames(arr)[[2]]
  dimnames_3 <- dimnames(arr)[[3]]
  mat <- aperm(arr, c(1, 3, 2)) |> array(dim=c(dims[1], dims[3]*dims[2]))
  if (!is.null(dimnames_2) && !is.null(dimnames_3)) {
    new_dimnames <- dimnames_2 |> purrr::map(\(x) paste0(x, dimnames_3)) |> unlist()
    colnames(mat) <- new_dimnames
  }
  rownames(mat) <- dimnames_1
  mat
}

clamp_cdf_values <- function(cdf_values, epsilon=1e-6) {
  # Densities are defined on open margin intervals (0,1) so we clamp the values to avoid numerical issues
  pmin(pmax(cdf_values, epsilon), 1 - epsilon) # avoid numerical issues with copula evaluation and fitting
}

assessments_to_decoupler_observations <- function(assessments, realizations, indep_fun) {
  checkmate::assert_array(assessments, d=3)
  checkmate::assert_numeric(realizations)
  checkmate::assert_function(indep_fun)
  stopifnot(dim(assessments)[1] == length(realizations))
  N = dim(assessments)[1]
  E = dim(assessments)[2]
  D = dim(assessments)[3]

  over_Q_estimates <- purrr::array_branch(assessments, 1) # list with ExD matrices

  errors <- purrr::map2(over_Q_estimates, realizations, \(estimates, realization) {
    error_values <- indep_fun(realization, estimates) # an 1xExD matrix
  }) |> abind::abind(along=1) # NxExD

  #errors <- errors |> aperm(c(3,1,2)) # NxExD

  assessment_dimnames <- dimnames(assessments)
  dimnames(errors)[1:2] <- assessment_dimnames[1:2]
  if (is.null(assessment_dimnames) || is.null(assessment_dimnames[[1]])) {
    dimnames(errors)[[1]] <- paste0("Q", 1:N)
  }
  if (is.null(assessment_dimnames) || is.null(assessment_dimnames[[2]])) {
    dimnames(errors)[[2]] <- paste0("E", 1:E)
  }
  if (is.null(assessment_dimnames) || is.null(assessment_dimnames[[3]])) {
    dimnames(errors)[[3]] <- paste0("D", 1:D)
  }

  errors
}

assessment_array_to_errors <- function(assessments, realizations, error_fun) {
  checkmate::assert_array(assessments, d=3)
  checkmate::assert_numeric(realizations)
  checkmate::assert_function(error_fun, nargs=2)
  N = dim(assessments)[1]
  E = dim(assessments)[2]
  D = dim(assessments)[3]
  checkmate::assert_set_equal(N, length(realizations))

  over_Q_estimates <- purrr::array_branch(assessments, 1) # list with Exd matrices

  errors <- purrr::map2(over_Q_estimates, realizations, \(estimates, realization) {
    q_repeated_for_E <- rep(realization, E) # same realisation for each expert
    error_values <- error_fun(estimates, q_repeated_for_E) # an Exd matrix
  }) |> abind::abind(along=3) # ExDxN

  errors <- errors |> aperm(c(3,1,2)) # NxExD

  assessment_dimnames <- dimnames(assessments)
  dimnames(errors) <- assessment_dimnames
  if (is.null(assessment_dimnames) || is.null(assessment_dimnames[[1]])) {
    dimnames(errors)[[1]] <- paste0("Q", 1:N)
  }
  if (is.null(assessment_dimnames) || is.null(assessment_dimnames[[2]])) {
    dimnames(errors)[[2]] <- paste0("E", 1:E)
  }
  if (is.null(assessment_dimnames) || is.null(assessment_dimnames[[3]])) {
    dimnames(errors)[[3]] <- paste0("D", 1:D)
  }

  errors
}



#' Calculates the intersection of a list of supports.
#'
#' Throws warning if the support is 0 wide
#'
#' @param support_list,nested_list a named list of elements with a 'support'
#' element being a numeric vector of length 2 if nested_list=TRUE. Otherwise
#' directly a list with the numeric vectors as elements.
#'
#' @returns
#' @export
#'
#' @examples
support_intersection <- function(support_list, nested_list=FALSE) {
  if (nested_list) {
    support_list <- support_list |> purrr::map(\(x) x$support)
  }

  support <- support_list |>
    purrr::reduce(\(s1, s2) {
      c(max(s1[1], s2[1]), min(s1[2], s2[2]))
    })

  if (support[1] > support[2]) {
    warning("Support is NULL")
    return(NULL)
  }

  support
}

# returns the smallest [a,b] that contains all given segments
support_interval_closure <- function(support_list, nested_list=FALSE) {
  if (nested_list) {
    support_list <- support_list |> purrr::map(\(x) x$support)
  }

  support <- support_list |>
    purrr::reduce(\(s1, s2) {
      c(min(s1[1], s2[1]), max(s1[2], s2[2]))
    })

  support
}

d_E_to_linear_index <- function(d, E, d_max) {
  d + d_max * (E-1)
}

linear_index_to_d_E <- function(linear_index, d_max) {
  d <- (linear_index - 1) %% d_max + 1
  e <- (linear_index - 1) %/% d_max + 1
  list(d=d, e=e)
}

add_d_e_to_list <- function(list_, D) {
  list_ |> purrr::imap(\(x, i) {
    dE <- linear_index_to_d_E(i, D)
    x$d <- dE$d
    x$expert_id <- dE$e
    x
  })
}




error_support_to_q_support <- function(error_support, m, d_i, error_metric) {
  # Calculate the support of q given the support of the
  # error distribution
  checkmate::assert_number(m)
  new_support <- error_metric$f_single_inverse_q(c(m, m), error_support, d_i)
  is_inc <- error_metric$f_increasing_q(m)
  if (is_inc == FALSE) {
    new_support <- c(new_support[2], new_support[1])
  }
  if (new_support[1] > new_support[2]) {
    warning("Error function is not monotonic or not continous over the provided support. Returning NULL.")
    return(NULL)
  }
  new_support
}

clamp_array_inside_support <- function(array_data, support, epsilon=1e-3) {
  new_support <- support + c(epsilon, -epsilon)
  pmin(pmax(array_data, new_support[1]), new_support[2])
}

error_supports_to_q_supports <- function(error_intervals, m_matrix, error_metric) {
  # Calculate the support of q given the support of the
  # error distribution
  checkmate::assert_matrix(m_matrix, "numeric")
  checkmate::assert_list(error_intervals, len = length(m_matrix), any.missing = FALSE, types = c("numeric"))
  d_max = ncol(m_matrix)

  q_supports <- error_intervals |>
    purrr::imap(\(error_min_max, i) {
      d_E <- linear_index_to_d_E(i, d_max)
      m_observed <- m_matrix[d_E$e, d_E$d]
      error_support_to_q_support(error_min_max,  m_observed, d_E$d, error_metric)
    })

}

q_support_to_error_support <- function(q_support, m, d_i, error_metric) {
  # Calculate the support of the error distribution given the support of q
  checkmate::assert_number(m)
  new_support <- error_metric$f_single(c(m, m), q_support, d_i)
  is_inc <- error_metric$f_increasing_q(m)
  if (is_inc == FALSE) {
    new_support <- c(new_support[2], new_support[1])
  }
  if (new_support[1] > new_support[2]) {
    warning("Error function is not monotonic or not continous over the provided support. Returning NULL.")
    return(NULL)
  }
  new_support
}


correlation_df_from_flat_errors <- function(flat_errors) {
  corr_matrix <- flat_errors |> cor(method="kendall") |> round(2)
  a <- tibble::as_tibble(corr_matrix)
  a_variable_order <- colnames(a)
  a$variable1 <- factor(colnames(a), levels = colnames(a))
  a_long <- a |> tidyr::pivot_longer(!variable1, names_to = "variable2", values_to="correlation")
  a_long$variable2 <- factor(a_long$variable2, levels = a_variable_order)
  a_long$variable1 <- factor(a_long$variable1, levels = a_variable_order)
  a_long
}

distributions_from_percentile_matrix <- function(m, overshoot=0.1, k_percentiles = c(5,50,95), scale="linear", support_restriction=NULL) {
  m <- typecheck_and_convert_matrix_vector(m, vector())
  N = nrow(m)
  checkmate::assert_numeric(k_percentiles, len=ncol(m))
  L = min(m)
  U = max(m)
  if (scale=="linear"){
    new_support <- widen_support(c(L, U), overshoot, support_restriction = support_restriction)
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
  list(
  distributions=distributions,
  support=c(L_star, U_star)
  )
}
