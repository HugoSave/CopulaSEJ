#' A class to hold a summarizing function
#'
#' @param name Descriptive name of summarizing function
#' @param f Function that takes a nx3 matrix of assessments and returns a nxd matrix
#' @param output_dim The d dimension of the summarizing function
#' @param short_name Short name of the summarizing function. Should have no spaces and be filesystem friendly.
#'
#'
new_summarizing_function <- function(f, name="", output_dim=NULL, short_name="") {
  # check that arguments are functions
  stopifnot(is.function(f))
  stopifnot(is.character(name))
  stopifnot(is.character(short_name))
  stopifnot((is.numeric(output_dim) && length(output_dim) == 1) || is.null(output_dim))
  return(structure(list(name=name, f = f, output_dim=output_dim, short_name=short_name), class = "summarizing_function"))
}

#' The median from a set of assessments
#'
#' @returns
#' @export
#'
#' @examples
get_median_summarizing_function <- function() {
  return(new_summarizing_function(m_median, "median", 1, "Md"))
}

#' Title
#'
#' @returns
#' @export
#'
#' @examples
get_three_quantiles_summarizing_function <- function() {
  return(new_summarizing_function(m_three_quantiles, "three quantiles", 3, "3Q"))
}

#' Title
#'
#' @returns
#' @export
#'
#' @examples
get_mean_summarizing_function <- function() {
  return(new_summarizing_function(m_mean_estimate, "mean", 1, "Mn"))
}

#' Extracts
#'
#' @param assessments - Should be a nx3 matrix or data frame
#'
#' @returns - An n long vector of the median values of the assessments
#'
m_median <- function(assessments) {
  if (is.data.frame(assessments)) {
    assessments = as.matrix(assessments)
  }
  # check that assessments has 3 cols
  if (ncol(assessments) != 3) {
    stop("Assessments must have 3 columns")
  }
  assessments[,2,drop=FALSE]
}


#' Extracts the three quantiles from a set of assessments. Essentially is the
#' identity function.
#'
#' @param assessments - Should be a nx3 matrix or data frame
#'
#' @returns - The same as the input
#'
m_three_quantiles <- function(assessments) {
  if (is.data.frame(assessments)) {
    assessments = as.matrix(assessments)
  }
  if (ncol(assessments) != 3) {
    stop("Assessments must have 3 columns")
  }
  assessments
}

m_mean_estimate <- function(assessments, quantiles_probs = c(0.05,0.50,0.95), overshoot=0.1, support_restriction=NULL) {
  stopifnot(overshoot > 0) # some overshoot is needed for the interpolation to make sense
  if (is.data.frame(assessments)) {
    assessments = as.matrix(assessments)
  }

  # check that assessments has 3 cols
  if (ncol(assessments) != 3) {
    stop("Assessments must have 3 columns")
  }

  quantiles <- add_0_and_100_percentiles_matrix(assessments, overshoot=overshoot, support_restriction=support_restriction)

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
