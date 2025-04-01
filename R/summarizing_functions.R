#' A class to hold a summarizing function
#'
#' @param name Descriptive name of summarizing function
#' @param f Function that takes a nx3 matrix of assessments and returns a nxd matrix
#' @param output_dim The d dimension of the summarizing function
#'
#'
new_summarizing_function <- function(f, name="", output_dim=NULL) {
  # check that arguments are functions
  stopifnot(is.function(f))
  stopifnot(is.character(name))
  stopifnot((is.numeric(output_dim) && length(output_dim) == 1) || is.null(output_dim))
  return(structure(list(name=name, f = f, output_dim=output_dim), class = "summarizing_function"))
}

get_median_summarizing_function <- function() {
  return(new_summarizing_function(m_median, "median", 1))
}

get_three_quantiles_summarizing_function <- function() {
  return(new_summarizing_function(m_three_quantiles, "three quantiles", 3))
}

#' Extracts the median from a set of assessments
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
