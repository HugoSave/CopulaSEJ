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
get_mean_summarizing_function <- function(support_restriction=NULL, overshoot=0.1) {
  name <- "mean"
  short_name <- "Mn"
  return(new_summarizing_function(\(assessments) m_mean_estimate(assessments=assessments,
                                                                 overshoot=overshoot,
                                                                 support_restriction=support_restriction,
                                                                 unified_support=FALSE),
                                  name,
                                  1,
                                  short_name)
         )
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


