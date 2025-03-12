
# Returns a list with the following elements:
# cdf: the CDF function
# pdf: the PDF function
# sample: a function to sample from the distribution
# cdf_inv: the inverse CDF function
# support: the support of the pdf (and range of the sample function and cdf_inv)
#' Distribution object
#'
#' @param cdf,pdf,cdf_inv of the distribution
#' @param sample function to sample from the distribution. Optional.
#' @param support of the distribution. An interval [a,b] such that the pdf is 0
#' outside of [a,b]. Optional
#'
#' @returns A new distribution object with class distribution
#' @export
#'
new_distribution <- function(cdf, pdf, support, sample=NULL, cdf_inv=NULL) {
  stopifnot(is.function(cdf))
  stopifnot(is.function(pdf))
  stopifnot(is.numeric(support) && length(support) == 2)
  stopifnot(is.null(sample) || is.function(sample))
  stopifnot(is.null(cdf_inv) || is.function(cdf_inv))
  return(structure(list(cdf=cdf, pdf=pdf, sample=sample, cdf_inv=cdf_inv, support=support), class = "distribution"))
}
