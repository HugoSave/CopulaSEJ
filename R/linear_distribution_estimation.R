
# Function to estimate a linear cdf distribution from a set of samples
# x: the fractile values
# cdf_values: the values of the CDF at x
# Returns a list with the following elements:
# cdf: the CDF function
# pdf: the PDF function
# sample: a function to sample from the distribution
# cdf_inv: the inverse CDF function
# support: the support of the pdf (and range of the sample function and cdf_inv)
linear_distribution_estimation <- function(x, cdf_values) {
  # Linear interpolation for the CDF
  cdf_function <- approxfun(x, cdf_values, method = "linear", ties = "ordered", yleft=0, yright=1)

  # cdf inverse
  cdf_inv <- approxfun(cdf_values, x, method = "linear", ties = "ordered", yleft=NA, yright=NA)
  sample_fun <- function(n) cdf_inv(runif(n))

  # Calculate PDF values as differences of CDF divided by differences of x
  pdf_values <- diff(cdf_values) / diff(x)

  # Define the step function for the PDF
  pdf_function <- stepfun(x, c(0, pdf_values, 0)) # Start and end with 0 before the first and after the last interval
  return(new_distribution(cdf_function, pdf_function, range(x), sample_fun, cdf_inv))
}

linear_distribution_estimation_matrix <- function(x_matrix, cdf_values) {
  # Linear interpolation for the CDF
  apply(x_matrix, 1, linear_distribution_estimation, cdf_values=cdf_values)
}
