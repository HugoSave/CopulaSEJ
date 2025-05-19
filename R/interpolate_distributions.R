#' Constructs a list of distributions for each expert in the data frame. Does
#' so by linear interpolation of the cdf distribution for each expert.
#'
#' @param data - a data frame containing quantile data for each expert. Should
#' only contain a single question in a single study and the 0 and 100 percentile
#' must be part of the data.
#'
#' @param k_percentiles
#'
#' @returns a list of one distribution for each expert in ascending order by
#' expert_id.
#' @export
#'
interpolate_distributions <- function(data, k_percentiles=c(0,5,50,95,100), interpolation="linear") {
  percentile_col_names = k_percentiles_to_colname(k_percentiles)
  N= nrow(data)
  cdf_probs = k_percentiles/100
  data_ordered <- data |> dplyr::arrange(expert_id)
  cdf_values <- percentiles_from_dataframe(data_ordered, k_percentiles)

  if (interpolation == "linear") {
    distributions <- linear_distribution_interpolation_matrix(cdf_values, cdf_probs)
  } else if (interpolation == "spline") {
    distributions <- lapply(seq(N), function(i) {
      cubic_spline_distribution_estimation(cdf_values[i,], cdf_probs)
    })
  }else if (interpolation == "spline clamp") {
    distributions <- lapply(seq(N), function(i) {
      cubic_poly_distribution_estimation_cics(cdf_values[i,], cdf_probs)
    })
  }
  distributions
}

log_linear_distribution_interpolation <- function(x, cdf_values) {
  checkmate::assert_numeric(x, finite=TRUE, unique=TRUE, sorted=TRUE)
  checkmate::assert_true(min(x) > 0)
  checkmate::assert_numeric(cdf_values, lower=0, upper=1, sorted=TRUE)
  support = range(x)
  x_log = log(x)
  x_log_dist = linear_distribution_interpolation(x_log, cdf_values)
  adjusted_pdf <- function(x_input) {
    x_in_support = x_input > support[1] & x_input < support[2]
    pdf_in_support = x_log_dist$pdf(log(x_in_support)) / x_in_support
    x_input[!x_in_support] <- 0
    x_input[x_in_support] <- pdf_in_support
    x_input
  }
  adjusted_cdf <- function(x_input) {
    x_left_outside <- x_input < support[1]
    x_right_outside <- x_input > support[2]
    x_in_support <- !(x_left_outside | x_right_outside)
    cdf_in_support <- x_log_dist$cdf(log(x_in_support))
    x_input[x_left_outside] <- 0
    x_input[x_right_outside] <- 1
    x_input[x_in_support] <- cdf_in_support
    x_input
  }

  adjusted_cdf_inv <- function(p_input) {
    exp(x_log_dist$cdf_inv(p_input))
  }

  sample_fun <- function(n) {
    exp(x_log_dist$sample(n))
  }

  return(
    new_distribution(
      adjusted_cdf,
      adjusted_pdf,
      range(x),
      sample_fun,
      adjusted_cdf_inv
    )
  )
}


# Function to estimate a linear cdf distribution from a set of samples
# x: the fractile values
# cdf_values: the values of the CDF at x
# Returns a list with the following elements:
# cdf: the CDF function
# pdf: the PDF function
# sample: a function to sample from the distribution
# cdf_inv: the inverse CDF function
# support: the support of the pdf (and range of the sample function and cdf_inv)
linear_distribution_interpolation <- function(x, cdf_values) {
  checkmate::assert_numeric(x, finite=TRUE, unique=TRUE, sorted=TRUE)
  checkmate::assert_numeric(cdf_values, lower=0, upper=1, sorted=TRUE)
  stopifnot(0 %in% cdf_values)
  stopifnot(1 %in% cdf_values)
  # Linear interpolation for the CDF
  cdf_function <- approxfun(x, cdf_values, method = "linear", ties = "ordered", yleft=0, yright=1)

  # cdf inverse
  cdf_inv <- approxfun(cdf_values, x, method = "linear", ties = "ordered", yleft=NA, yright=NA)
  sample_fun <- function(n) cdf_inv(runif(n))

  # Calculate PDF values as differences of CDF divided by differences of x
  pdf_values <- diff(cdf_values) / diff(x)

  # Define the step function for the PDF
  pdf_function <- stepfun(x, c(0, pdf_values, 0)) # Start and end with 0 before the first and after the last interval
  return(new_distribution(cdf=cdf_function, pdf=pdf_function, support=range(x),
                          sample=sample_fun,
                          cdf_inv=cdf_inv,
                          mean=mean_of_linear_distribution(x, cdf_values),
                          median=cdf_inv(0.5)))
}

mean_of_linear_distribution <- function(x, cdf_values) {
  # mean of linear CDF
  N = length(x)
  x_lag_sum = x[1:(N-1)] + x[2:N]
  sum(diff(cdf_values) * x_lag_sum)/2
}

linear_distribution_interpolation_matrix <- function(x_matrix, cdf_values) {
  # Linear interpolation for the CDF
  apply(x_matrix, 1, linear_distribution_interpolation, cdf_values=cdf_values)
}
