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
  cdf_values <- as.matrix(data_ordered[,percentile_col_names])
  if (interpolation == "linear") {
    distributions <- linear_distribution_estimation_matrix(cdf_values, cdf_probs)
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
