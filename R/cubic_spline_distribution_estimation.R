
cubic_spline_distribution_estimation <- function(x, cdf_values) {
  support = range(x)
  # Linear interpolation for the CDF
  # monoH.FC seem to have a bug to not actually always produce monotonic hermite
  # functions
  #spline_fit <- splinefun(x, cdf_values, method = "monoH.FC")
  spline_fit <- splinefun(x, cdf_values, method = "hyman")

  cdf_function <- function(values) {
    intervals = findInterval(values, support)
    values[intervals==0] <- 0
    values[intervals==2] <- 1
    values[intervals==1] <- spline_fit(values[intervals==1])
    values
  }


  pdf_function <- function(values) {
    intervals = findInterval(values, support)
    outside_support <- which(values <= support[1] | values >= support[2], arr.ind=TRUE)
    inside_support <- which(values > support[1] & values < support[2], arr.ind=TRUE)
    # return spline fit if inside support otherwise 0
    values[intervals == 0 | intervals == 2] <- 0
    values[intervals==1] <- spline_fit(values[intervals==1], deriv=1)
    values
  }

  return(list(cdf = cdf_function, pdf = pdf_function, support=range(x)))
}

cubic_spline_distribution_estimation_matrix <- function(x_matrix, cdf_values) {
  # Linear interpolation for the CDF
  apply(x_matrix, 1, cubic_poly_distribution_estimation, cdf_values=cdf_values)
}

cubic_poly_distribution_estimation_cics <- function(x, cdf_values) {
  # Linear interpolation for the CDF
  pdf(file = NULL)
  cdf_params <- ICSsmoothing::cics_explicit(x, cdf_values, c(0,0))
  dev.off()
  cdf_polynoms <- cdf_params$spline_polynomials

  pdf_polynoms <- map(cdf_polynoms, function(cdf) as.function(deriv(cdf)))
  cdf_polynoms <- map(cdf_params$spline_polynomials, as.function)

  nr_polynomials = length(cdf_polynoms)


  cdf_func <- function(values) {
    polynomial_index = findInterval(values, x)
    values[polynomial_index == 0] <- 0
    values[polynomial_index == (nr_polynomials+1)] <- 1

    for (i in 1:nr_polynomials) {
      values[polynomial_index == i] <- cdf_polynoms[[i]](values[polynomial_index == i])
    }
    values
  }

  pdf_func <- function(values) {
    ret = values
    polynomial_index = findInterval(values, x)
    ret[polynomial_index == 0 | polynomial_index == (nr_polynomials+1)] <- 0

    for (i in 1:nr_polynomials) {
      ret[polynomial_index == i] <- pdf_polynoms[[i]](values[polynomial_index == i])
    }
    ret
  }

  return(list(cdf = cdf_func, pdf = pdf_func, support=range(x)))
}

test_cubic_poly_distribution_estimation_cics <- function() {
  # test the find median function
  expert_belief_1 = list(cumprob = c( 0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                         fractile = c(-5.5, -3.5, -1.5, -0.7, 0.5, 2.5, 4.5))

  cubic_poly_distribution_estimation_cics(expert_belief_1$fractile, expert_belief_1$cumprob)
  #plot
  x = seq(-6, 8, by = 0.01)
  plot(x, cubic_poly_distribution_estimation_cics(expert_belief_1$fractile, expert_belief_1$cumprob)$pdf(x), type = "l", col = "red", xlab = "x", ylab = "Density")
}
