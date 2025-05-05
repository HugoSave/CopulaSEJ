test_that("estimate_mean_from_quantiles works", {
  quantiles <- matrix(
    c(-2, -1, 0, 1,
      -3, 0, 1, 2),
    nrow=2,
    ncol=4,
    byrow = TRUE
    )
  cdf_values <- c(0.3, 0.4)

  # to help visualize the distributions and calculate by hand the expected means
  # cdf_values_with_01 <- c(0, cdf_values, 1)
  # dists <- linear_distribution_interpolation_matrix(
  #   quantiles,
  #   cdf_values_with_01
  # )
  # plot_distributions(dists)
  expected_means <- c(-0.2, 0.5)

  means <- estimate_mean_from_quantiles(
    quantiles,
    cdf_values
  )

  expect_equal(
    means,
    expected_means
  )
})
