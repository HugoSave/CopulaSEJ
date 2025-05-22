

mock_training_estimates <- array(runif(60), dim = c(5, 4, 3)) # Q=5, E=4, D=3
mock_realizations <- runif(5)


test_that("reject_experts returns valid structure for classical_calibration", {
  result <- reject_experts(mock_training_estimates, mock_realizations,
                           rejection_level = 0.2, test = "classical_calibration")

  expect_type(result, "list")
  expect_named(result, c("accepted_assessments", "p_values", "accepted_experts", "rejected_experts"))
  expect_true(is.numeric(result$p_values))
  expect_true(is.numeric(result$accepted_experts))
})

test_that("reject_experts uses distance_correlation test correctly", {
  result <- reject_experts(mock_training_estimates[,,2,drop=FALSE], mock_realizations,
                           rejection_level = 0.2, test = "distance_correlation", decoupler = get_linear_decoupler(D_tilde = 3))

  expect_type(result$p_values, "double")
  expect_true(all(result$p_values >= 0 & result$p_values <= 1))
})

test_that("adaptive_p_value_rejection keeps min_nr_experts", {
  p_values <- c(0.01, 0.02, 0.5, 0.9)
  min_nr_experts <- 3
  rejection_level <- 0.05

  result <- adaptive_p_value_rejection(p_values, rejection_level, min_nr_experts)
  expect_length(result, min_nr_experts)
  expect_true(all(result %in% 1:4))
})

test_that("p_values_test returns valid results for classical_calibration", {
  result <- p_values_test(mock_training_estimates, mock_realizations, test = "classical_calibration")
  expect_true(is.numeric(result))
  expect_length(result, dim(mock_training_estimates)[2])
})

test_that("p_values_test throws for wrong test type", {
  expect_error(p_values_test(mock_training_estimates, mock_realizations, test = "unknown_test"))
})

test_that("a", {
  # example here comes from study 33 and expert 10 with the CDF decoupler
  percentiles <- matrix(c(
    0.0e+00, 1.5e+00, 3.5e+00,
    5.5e+01, 6.8e+01, 7.5e+01,
    5.0e+01, 1.03e+02, 1.15e+02,
    3.0e+01, 7.0e+01, 9.0e+01,
    0.0e+00, 3.0e+00, 5.0e+00,
    5.0e+03, 1.00e+05, 2.50e+05,
    1.9e+02, 3.0e+02, 5.0e+02,
    1.0e+06, 1.00e+09, 1.00e+14,
    5.0e+00, 3.0e+01, 6.0e+01,
    5.0e+01, 7.0e+01, 9.0e+01,
    1.0e-01, 1.0e+00, 5.0e+00
  ), nrow = 11, byrow = TRUE,
  dimnames = list(
    paste0("Q", 1:11),
    c("5th", "50th", "95th")
  ))

  decoupler_values <- matrix(c(
    0.96320755,
    0.50000000,
    0.50000000,
    0.01000000,
    0.95771605,
    0.95342466,
    0.50000000,
    0.05000000,
    0.04233333,
    0.27500000,
    0.95256410
  ), ncol = 1, dimnames = list( paste0("Q", 1:11),NULL))

  # Bug in the dcorT.test, it does not implement the special case if the product of bias corrected variances becomes negative
  # expect_warning(energy::dcorT.test(percentiles, decoupler_values),
  #                regexp = "NaNs produced")
  # This in call sqrt(XX * YY)

  expect_no_warning(fixed_dcorT_test(percentiles, decoupler_values))
})
