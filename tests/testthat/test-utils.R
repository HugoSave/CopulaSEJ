
test_that("split_dataframe_to_error_observations works", {
  m_list = list(
    c(70.84, 21.18, 15.00, 5.00, 65.00, 55.00),
     c(118.06, 35.30, 25.00, 10.00, 75.00, 65.00),
    c(177.10, 52.95, 27.00, 13.00, 77.00, 67.00)
  )
  m = matrix(unlist(m_list), ncol=3, byrow=FALSE)
  realizations = c(213.78, -0.74, 12.00, 4.00, 49.00, 55.69)
  # a deliberately 'unordered' data frame
  training_set = data.frame(
    `5th percentile` = m_list[[1]],
    `50th percentile` = m_list[[2]],
    `95th percentile` = m_list[[3]],
    realization = realizations,
    expert_id = c(2, 1, 2, 1, 2, 1),
    question_id = c(1, 4, 3, 1, 4, 3),
    study_id = c(1, 1, 1, 1, 1, 1),
    check.names = FALSE
  )
  error_metric <- get_ratio_error_metric()
  errors <- error_metric$f(m, realizations)
  question_1_rows <- c(1, 4)
  question_3_rows <- c(3, 6)
  question_4_rows <- c(2, 5)

  expert_1_rows <- c(2, 4, 6)
  expert_2_rows <- c(1, 3, 5)

  expected_output = list(
    c(errors[intersect(question_1_rows, expert_1_rows),], errors[intersect(question_1_rows, expert_2_rows),]),
    c(errors[intersect(question_3_rows, expert_1_rows),], errors[intersect(question_3_rows, expert_2_rows),]),
    c(errors[intersect(question_4_rows, expert_1_rows),], errors[intersect(question_4_rows, expert_2_rows),])
  ) |> do.call(what=rbind)
  dimnames(expected_output) <- list(c(1, 3, 4), NULL) # row dimnames correspond to question_id
  output<-split_dataframe_to_error_observations(training_set, error_metric, m_three_quantiles, c(5, 50, 95))
  expect_equal(output, expected_output)
})

test_that("support_intersection works", {

  support_list = list(
    list(support=c(0, 1)),
    list(support=c(0.5, 1.5)),
    list(support=c(0.25, 1.75))
  )
  expect_equal(support_intersection(support_list, nested_list = TRUE), c(0.5, 1))

  warning_support <- list(
    c(0.5, 1.5),
    c(-1, 0.4)
  )

  expect_warning(output <- support_intersection(warning_support), "Support is NULL")
  expect_null(output)

  support_single_value <- list(
    c(-1, 0.5),
    c(0.5, 1.5)
  )
  output <- expect_equal(support_intersection(support_single_value), c(0.5, 0.5))

})

test_that("linear_index_to_d_E works", {
  d_max = 2
  linear_index = 5
  d = 1
  e = 3
  expect_equal(linear_index_to_d_E(linear_index, d_max), list(d=d, e=e))

  expect_equal(linear_index_to_d_E(1, 3), list(d=1, e=1))
})

test_that("d_E_to_linear_index works", {
  d_max = 2
  d = 1
  E = 3
  linear_index = 5
  expect_equal(d_E_to_linear_index(d, E, d_max), linear_index)

  expect_equal(d_E_to_linear_index(1,1, 3), 1)
})

test_that("study_df_to_assessment_array works with 3 quantiles", {
  training_set <- get_scrambled_training_data()
  q_vals <- unique(training_set$question_id) |> sort()
  e_vals <- unique(training_set$expert_id) |> sort()
  expected_out <- array(NA, dim=c(length(q_vals), length(e_vals), 3))
  dimnames(expected_out) <- list(default_Q_names(q_vals), default_E_names(e_vals), c("5th percentile", "50th percentile", "95th percentile"))
  for (q in seq_along(q_vals)) {
    for (e in seq_along(e_vals)) {
      expected_out[q, e, ] <- training_set[training_set$question_id == q_vals[q] & training_set$expert_id == e_vals[e], 1:3] |> unlist()
    }
  }

  output <- study_df_to_assessment_array(training_set, m_three_quantiles, c(5, 50, 95))
  expect_equal(output, expected_out)
})

test_that("study_df_to_assessment_array works with median", {
  training_set <- get_scrambled_training_data()
  q_vals <- unique(training_set$question_id) |> sort()
  e_vals <- unique(training_set$expert_id) |> sort()
  expected_out <- array(NA, dim=c(length(q_vals), length(e_vals), 1))
  dimnames(expected_out) <- list(default_Q_names(q_vals), default_E_names(e_vals), c("50th percentile"))
  for (q in seq_along(q_vals)) {
    for (e in seq_along(e_vals)) {
      expected_out[q, e, ] <- training_set[training_set$question_id == q_vals[q] & training_set$expert_id == e_vals[e], 2] |> unlist()
    }
  }

  output <- study_df_to_assessment_array(training_set, m_median, c(5, 50, 95))
  expect_equal(output, expected_out)
})


test_that("study_df_single_question_to_assessment_matrix returns right dims with 3 quantiles", {
  test_set <- get_test_data()
  E <- length(unique(test_set$expert_id))
  d <- 3
  output <- study_df_single_question_to_assessment_matrix(test_set, m_three_quantiles, c(5, 50, 95))
  expect_equal(dim(output), c(E,d))
})

test_that("study_df_single_question_to_assessment_matrix returns right dims with median", {
  test_set <- get_test_data()
  E <- length(unique(test_set$expert_id))
  d <- 1
  output <- study_df_single_question_to_assessment_matrix(test_set, m_median, c(5, 50, 95))
  expect_equal(dim(output), c(E,d))
})


# Dummy error function for testing
mock_linear_error_fun <- function(estimates, realizations) {
  return(realizations-estimates)  # Simple absolute error computation
}

test_that("assessment_array_to_flattened_errors and assessment_array_to_errors works correctly", {

  # Explicit test case
  assessments <- abind::abind(matrix(1:9, nrow=3, byrow=TRUE), matrix(10:18, nrow=3, byrow=TRUE), along=3) |> aperm(c(3, 1, 2))

  realizations <- c(2, 8)  # Two instances

  # Expected calculations
  # First instance:
  # Estimates:     [1  2  3]
  #                [4  5  6]
  #                [7  8  9]
  # Realizations:  [2  2  2]
  #                [2  2  2]
  #                [2  2  2]
  # Errors:        [1  0 -1]
  #                [-2 -3 -4]
  #                [-5 -6 -7]
  expected_1 <- c(1, 0, -1, -2, -3, -4, -5, -6, -7)  # Flatten row-wise

  # Second instance:
  # Estimates:     [10 11 12]
  #                [13 14 15]
  #                [16 17 18]
  # Realizations:  [8  8  8]
  #                [8  8  8]
  #                [8  8  8]
  # Errors:        [-2 -3 -4]
  #                [-5 -6 -7]
  #                [-8 -9 -10]
  expected_2 <- c(-2, -3, -4, -5, -6, -7, -8, -9, -10)  # Flatten row-wise

  expected_output <- rbind(expected_1, expected_2)
  rownames(expected_output) <-  c("Q1", "Q2")
  colnames(expected_output) <- c("E1D1", "E1D2", "E1D3", "E2D1", "E2D2", "E2D3", "E3D1", "E3D2", "E3D3")

  result <- assessment_array_to_flattened_errors(assessments, realizations, mock_linear_error_fun)
  #result <- assessment_array_to_flattened_errors(assessments, realizations, get_linear_error_metric()$f)

  expect_equal(result, expected_output)

  result_2 <- assessment_array_to_errors(assessments, realizations, mock_linear_error_fun) |> flatten_3d_array_to_matrix()

  expect_equal(result_2, expected_output)
})




test_that("q_support_to_error_support works", {
  q1 <- c(2/10, 2/5)
  m1 <- 2
  q2 <- c(-2/5, -2/10)
  m2 <- -2
  q3 <- c(-1, 1)
  m3 <- 1

  error_metric <- get_ratio_error_metric()
  d_i = 1 # dim does not matter for ratio error metric
  expect_equal(q_support_to_error_support(q1, m1, d_i, error_metric), c(5,10))
  expect_equal(q_support_to_error_support(q2, m2, d_i, error_metric), c(5,10))
  expect_warning(out <- q_support_to_error_support(q3, m3, d_i, error_metric), "Error function is not monotonic or not continous over the provided support.")
  expect_null(out)
})

test_that("error_support_to_q_support works", {
  s1 <- c(5,10)
  m1 <- 2
  s2 <- c(5,10)
  m2 <- -2
  s3 <- c(-1,1)
  m3 <- 1

  q_support1 <- c(2/10, 2/5)
  q_support2 <- c(-2/5, -2/10)

  error_metric <- get_ratio_error_metric()
  d_i = 1 # dim does not matter for ratio error metric
  expect_equal(error_support_to_q_support(s1, m1, d_i, error_metric), q_support1)
  expect_equal(error_support_to_q_support(s2, m2, d_i, error_metric), q_support2)
  expect_warning(out <- error_support_to_q_support(s3, m3, d_i, error_metric), "Error function is not monotonic or not continous over the provided support.")
  expect_null(out)
})


test_that("error_supports_to_q_supports works", {

  elim11 <- 0.5
  elim12 <- 10
  elim21 <- 0.4
  elim22 <- 8
  error_supports <- list(
    c(elim11, elim12),
    c(elim21, elim22)
  )
  error_metric <- get_ratio_error_metric() # m/q. So inverse is m/e
  m1 = 1
  m2 = 3

  # Given m1=1 what q gives rise to erros within 0.5-10?
  support1 = c(0.1, 2)
  support2 = c(3/8, 3/0.4)
  expected_out <- list(support1, support2)

  output <- error_supports_to_q_supports(error_supports, matrix(c(m1, m2), nrow=1), error_metric)
  expect_equal(output, expected_out)
})

test_that("matrix_to_nested_list_by_row works", {
  m <- matrix(1:6, nrow=2, byrow=TRUE)
  expected_output <- list(
    list(1, 2, 3),
    list(4, 5, 6)
  )
  output <- matrix_to_nested_list_by_row(m)
  expect_equal(output, expected_output)
})
