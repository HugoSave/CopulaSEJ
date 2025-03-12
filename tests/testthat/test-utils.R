
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
