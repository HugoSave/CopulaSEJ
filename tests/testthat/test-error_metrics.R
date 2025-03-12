test_error_metric <- function(e_metric) {
  assessments = matrix(c(-7:-1, 1:8), ncol=3, byrow=TRUE)
  q = -2:2
  number_r = 5
  m_all = assessments
  d = ncol(assessments)
  specific_expert = 3 # arbitrary
  m_all_specific_expert = matrix(assessments[specific_expert,], nrow=1)
  number_c = 3
  number_c_median = 1
  m_median = assessments[,2]
  m_median_specific_expert = m_median[specific_expert]


  # Do output dimension checks
  for (f_name in c("f", "f_prime_m")) {
    f = e_metric[[f_name]]
    output_size = dim(f(m_all, q))
    expect_equal(output_size, c(number_r, number_c))
    output_size = dim(f(m_median, q))
    stopifnot(output_size == c(number_r, 1))
  }

  errors = e_metric$f(m_all, q)
  errors_median = e_metric$f(m_median, q)
  expect_equal(matrix(errors[,2], ncol=1), errors_median)

  for (f_name in c("f_inverse_q", "f_prime_inverse_q")) {
    f = e_metric[[f_name]]
    output_size = dim(f(m_all, errors))
    expect_equal(output_size, c(number_r, number_c))
    output_size = dim(f(m_median, errors_median))
    expect_equal(output_size, c(number_r, 1))
  }
  f_inc = e_metric$f_increasing_q
  stopifnot(all(dim(f_inc(m_all)) == dim(m_all)))
  stopifnot(dim(f_inc(m_median)) == c(number_r, 1))

  # check for inverse propoerty
  err_inv <- e_metric$f_inverse_q(m_all, errors) # should yield q regardless of column choice
  # every column should be equal to q.
  expect_equal(err_inv, matrix(q, nrow=length(q), ncol=ncol(m_all)))
}

test_that("dimensions of linear error metric are correct", {
  test_error_metric(get_linear_error_metric())
})

test_that("dimensions of ratio error metric are correct", {
  test_error_metric(get_ratio_error_metric())
})

test_that("ratio error is invaraint of m", {
  e_metric = get_ratio_error_metric()
  assessments = matrix(-7:7, ncol=3, byrow=TRUE)
  q = rep(10, nrow(assessments))
  values <- e_metric$f_prime_m(assessments, q)
  value = 1/10
  expect_equal(values, matrix(value, nrow=nrow(assessments), ncol=ncol(assessments)))
})

# test ratio error is increasing and decreasing correctly
test_that("ratio error is increasing and decreasing correctly", {
  e_metric = get_ratio_error_metric()
  assessments = matrix(-2:3, ncol=3, byrow=TRUE)
  values <- e_metric$f_increasing_q(assessments)
  expected_outputs = matrix(c(T, T, NA, F, F, F), ncol=3, byrow = TRUE)
  expect_equal(values, expected_outputs)
})


