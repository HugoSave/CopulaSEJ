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

  # check f_single_*. Only first dimension
  single_errors <- e_metric$f_single(m_median, q, 1)
  expect_vector(single_errors, size=length(q))
  output <- e_metric$f_single_prime_m(m_median, q, 1)
  expect_vector(output, size=length(q))
  output <- e_metric$f_single_increasing_q(m_median, 1)
  expect_vector(output, size=length(m_median))

  # check inversion
  inverted_vals <- e_metric$f_single_inverse_q(m_median, single_errors, 1)
  expect_equal(inverted_vals, q)

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

test_that("decoupler linear error metric is correct", {
  m_test = matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=FALSE)
  q = c(1, 2)
  expected_out = array(NA, dim=c(2,2,2))
  expected_out[1,,] = q[1] - m_test
  expected_out[2,,] = q[2] - m_test
  expect_equal(decoupler_linear(q, m_test), expected_out)
})

test_that("dimensions of linear error metric are correct", {
  test_error_metric(get_linear_error_metric())
})

test_that("dimensions of ratio error metric are correct", {
  test_error_metric(get_ratio_error_metric())
})

test_that("dimensions of q/m error metric are correct", {
  test_error_metric(get_q_over_m_error_metric())
})

test_that("dimensions of relative error metric are correct", {
  test_error_metric(get_relative_error_metric())
})

# test ratio error is increasing and decreasing correctly
test_that("q/m error metric is increasing and decreasing correctly", {
  e_metric = get_q_over_m_error_metric()
  assessments = matrix(-2:3, ncol=3, byrow=TRUE)
  values <- e_metric$f_increasing_q(assessments)
  expected_outputs = matrix(c(F, F, NA, T, T, T), ncol=3, byrow = TRUE)
  expect_equal(values, expected_outputs)
})


df_from_m_columns <- function(m_mat, y_mat, q_vec) {
  m_vals_flat = as.vector(m_mat)
  y_flat = as.vector(y_mat)
  df <- tibble::tibble(q=rep(q_vec, times=3), m=m_vals_flat, y=y_flat)
  df$m <- factor(df$m, levels=unique(df$m))
  df
}

test_that("get_sigmoid_q_over_m_error_metric works", {
  n = 100
  x_start = -10
  x_end = 10
  m = matrix(c(-1,1,3), nrow=n, ncol=3, byrow=TRUE)
  q = seq(x_start, x_end, length.out=n)
  e_metric <- get_sigmoid_q_over_m_error_metric(k=0.5)
  y <- e_metric$f(m, q)

  df_f <- df_from_m_columns(m, y, q)

  p <- ggplot2::ggplot(df_f, aes(x=q, y=y, color=m)) + geom_line()
  path <- tempfile(fileext=".png")
  ggsave(path, p)

  expect_snapshot_file(path, "sigmoid_q_over_m_error_metric.png")
})

test_that("get_sigmoid_q_over_m_error_metric prime works", {
  m = matrix(c(-1,1, -2, 2), nrow=2, ncol=2, byrow=TRUE)
  q = c(3,-1)
  e_metric <- get_sigmoid_q_over_m_error_metric(k=0.5)
  sigmoid_input <- matrix(c(3/-1, 3/1, -1/-2, -1/2), nrow=2, ncol=2, byrow=TRUE)
  sig_prime_val <- sigmoid_centered_prime(sigmoid_input, L=2, x_0=1, k=0.5)
  # prime is -q/m^2
  e_metric_prime_vals <- matrix(c(-3/(-1)^2, -3/1^2, -(-1)/(-2)^2, -(-1)/2^2), nrow=2, ncol=2, byrow=TRUE)
  expected_out <- sig_prime_val * e_metric_prime_vals

  expect_equal(e_metric$f_prime_m(m, q), expected_out)
})

test_that("get_sigmoid_q_over_m_error_metric inverse prime works", {
  m = matrix(c(-1,1, -2, 2), nrow=2, ncol=2, byrow=TRUE)
  epsilon = matrix(c(0.9,-0.9, -0.3, 0), nrow=2, ncol=2, byrow=TRUE)
  L=2
  x_0=1
  k=0.5
  sigmoid_inverse = x_0 - log((L - 2*epsilon) / (L + 2*epsilon)) / k
  # prime inverse of sigmoid function
  sigmoid_inverse_prime <- 4 * L /(k*(L- 2 * epsilon) * (L + 2 * epsilon))
  # inverse of e=q/m is q=m*e. Whose derivative wrt e is m
  # then from chain rule we get
  q_exptected = (m) * sigmoid_inverse_prime

  e_metric <- get_sigmoid_q_over_m_error_metric(k=k)

  expect_equal(e_metric$f_prime_inverse_q(m, epsilon), q_exptected)
})



test_that("ratio error derivative is invaraint of m", {
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


test_that("error_to_3d_array works", {
  e_metric = get_linear_error_metric()
  E = 4
  d=3
  n=2
  assessments = matrix(1:12, nrow=4, ncol=3, byrow=TRUE)
  q = c(1,2)
  # Exd output for q=1
  expected_output_q_1 <- matrix(c(0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11), nrow=4, ncol=3, byrow=TRUE)
  # Exd output for q=2
  expected_output_q_2 <- matrix(c(1,0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10), nrow=4, ncol=3, byrow=TRUE)

  output <- error_to_3d_array(e_metric$f, assessments, q)
  expect_equal(output[1,,], expected_output_q_1)
  expect_equal(output[2,,], expected_output_q_2)

  # check that the dimensions are correct
  expect_equal(dim(output), c(n,E,d))
})
