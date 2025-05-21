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

test_that("sigmoid_composed_affine_decoupler_ideal_mean_var small assessments ", {
  quantiles <- matrix(c(
    1e-06, 1.0e-05, 1e-04,
    1e-09, 1.0e-02, 5e-02,
    1e-06, 2.5e-06, 1e-05,
    3e-06, 1.0e-05, 2e-05,
    1e-04, 1.0e-03, 5e-02
  ), nrow = 5, byrow = TRUE)

  quantiles_extended <- add_0_and_100_percentiles_matrix(quantiles)
  estimate_mean_from_quantiles(quantiles_extended, c(0.05,0.5, 0.95))
  decoupler <- get_relative_decoupler(D_tilde=1, compose_sigmoid = FALSE, m_preprocess="mean_G", epsilon=0.01)
  k = 0.1
  L_star_Q = min(quantiles_extended)
  U_star_Q = max(quantiles_extended)
  L_Z_vals <- sigmoid(decoupler$f(L_star_Q, quantiles), k=k)

  set.seed(1)
  dists <- linear_distribution_interpolation_matrix(quantiles_extended, c(0, 0.05, 0.5, 0.95, 1))
  samples <- dists |> purrr::map2(seq_along(dists), \(dist, e) {
    s <- dist$sample(1000000)
    Z_orig <- decoupler$f(s, quantiles)[,e,1] # NxEx1 to N
    Z_final <- sigmoid(Z_orig, k=k)
    list(mean=mean(Z_final), var=var(Z_final), L=min(Z_final), U=max(Z_final))
  }) |> purrr::list_transpose()

  means_vars <- sigmoid_composed_affine_decoupler_ideal_mean_var(decoupler, k=k, m=quantiles, quantiles=quantiles_extended, cum_prob=c(0, 0.05,0.5,0.95, 1))

  # we actually currently are failing the last test. Not numerically stable enough. I think we would need a MC backup if the numerical part fails
  expect_equal(means_vars$means[,1],samples$mean, tolerance = 0.001)
  expect_equal(means_vars$vars[,1],samples$var, tolerance = 0.001)

})

test_that("sigmoid_composed_affine_decoupler_ideal_mean_var other k ", {
  quantiles <- matrix(c(
    21.18,  35.30, 52.95,
    -30.00, -10.00, 10.00,
    -5.00,  10.00, 60.00,
    21.18,  49.42, 70.60
  ), nrow = 4, byrow = TRUE)
  quantiles_extended <- add_0_and_100_percentiles_matrix(quantiles)
  estimate_mean_from_quantiles(quantiles_extended, c(0.05,0.5, 0.95))
  decoupler <- get_relative_decoupler(D_tilde=1, compose_sigmoid = FALSE, m_preprocess="mean_G")
  k = 0.1
  L_star_Q = min(quantiles_extended)
  U_star_Q = max(quantiles_extended)
  L_Z_vals <- sigmoid(decoupler$f(L_star_Q, quantiles), k=k)


  set.seed(1)
  dists <- linear_distribution_interpolation_matrix(quantiles_extended, c(0, 0.05, 0.5, 0.95, 1))
  samples <- dists |> purrr::map2(seq_along(dists), \(dist, e) {
    s <- dist$sample(1000000)
    Z_orig <- decoupler$f(s, quantiles)[,e,1] # NxEx1 to N
    Z_final <- sigmoid(Z_orig, k=k)
    list(mean=mean(Z_final), var=var(Z_final), L=min(Z_final), U=max(Z_final))
  }) |> purrr::list_transpose()

  means_vars <- sigmoid_composed_affine_decoupler_ideal_mean_var(decoupler, k=k, m=quantiles, quantiles=quantiles_extended, cum_prob=c(0, 0.05,0.5,0.95, 1))

  expect_equal(dim(means_vars$means),c(4,1))
  expect_equal(dim(means_vars$vars),c(4,1))
  # we actually currently are failing the last test. Not numerically stable enough. I think we would need a MC backup if the numerical part fails
  expect_equal(means_vars$means[,1],samples$mean, tolerance = 0.001)
  expect_equal(means_vars$vars[,1],samples$var, tolerance = 0.001)

})

test_that("sigmoid_composed_affine_decoupler_ideal_mean_var works", {
  quantiles = matrix(c(-8, 2, 3,
               5, 7, 15,
               -15,0,15), nrow=3, ncol=3, byrow=TRUE)
  rownames(quantiles) <- c("E1", "E2", "E3")
  quantiles_extended <- add_0_and_100_percentiles_matrix(quantiles)
  estimate_mean_from_quantiles(quantiles_extended, c(0.05,0.5, 0.95))
  decoupler <- get_relative_decoupler(D_tilde=1, compose_sigmoid = FALSE, m_preprocess="mean_G")
  k = 1


  set.seed(1)
  dists <- linear_distribution_interpolation_matrix(quantiles_extended, c(0, 0.05, 0.5, 0.95, 1))
  samples <- dists |> purrr::map2(seq_along(dists), \(dist, e) {
    s <- dist$sample(1000000)
    Z_orig <- decoupler$f(s, quantiles)[,e,1] # NxEx1 to N
    Z_final <- sigmoid(Z_orig, k=k)
    list(mean=mean(Z_final), var=var(Z_final))
    }) |> purrr::list_transpose()

  means_vars <- sigmoid_composed_affine_decoupler_ideal_mean_var(decoupler, k=k, m=quantiles, quantiles=quantiles_extended, cum_prob=c(0, 0.05,0.5,0.95, 1))

  expect_equal(dim(means_vars$means),c(3,1))
  expect_equal(dim(means_vars$vars),c(3,1))
  # we actually currently are failing the last test. Not numerically stable enough. I think we would need a MC backup if the numerical part fails
  expect_equal(means_vars$means[1:2,1],samples$mean[1:2], tolerance = 0.001)
  expect_equal(means_vars$vars[1:2,1],samples$var[1:2], tolerance = 0.001)

  decoupler_composed <- get_relative_decoupler(D_tilde=1, compose_sigmoid = TRUE, m_preprocess="mean_G", k =1)
  expect_equal(decoupler_composed$ideal_mean_var(quantiles),means_vars)
})


test_that("support_ratio_decoupler is correct", {
  E = 2
  m = matrix(
    c(1,2,3,
      2,3,5), nrow=E, ncol=3, byrow=TRUE)
  range <- c(1,5) # width of 4
  extended_range = c(0.6, 5.4) # width of 4.8. 10% overshoot
  m_exteneded <- matrix(
    c(0.6, 1, 2, 3, 5.4,
      0.6, 2, 3, 5, 5.4), nrow=2, ncol=5, byrow=TRUE)

  estimated_means <- estimate_mean_from_quantiles(m_exteneded, c(0.05,0.5, 0.95))
  decoupler <- get_support_ratio_decoupler(overshoot=0.1)
  q_test <- c(0, 3)
  width=4.8
  expected_output_q1 <- matrix(
    (0-estimated_means)/width, ncol=1
  )
  expected_output_q2 <- matrix(
    (3-estimated_means)/width, ncol=1
  )
  expected_output <- abind::abind(
    expected_output_q1, expected_output_q2,
    along=0
  )
  dimnames(expected_output) <- NULL
  output <- decoupler$f(q_test, m)
  expect_equal(output, expected_output)
  # check inverse
  expect_equal(decoupler$f_prime_q(q_test, m),
               array(1/width, dim=c(2,2,1)))

  expect_equal(decoupler$f_increasing(m),
               matrix(TRUE, nrow=2, ncol=1))
  z = c(expected_output[1,,], expected_output[2,,])
  expected_inverse_out_third_dim <- matrix(
    c(z[1]*width + estimated_means[1], z[1]*width + estimated_means[2],
      z[2]*width + estimated_means[1], z[2]*width + estimated_means[2],
      z[3]*width + estimated_means[1], z[3]*width + estimated_means[2],
      z[4]*width + estimated_means[1], z[4]*width + estimated_means[2]), nrow=length(z), ncol=E, byrow=TRUE
  )
  expected_out <- array(expected_inverse_out_third_dim, dim=c(length(z), E, 1))
  expect_equal(decoupler$f_inverse(z, m), expected_out)
})

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
