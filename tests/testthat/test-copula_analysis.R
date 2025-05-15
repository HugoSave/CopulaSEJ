
test_that("estimate_margins returns correct objects", {
  error_obs <- matrix(rnorm(40*4), nrow=40, ncol=4)
  margins <- estimate_margins(error_obs)
  expect_equal(length(margins), 4)
  for (i in 1:length(margins)) {
    # check that pdf and cdf and support are objects
    expect_type(margins[[i]]$pdf, "closure")
    expect_type(margins[[i]]$cdf, "closure")
    expect_type(margins[[i]]$support, "double")
    expect_equal(length(margins[[i]]$support), 2)
  }
})

test_that("estimate_margins fulfill basic probability properties", {
  error_obs <- matrix(c(0.3313687,-28.6216216,1.2500000,1.2500000,1.3265306,0.9876100,1.0080645,0.8939394,0.7341772), nrow=9, ncol=1)
  margin_dist <- estimate_margins(error_obs)[[1]]
  x = seq(-40, 10, length.out=1000)
  n = length(x)
  outs_pdf <- margin_dist$pdf(x)
  outs_cdf <- margin_dist$cdf(x)

  checkmate::test_numeric(outs_pdf, lower=0, finite=TRUE, any.missing = FALSE, len = n)
  checkmate::test_numeric(outs_cdf, lower=0, upper=1, finite=TRUE, any.missing = FALSE, len = n)
  expect_equal(outs_cdf[1], 0)
  expect_equal(outs_cdf[n], 1)
  expect_equal(outs_pdf[1], 0)
  expect_equal(outs_pdf[n], 0)

  #estimate_margins(error_obs, list(c(-5, 10)))
})

test_that("calculate_marginal_posteriors works out of boundary", {
  error_dists <- example_error_distributions()
  m_matrix <- example_test_matrix()
  error_metric <- get_ratio_error_metric()
  posteriors <- calculate_marginal_posteriors(error_dists, m_matrix, error_metric)
  error_supports <- purrr::map(error_dists, \(x) x$support)
  q_supports <- error_supports_to_q_supports(error_supports, m_matrix, error_metric)
  nr_supports <- length(q_supports)
  pdf_inside <- matrix(NA, nrow=nr_supports, ncol=98)
  cdf_inside <- matrix(NA, nrow=nr_supports, ncol=98)
  pdf_left <- matrix(NA, nrow=nr_supports, ncol=99)
  pdf_right <- matrix(NA, nrow=nr_supports, ncol=99)
  cdf_left <- matrix(NA, nrow=nr_supports, ncol=100) #cdf is well defined on the boundary in contrast to the pdf
  cdf_right <- matrix(NA, nrow=nr_supports, ncol=100)

  for (i in seq_along(q_supports)) {
    q_support <- q_supports[[i]]
    posterior <- posteriors[[i]]
    expect_equal(posterior$support, q_support)
    q_inside <- seq(q_support[1], q_support[2], length.out=100)
    pdf_inside[i,] <- posterior$pdf(q_inside)[2:99]
    cdf_inside[i,] <- posterior$cdf(q_inside)[2:99]

    q_left <- seq(q_support[1]-10, q_support[1], length.out=100)
    pdf_left[i,] <- posterior$pdf(q_left)[1:99]
    cdf_left[i,] <- posterior$cdf(q_left)
    q_right <- seq(q_support[2], q_support[2]+10, length.out=100)
    pdf_right[i,] <- posterior$pdf(q_right)[2:100]
    cdf_right[i,] <- posterior$cdf(q_right)
  }

  expect_true(all(pdf_inside > 0))
  expect_true(all(cdf_inside > 0))
  expect_true(all(cdf_inside < 1))

  expect_equal(pdf_left, matrix(0, nrow=nr_supports, ncol=99))
  expect_equal(pdf_right, matrix(0, nrow=nr_supports, ncol=99))

  expect_equal(cdf_left, matrix(0, nrow=nr_supports, ncol=100))
  expect_equal(cdf_right, matrix(1, nrow=nr_supports, ncol=100))
})

test_that("estimate_margin_beta works for single value", {
  obs_vec_single_values <- rep(1,20)
  support <- c(2,3) # values not inside support
  estimate_margin_beta(obs_vec_single_values, support=support, out_of_boundary = "clamp")
})

test_that("find_target_q_support works", {
  error_supports <- list(
    c(5,10),
    c(5,10),
    c(0.5, 2),
    c(-1,-0.2)
  )
  m_values <- matrix(c(8,5,2,-1), nrow=2, ncol=2, byrow=TRUE)
  error_metric <- get_ratio_error_metric()
  individual_supports <- list(
    c(8/10, 8/5),
    c(5/10, 5/5),
    c(2/2, 2/(0.5)),
    c(-1/(-1), -1/(-0.2))
  )
  target <- c(5/10, 5) # smallest support that still contain all supports
  output <- find_target_q_support(error_supports, m_values, error_metric)
  expect_equal(output, target)
})

test_that("target_q_support_to_error_supports works", {
  target <- c(0.5, 5)
  m_values <- matrix(c(8,5,2,-1), nrow=2, ncol=2, byrow=TRUE)
  error_metric <- get_ratio_error_metric()
  expected_error_supports <- list(
    c(8/5,8/0.5),
    c(5/5,5/0.5),
    c(2/5, 2/0.5),
    c(-1/0.5,-1/5)
  )

  output<-target_q_support_to_error_supports(target, m_values, error_metric)
  expect_equal(output, expected_error_supports)
})



