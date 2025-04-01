get_example_expert_distributions <- function() {
  # return one normal(mean=1, sd=1) and one uniform(1, 3) distribution
  list(
    list(
      cdf = function(x) pnorm(x, mean=1, sd=1),
      pdf = function(x) dnorm(x, mean=1, sd=1),
      support = c(-Inf, Inf)
    ),
    list(
      cdf = function(x) punif(x, min=1, max=3),
      pdf = function(x) dunif(x, min=1, max=3),
      support = c(1, 3)
    )
  )
}


# test calc_is_increasing_flat
test_that("calc_is_increasing_flat flattens in right order", {
  error_metric = get_ratio_error_metric()
  # two experts and d=3
  m = matrix(c(-1, -2, 3, 0, 1, 2), nrow=2, ncol=3, byrow = TRUE)
  expected_output = c(T, T, F, NA, F, F)

  received_output <- calc_is_increasing_flat(error_metric, m)

  expect_equal(received_output, expected_output)
})

# test calc_cdf_values
test_that("calc_cdf_values_JC flattens in right order", {
  expert_distributions = get_example_expert_distributions()
  cdf1 = expert_distributions[[1]]$cdf
  cdf2 = expert_distributions[[2]]$cdf
  # should mean that first expert is inc, desc, desc and second expert is
  # inc, inc, desc
  is_increasing = c(T, F, F, T,T, F)
  d=3 # number of quantiles
  q = c(1, 2)
  expected_outputs = list(
    c(cdf1(1), 1-cdf1(1), 1-cdf1(1), cdf2(1), cdf2(1), 1-cdf2(1)),
    c(cdf1(2), 1-cdf1(2), 1-cdf1(2), cdf2(2), cdf2(2), 1-cdf2(2))
  ) |> do.call(what = rbind)

  received_output <- calc_cdf_values_JC(expert_distributions, q, is_increasing, d)

  expect_equal(received_output, expected_outputs)
})

test_that("calc_sum_log_pdf_values_JC is correct", {
  expert_distributions = get_example_expert_distributions()
  pdf1 = expert_distributions[[1]]$pdf
  pdf2 = expert_distributions[[2]]$pdf
  d=3
  q = c(1, 2, 3)
  expected_outputs = c(
    d*log(pdf1(1) * pdf2(1)),
    d*log(pdf1(2) * pdf2(2)),
    d*log(pdf1(3) * pdf2(3))
  )

  received_output <- calc_sum_log_pdf_values_JC(expert_distributions, q, d)

  expect_equal(received_output, expected_outputs)
})

test_that("posterior support is correct", {
  expert_distributions = get_example_expert_distributions()
  E=2
  d=3
  m_samples = matrix(c(-1, -2, -3, 1, 2, 3), nrow=E, ncol=d)
  error_metric = get_ratio_error_metric()
  c = copula::frankCopula(dim = E*d, param = 4.161)
  posterior <- create_log_unnormalized_posterior_JC(c, expert_distributions,error_metric, m_samples)

  expect_equal(posterior$support, c(1, 3))
})


test_that("log ratio error terms are correct", {
  # for the ratio error metric, the log error term is log(1/abs(q)) + log(abs(m_i^e * q^2))
  # for each q.
  error_metric <- get_ratio_error_metric()
  d = 3
  E = 2
  m_samples = matrix(c(-1, -2, -3, 1, 2, 3), nrow=E, ncol=d)
  q = c(1, -2, 3, -4)

  expected_outputs = vector("double", length = length(q))
  expected_outputs[] = 0
  for (q_index in 1:length(q)) {
    q_j = q[q_index]
    for (i in 1:d) {
      for (e in 1:E) {
        error_term = log(1/abs(q_j)) + log(abs(m_samples[e, i] * q_j^2))
        expected_outputs[q_index] = expected_outputs[q_index] + error_term
      }
    }
  }

  received_output <- calc_sum_log_error_terms_JC(error_metric, m_samples, q)

  expect_equal(received_output, expected_outputs)
})

test_that("calc_pdf_and_cdf_values_flat flattens and works correctly", {
  linear_e <- get_linear_error_metric() #q-m
  dists <- example_error_distributions() # 6 distributions
  D=3
  E=2
  m = matrix(-2:3, nrow=E, ncol=D, byrow = TRUE)
  q = seq(-5,5, length.out=10)
  n = length(q)
  e_vals <- error_to_3d_array(linear_e$f, m, q)

  cdf_values <- matrix(NA, nrow=n, ncol=6)
  pdf_values <- matrix(NA, nrow=n, ncol=6)

  for (e in seq_len(2)) {
    for (d in seq_len(3)){
      flatt_i <- d_E_to_linear_index(d,e, 3)
      pdf_values[,flatt_i] <- dists[[flatt_i]]$pdf(e_vals[,e,d])
      cdf_values[,flatt_i] <- dists[[flatt_i]]$cdf(e_vals[,e,d])
    }
  }

  output <- calc_pdf_and_cdf_values_flat(e_vals, dists)
  expect_equal(output$pdf_values, pdf_values)
  expect_equal(output$cdf_values, cdf_values)
})

test_that("calc_sum_log_e_prime_m works", {
  e <- get_ratio_error_metric() # the prime is equal to 1/q
  q= c(-1,1,2)
  m= matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
  expected_output <- c(
    log(1/abs(-1))*4,
    log(1/abs(1)) *4,
    log(1/abs(2)) * 4
  )
  output <- calc_sum_log_e_prime_m(e, m, q)
  expect_equal(output, expected_output)
})

test_that("calculate_marginal_posteriors works", {
  error_metric <- get_ratio_error_metric()
  shape1=1.3
  shape2=2.5
  error_marginals <- list(
    list(
      pdf = function(x) extraDistr::dnsbeta(x, shape1=shape1, shape2=shape2, min=0.3, max=3),
      cdf = function(x) extraDistr::pnsbeta(x, shape1=shape1, shape2=shape2, min=0.3, max=3),
      support = c(0.3,3)
    ),
    list(
      pdf = function(x) extraDistr::dnsbeta(x, shape1=shape1, shape2=shape2, min=-2, max=-0.8),
      cdf = function(x) extraDistr::pnsbeta(x, shape1=shape1, shape2=shape2, min=-2, max=-0.8),
      support = c(-2,-0.8)
    )
  )

  m_samples <- matrix(c(1,-2),nrow=2, ncol=1)
  is_increasing_flat <- c(FALSE, TRUE) # m/q is decreasing wrt q when m>0

  q_values <- c(0, 0.3, 0.8, 1.1, 2.5, 4)
  n=length(q_values)
  pdf_matrix <- matrix(NA, nrow=n, ncol=2)
  cdf_matrix <- matrix(NA, nrow=n, ncol=2)
  for (e_i in 1:2) {
    m_rep <- matrix(rep.int(m_samples[e_i,], n), nrow=n)
    error_vals <- error_metric$f(m_rep, q_values)
    pdf_raw_vals <- error_marginals[[e_i]]$pdf(error_vals)
    e_deriv <- error_metric$f_prime_m(m_rep, q_values)

    pdf_matrix[,e_i] <- pdf_raw_vals * abs(e_deriv)
    # is_increasing_flat is false
    if (is_increasing_flat[[e_i]]) {
      cdf_matrix[,e_i] <- error_marginals[[e_i]]$cdf(error_vals)
    } else{
      cdf_matrix[,e_i] <- 1 - error_marginals[[e_i]]$cdf(error_vals)
    }
  }
  pdf_matrix[1,] <- 0 # with ratio error q=0 is outside of the support.

  posterior_marginals  <- calculate_marginal_posteriors(error_marginals, m_samples, error_metric)

  expect_true(posterior_marginals[[1]]$support[2] > posterior_marginals[[1]]$support[1])
  expect_true(posterior_marginals[[2]]$support[2] > posterior_marginals[[2]]$support[1])

  pdf_values_out <- cbind(posterior_marginals[[1]]$pdf(q_values),
                          posterior_marginals[[2]]$pdf(q_values))
  cdf_values_out <- cbind(posterior_marginals[[1]]$cdf(q_values),
                          posterior_marginals[[2]]$cdf(q_values))


  expect_equal(pdf_values_out, pdf_matrix)
  expect_equal(cdf_values_out, cdf_matrix)
})

test_that("create_log_unnormalized_posterior support is correct", {
  E=2
  d=3
  m_samples = matrix(1:6, nrow=E, ncol=d)
  error_metric = get_ratio_error_metric()
  q_support <- c(min(m_samples), max(m_samples))
  e_supports <- target_q_support_to_error_supports(q_support, m_samples, error_metric)
  error_marginals = example_error_distributions(e_supports)
  # generate copula params from ruinf
  set.seed(1)
  params <- runif((E*d - 1) * (E*d)/2) / (E*d) # need to sum to less than 1 to guarantee positive definite
  c = copula::normalCopula(dim=E*d, param = params, dispstr = "un")
  posterior <- create_log_unnormalized_posterior(c, error_marginals, error_metric, m_samples)

  q_outside <- c(seq(0.01, min(m_samples), length.out=10), seq(6.1, 7, length.out=10))
  expect_equal(posterior$logDM(q_outside), rep(-Inf, length(q_outside)))
  q_inside <- seq(q_support[1] + 1e-3, q_support[2] - 1e-3, length.out=20)

  checkmate::expect_numeric(posterior$logDM(q_inside), finite=TRUE, any.missing=FALSE)
})

test_that("create_log_unnormalized_posterior is equivelent with calculate_marginal_posteriors", {
  E=2
  d=3
  m_samples = matrix(1:6, nrow=E, ncol=d)
  error_metric = get_ratio_error_metric()
  q_support <- c(min(m_samples), max(m_samples))
  e_supports <- target_q_support_to_error_supports(q_support, m_samples, error_metric)
  error_marginals = example_error_distributions(e_supports)
  # generate copula params from ruinf
  set.seed(1)
  params <- runif((E*d - 1) * (E*d)/2)
  c = copula::normalCopula(dim = E*d, param = params, dispstr = "un")
  posterior <- create_log_unnormalized_posterior(c, error_marginals,error_metric, m_samples)

  posteriror_marginals <- calculate_marginal_posteriors(error_marginals, m_samples, error_metric)
  posterior2 <- Vectorize(function(x) {
    cdf_values <- posteriror_marginals |>
      purrr::map_dbl(function(marginal) marginal$cdf(x))
    cop_value <- copula::dCopula(cdf_values, c)
    marginal_value <- posteriror_marginals |>
      purrr::map_dbl(function(marginal) marginal$pdf(x)) |>
      prod()
    return(cop_value * marginal_value)
  })

  q_vals <- seq(q_support[1], q_support[2], length.out=10)
  posterior1_vals <- posterior$logDM(q_vals)

  expect_equal(posterior$logDM(q_vals), log(posterior2(q_vals)))

  q_vals_left <- seq(q_support[1]-100, q_support[1], length.out=10)
  expect_equal(posterior$logDM(q_vals_left), rep(-Inf, length(q_vals_left)))
  expect_equal(posterior2(q_vals_left), rep(0, length(q_vals_left)))

  q_vals_right <- seq(q_support[2], q_support[2]+100, length.out=10)
  expect_equal(posterior$logDM(q_vals_right), rep(-Inf, length(q_vals_right)))
  expect_equal(posterior2(q_vals_right), rep(0, length(q_vals_right)))

})
