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
test_that("calc_cdf_values flattens in right order", {
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

  received_output <- calc_cdf_values(expert_distributions, q, is_increasing, d)

  expect_equal(received_output, expected_outputs)
})

test_that("calc_sum_log_pdf_values is correct", {
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

  received_output <- calc_sum_log_pdf_values(expert_distributions, q, d)

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
