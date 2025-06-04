# file for producing estimates from log unnormalized posteriors

performance_metrics_list <- function(mean=NA, median=NA, neg_log_lik=NA, post_neg_log_lik=NA, sd=NA, cum_prob=cum_prob) {
  return(list(mean=mean, median=median, neg_log_lik=neg_log_lik, post_neg_log_lik=post_neg_log_lik, sd=sd, cum_prob=cum_prob))
}

# compute simultaneous because they share the same samples
posterior_performance_metrics <- function(log_unnormalized_posterior, support, true_value, num_samples=1000, mean_value=NULL, median_value=NULL, sample_prior=NULL) {
  #if (is.null(mean_value) || is.null(median_value)) {

  samples <- sample_log_unnormalized_density(log_unnormalized_posterior, support, num_samples, sample_prior=sample_prior)
  mean_value <- if(is.null(mean_value)) mean(samples) else mean_value
  median_value <- if(is.null(median_value)) median(samples) else median_value
  sd_value = sd(samples)
  # get what percentile true_value is in
  cum_prob <- ecdf(samples)(true_value)



  d <- density(samples)
  likelihood <- approx(d$x, d$y, xout=true_value, method="linear", yleft =0, yright=0)$y
  neg_log_like <- -log(likelihood)


  #} else {
  #  sd_value = NA
  #  neg_log_like <- NA
  #}

  # is unnormalized so not really useful but potentially interesting
  post_neg_log_lik <- -log_unnormalized_posterior(true_value)

  #f <- approxfun(d$x, d$y, yleft=0, yright=0)

  # This is fancy but seems a bit less robust.
  # fit <- ks::kde.boundary(x=samples, xmin=support[1], xmax=support[2], boundary.kernel="beta")
  # f<- approxfun(fit$eval.points, fit$estimate, yleft=0, yright=0)

  # below is not robust enough for some samples.
  # try integration
  # tryCatch(
  #   cumulative_probability<-integrate(f, min(fit$eval.points), true_value, subdivisions = 100000, stop.on.error = FALSE)$value,
  #   error = function(e) {
  #     # If integration fails, return NA
  #     cumulative_probability <- NA
  #   }
  # )

  return(performance_metrics_list( mean=mean_value, median=median_value, neg_log_lik=neg_log_like, post_neg_log_lik=post_neg_log_lik, sd=sd_value, cum_prob=cum_prob))
}


#' Title
#'
#' @param log_density
#' @param support
#' @param num_samples
#' @param start_point Needs to be provided if the support is infinite
#'
#' @returns
#' @export
#'
#' @examples
sample_log_unnormalized_density <- function(log_density, support, num_samples, start_point=NULL, method="BayesianTools", sample_prior=NULL) {
  #MCMCpack::MCMCmetrop1R(log_density, (support[1] + support[2])/2, mcmc=num_samples)$batch
  #as.vector(mcmc::metrop(log_density, (support[1] + support[2])/2, num_samples)$batch)
  if (is.infinite(support[1]) || is.infinite(support[2])) {
    stop("not yet implemented")
    checkmate::expect_number(start_point)
    mcmc::metrop(log_density, start_point, num_samples)$batch
  }  else if (method=="BayesianTools") {
    bayesian_setup <- BayesianTools::createBayesianSetup(log_density, lower=support[1], upper=support[2])

    if (!is.null(sample_prior)) {
      # if we have a prior, we can use it to get starting values
      checkmate::assert_function(sample_prior)
      starting_values <- sample_prior(50)
      Z <- sample_prior(100)
    } else {
      # otherwise, we use the support to get starting values
      starting_values <- get_starting_values(support, N=50, max_width_for_uniform = 1000)
      Z <- get_starting_values(100)
    }

    bay_settings <- list(iterations=num_samples, Z=matrix(Z, length(Z), ncol=1), startValue=matrix(starting_values, nrow=length(starting_values), ncol=1))
    bayesian_tools_samples <- BayesianTools::runMCMC(bayesian_setup, sampler="DEzs", settings=bay_settings)
    bayesian_tools_samples$Z
  }
  else{
    armspp::arms(num_samples, log_density, support[1] , support[2])
  }
}

get_starting_values <- function(support, N=1000, max_width_for_uniform = 1000) {
  if (is.infinite(support[1]) || is.infinite(support[2])) {
    stop("not yet implemented")
  }
  width = support[2] - support[1]
  if (width < max_width_for_uniform) {
    # if the support is very small, just return the middle point
    return(runif(N, support[1], support[2]))
  }
  if (support[1] > 0 & support[2] >0) {
    return(exp(runif(N, log(support[1]), log(support[2]))))
  } else if (support[1] < 0 & support[2] < 0) {
    return(-exp(runif(N, log(-support[2]), log(-support[1]))))
  } else {
    return(runif(N, support[1], support[2]))
  }
}
