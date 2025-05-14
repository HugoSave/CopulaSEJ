# file for producing estimates from log unnormalized posteriors

posterior_to_mean <- function(log_unnormalized_posterior, support, num_samples=1000) {
  samples <- sample_log_unnormalized_density(log_unnormalized_posterior, support, num_samples)
  mean(samples)
}

performance_metrics_list <- function(mean=NA, median=NA, neg_log_lik=NA, post_neg_log_lik=NA, sd=NA) {
  return(list(mean=mean, median=median, neg_log_lik=neg_log_lik, post_neg_log_lik=post_neg_log_lik, sd=sd))
}

# compute simultaneous because they share the same samples
posterior_performance_metrics <- function(log_unnormalized_posterior, support, true_value, num_samples=1000) {
  samples <- sample_log_unnormalized_density(log_unnormalized_posterior, support, num_samples)
  mean_value <- mean(samples)
  median_value <- median(samples)
  sd_value = sd(samples)

  d <- density(samples)
  likelihood <- approx(d$x, d$y, xout=true_value, method="linear", yleft =0, yright=0)$y
  post_neg_log_lik <- -log_unnormalized_posterior(true_value) # should be the same as looking at the samples but I don't think it is.
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

  neg_log_like <- -log(likelihood)
  return(performance_metrics_list( mean=mean_value, median=median_value, neg_log_lik=neg_log_like, post_neg_log_lik=post_neg_log_lik, sd=sd_value))
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
sample_log_unnormalized_density <- function(log_density, support, num_samples, start_point=NULL, method="BayesianTools") {
  #MCMCpack::MCMCmetrop1R(log_density, (support[1] + support[2])/2, mcmc=num_samples)$batch
  #as.vector(mcmc::metrop(log_density, (support[1] + support[2])/2, num_samples)$batch)
  if (is.infinite(support[1]) || is.infinite(support[2])) {
    stop("not yet implemented")
    checkmate::expect_number(start_point)
    mcmc::metrop(log_density, start_point, num_samples)$batch
  }  else if (method=="BayesianTools") {
    bayesian_setup <- BayesianTools::createBayesianSetup(log_density, lower=support[1], upper=support[2])

    bay_settings <- list(iterations=num_samples, startValue=matrix((support[1]+support[2])/2, nrow=2, ncol=1)) # (support[1]+support[2])/2)
    bayesian_tools_samples <- BayesianTools::runMCMC(bayesian_setup, sampler="DEzs", settings=bay_settings)
    bayesian_tools_samples$Z
  }
  else{
    armspp::arms(num_samples, log_density, support[1] , support[2])
  }
}
