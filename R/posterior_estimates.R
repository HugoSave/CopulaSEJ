# file for producing estimates from log unnormalized posteriors

posterior_to_mean <- function(log_unnormalized_posterior, support, num_samples=1000) {
  samples <- sample_log_unnormalized_density(log_unnormalized_posterior, support, num_samples)
  mean(samples)
}

# compute simultaneous because they share the same samples
posterior_mean_and_likelihood <- function(log_unnormalized_posterior, support, true_value, num_samples=1000) {
  samples <- sample_log_unnormalized_density(log_unnormalized_posterior, support, num_samples)
  mean_value <- mean(samples)
  likelihood <- estimate_margin_kde(samples, support)$pdf(true_value)
  return(list(mean=mean_value, likelihood=likelihood))
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

    bay_settings <- list(iterations=num_samples, startValue=matrix(0.2, nrow=2, ncol=1)) # (support[1]+support[2])/2)
    bayesian_tools_samples <- BayesianTools::runMCMC(bayesian_setup, sampler="DEzs", settings=bay_settings)
    bayesian_tools_samples$Z
  }
  else{
    armspp::arms(num_samples, log_density, support[1] , support[2])
  }
}
