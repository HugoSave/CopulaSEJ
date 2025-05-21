library(rstan)
options(mc.cores = parallel::detectCores()-2)
rstan_options(auto_write = TRUE)

# priors. Gives a gamma with mean 1 and sd k/k^2=1/k

D = 3
rho_params <- c(0.3, -0.2, 0.5)
norm_cop <- copula::normalCopula(rho_params, dim = D, dispstr="un")
normal_cop_1D_samp <- copula::rCopula(100, norm_cop)
Z_latent <- qnorm(normal_cop_1D_samp) # This is a bit of a weird interface to the copula function but it is what it is

N = nrow(normal_cop_1D_samp)
D = ncol(normal_cop_1D_samp)
eta = 10

data_list <- list(
  N = N,
  D = D,
  eta = eta,
  normal_latent=Z_latent
)

# load copula_hiarch_model.stan
mod <- load_stan_copula_model(force_compilation=FALSE)

fit_map <- rstan::optimizing(
  mod,
  data = data_list,
  verbose = TRUE
)

fit_map$par

# sigma values are the second half of returned parameters

nr_params = D*D*2
sigma_values <- fit_map$par[(nr_params/2+1):nr_params]
sigma_matrix <- matrix(sigma_values, nrow = D, ncol = D)

fitted_norm_cop <- copula::normalCopula(copula::P2p(sigma_matrix), dim = D, dispstr="un")
