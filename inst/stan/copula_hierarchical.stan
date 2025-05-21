# Function definition credit to https://spinkney.github.io/helpful_stan_functions/group__normal.html
functions {
  real multi_normal_cholesky_copula_lpdf(matrix U, matrix L) {
    int N = rows(U);
    int J = cols(U);
    matrix[J, J] Gammainv = chol2inv(L);
    return -N * sum(log(diagonal(L))) - 0.5 * sum(add_diag(Gammainv, -1.0) .* crossprod(U));
  }
}

data {
  int<lower=1> N;                // number of observations
  int<lower=2> D;                // number of dimensions
  matrix[N, D] normal_latent;                // copula-transformed data that is then transformed to latent normal space with e.g. qnorm
  real<lower=0> eta;            // fixed hyperparameter for LKJ prior
}

parameters {
  cholesky_factor_corr[D] L;    // Cholesky factor of correlation matrix
}

model {
  // LKJ prior with fixed eta (reflects strength of belief in independence)
  L ~ lkj_corr_cholesky(eta);

  // Copula likelihood
  target += multi_normal_cholesky_copula_lpdf(normal_latent | L);
}

 generated quantities {
   corr_matrix[D] Sigma;
   Sigma = multiply_lower_tri_self_transpose(L);
 }
