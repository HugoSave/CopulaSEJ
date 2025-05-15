# Credit to andrjohns on the mc-stan forum for this 4 parameter beta function
functions {
  real beta_4p_lpdf(real y, real alpha, real beta, real lower, real upper) {
    // Scale 4-parameter Beta RV to 2-parameter Beta RV
    real x = (y - lower) / (upper - lower);

    // Return scaled 2-parameter beta lpdf
    return beta_lpdf(x | alpha, beta) - log(upper - lower);
  }
}

data {
  int<lower=1> N;
  real Z[N];
  real L;
  real U;
  real<lower=0> alpha_a;
  real<lower=0> beta_a;
  real<lower=0> alpha_b;
  real<lower=0> beta_b;
}
parameters {
  real<lower=0> A;
  real<lower=0> B;
}
model {
  A ~ gamma(alpha_a, beta_a);
  B ~ gamma(alpha_b, beta_b);
  for (n in 1:N)
    Z[n] ~ beta_4p(A, B, L, U);
}
