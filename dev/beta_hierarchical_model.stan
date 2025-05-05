data {
  int<lower=1> N;
  real<lower=0,upper=1> Z[N];
  real<lower=0> alpha_a;
  real<lower=0> beta_a;
  real<lower=0> alpha_b;
  real<lower=0> beta_b;
}
parameters {
  real<lower=0> a;
  real<lower=0> b;
}
model {
  a ~ gamma(alpha_a, beta_a);
  b ~ gamma(alpha_b, beta_b);
  for (n in 1:N)
    Z[n] ~ beta(a, b);
}
