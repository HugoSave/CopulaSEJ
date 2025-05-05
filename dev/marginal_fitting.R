install.packages("rstan")
library(rstan)
options(mc.cores = parallel::detectCores()-2)
rstan_options(auto_write = TRUE)

# priors. Gives a gamma with mean 1 and sd k/k^2=1/k
k = 10
alpha_a <- k # shape
beta_a <- k # rate
alpha_b <- k
beta_b <- k

N=17
Z_almost_one <- rep(0.95, N)
Z_one <- rep(1, N)

data_list <- list(
  N = N,
  Z = Z_one,
  alpha_a = alpha_a,
  beta_a = beta_a,
  alpha_b = alpha_b,
  beta_b = beta_b
)

mod <- CopulaSEJ:::load_stan_model()

fit_map <- rstan::optimizing(
  mod,
  data = data_list,
  verbose = TRUE
)
fit_map

map_a <- fit_map$par[[1]]
map_b <- fit_map$par[[2]]
# plot the results
samples <- rbeta(1000000, map_a, map_b)
density(samples) |> plot()

x <- seq(0, 1, length.out = 100)
x <- x[x > 0 & x < 1]
y <- dbeta(x, map_a, map_b)
# plot
df <- tibble(x = x, y = y)
ggplot(df, aes(x = x, y = y)) +
  geom_line() +
  labs(title = "Beta distribution with MAP estimates",
       x = "x",
       y = "Density") +
  theme_minimal() +
  ylim(0,NA)

