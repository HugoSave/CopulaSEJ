library(ICSsmoothing)

expert_belief_1 = list(cumprob = c( 0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                       fractile = c(-5.5, -3.5, -1.5, -0.7, 0.5, 2.5, 4.5))

expert_belief_2 = list(cumprob = c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                       fractile = c(-2.7, -1.8, -0.9, -0.5, 0.5, 1.8, 3.1))
                       
expert_belief_3 = list(cumprob = c(0.0, 0.05, 0.17, 0.5, 0.83, 0.95, 1),
                       fractile = c(-3.8, -2.0, -0.5, 1.0, 2.5, 5.0, 7.8))

source("Functions/linear_distribution_estimation.R")
source("Functions/create_copula_posterior.R")
source("Functions/cubic_spline_distribution_estimation.R")

# p_1 <- linear_distribution_estimation(expert_belief_1$fractile, expert_belief_1$cumprob)
# p_2 <- linear_distribution_estimation(expert_belief_2$fractile, expert_belief_2$cumprob)
# p_3 <- linear_distribution_estimation(expert_belief_3$fractile, expert_belief_3$cumprob)
# p_cdf_1 <- splinefun(expert_belief_1$fractile, expert_belief_1$cumprob, method = "monoH.FC")
# p_cdf_2 <- splinefun(expert_belief_2$fractile, expert_belief_2$cumprob, method = "monoH.FC")
# p_cdf_3 <- splinefun(expert_belief_3$fractile, expert_belief_3$cumprob, method = "monoH.FC")

p_1 <- cubic_poly_distribution_estimation_cics(expert_belief_1$fractile, expert_belief_1$cumprob)
p_2 <- cubic_poly_distribution_estimation_cics(expert_belief_2$fractile, expert_belief_2$cumprob)
p_3 <- cubic_poly_distribution_estimation_cics(expert_belief_3$fractile, expert_belief_3$cumprob)

cics_explicit(expert_belief_3$fractile, expert_belief_3$cumprob, c(0,0))
cics_explicit(expert_belief_2$fractile, expert_belief_2$cumprob, c(0,0))


# plot pdf of p_1
x = seq(-6, 8, by = 0.01)
# plot over p_2
plot(x, p_2$pdf(x), type = "l", col = "blue")
lines(x, p_3$pdf(x), type = "l", col = "green")
lines(x, p_1$pdf(x), type = "l", col = "red", xlab = "x", ylab = "Density")

# plot cdf
x = seq(-2.7, 3.1, by = 0.01)
plot(x, p_3$cdf(x), type = "l", col = "red", xlab = "x", ylab = "cdf")
plot(x, p_2$cdf(x), type = "l", col = "red", xlab = "x", ylab = "cdf")

# fit cubic spline to data
#p_cdf_spline_1 <- splinefun(expert_belief_1$fractile, expert_belief_1$cumprob, method = "monoH.FC")
#cics_explicit(expert_belief_1$fractile, expert_belief_1$cumprob, c(0,0))
#
#p_cdf_H_1 <- splinefunH(expert_belief_1$fractile, expert_belief_1$cumprob, m=c(0,NA,NA,NA,NA,NA,0))

## Define coupla
tau_paper = 0.4
alpha = exp(-4.161) # Also from paper
frank_cop <- frankCopula(dim = 3, param = 4.161)

posterior <- create_copula_posterior_numerical_integration(frank_cop, list(p_1, p_2, p_3))

x = seq(-3,3, by = 0.01)
plot(x, posterior$DM(x), type = "l", col = "red", xlab = "Fractile", ylab = "Density Marginal")

tau_paper = 0.5
alpha = exp(-5.737) # Also from paper
frank_cop <- frankCopula(dim = 3, param = 5.737)

posterior <- create_copula_posterior(frank_cop, list(p_1, p_2, p_3))
lines(x, posterior$DM(x), type = "l", col = "blue", xlab = "Fractile", ylab = "Density Marginal")

tau_paper = 0.6
alpha = exp(-7.930) # Also from paper
frank_cop <- frankCopula(dim = 3, param = 7.930)

posterior <- create_copula_posterior(frank_cop, list(p_1, p_2, p_3))
lines(x, posterior$DM(x), type = "l", col = "green", xlab = "Fractile", ylab = "Density Marginal")

tau_paper = 0.9
alpha = exp(-38.333) # Also from paper
frank_cop <- frankCopula(dim = 3, param = 38.333)

posterior <- create_copula_posterior(frank_cop, list(p_1, p_2, p_3))
lines(x, posterior$DM(x), type = "l", col = "green", xlab = "Fractile", ylab = "Density Marginal")

tau_paper = 0.0
alpha = exp(-0) # Also from paper
frank_cop <- frankCopula(dim = 3, param = alpha)

posterior <- create_copula_posterior(frank_cop, list(p_1, p_2, p_3))
lines(x, posterior$DM(x), type = "l", col = "green", xlab = "Fractile", ylab = "Density Marginal")

