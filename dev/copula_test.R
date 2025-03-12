# Copula package
library(copula)
# ggplot2
library(ggplot2)
set.seed(235)

joeCop <- joeCopula(param=32, dim =2)
joeCop
s<-dMvdc(c(0.99999,0.9999), mvdc(copula = joeCop, margins = c("norm", "norm"), paramMargins = list(list(mean = 2, sd=3), list(mean = 3, sd=2))))


# Generate a bivariate normal copula with rho = 0.7
normal <- normalCopula(param = 0.7, dim = 2)
# Generate a bivariate t-copula with rho = 0.8 and df = 2
stc <- tCopula(param = 0.8, dim = 2, df = 2)

## Archimedan copulas:
# Build a Frank, a Gumbel and a Clayton copula
frank <- frankCopula(dim = 2, param = 8)
gumbel <- gumbelCopula(dim = 3, param = 5.6)
clayton <- claytonCopula(dim = 4, param = 19)

# Print information on the Frank copula
print(frank)


# Select the copula
cp <- claytonCopula(param = c(3.4), dim = 2)

# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
multivariate_dist <- mvdc(copula = cp,
                          margins = c("norm", "t"),
                          paramMargins = list(list(mean = 2, sd=3),
                                              list(df = 2)) )

print(multivariate_dist)

samples <- rMvdc(2000, multivariate_dist)

as.data.frame(samples) |> ggplot(aes(x = V1, y = V2)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  labs(title = "Clayton Copula with Normal and t Marginals") +
  theme(plot.title = element_text(hjust = 0.5))

densityplot(samples[,1], display = "persp", col = "lightblue", main = "Normal Margin")
densityplot(samples[,2], display = "persp", col = "lightblue", main = "t Margin")

# pseudo observations:
pseudo_obs <- pobs(samples)
cor(pseudo_obs)
as.data.frame(pseudo_obs) |> ggplot(aes(x = V1, y = V2)) +
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  geom_point(alpha=0.6) +
  theme_minimal() +
  labs(title = "Pseudo obs - Clayton Copula with Normal and t Marginals") +
  theme(plot.title = element_text(hjust = 0.5))

# Try to fit copula to data
cop_model = claytonCopula(dim=2)
#cop_model = frankCopula(dim=2)
#cop_model = joeCopula(dim=2)
fit <- fitCopula(cop_model, pseudo_obs, method = "mpl")

fitted_copula <- claytonCopula(param = coef(fit))
tau(fitted_copula)
# So kendall tau is invariant of the pseudo transform but not pearson 
# correlation
cor(pseudo_obs, method = "kendall")
cor(pseudo_obs)
cor(samples, method = "kendall")
cor(samples)
