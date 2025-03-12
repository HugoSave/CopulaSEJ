# Function to calculate the median of a Beta distribution
beta_median <- function(shape1, shape2) {
  return(qbeta(0.5, shape1, shape2))
}

# Function to sample from a Beta distribution
sample_beta <- function(shape1, shape2, n_samples) {
  return(rbeta(n_samples, shape1, shape2))
}

# Example usage
shape1 <- 3  # Alpha parameter
shape2 <- 1.3  # Beta parameter
n_samples <- 10000  # Number of samples

# Calculate median
median_value <- beta_median(shape1, shape2)
cat("Median of Beta(", shape1, ",", shape2, "):", median_value, "\n")

# Generate samples
samples <- sample_beta(shape1, shape2, n_samples)

error_paper <- function(m, q) {
  q - m
}

error_normal <- function(m, q) {
  m - q
}

rel_error <- function(m, q) {
  (m - q) / q
}

# plot all histograms
par(mfrow=c(2,2))
hist(samples, breaks = 30, col = "skyblue", main = "Histogram of Beta Samples", xlab = "Value", probability = TRUE)
lines(density(samples), col = "red", lwd = 2)
abline(v = median_value, col = "blue", lwd = 2)
hist(error_paper(median_value, samples), breaks = 30, col = "skyblue", main = "q-median", xlab = "Error", probability = TRUE)
abline(v = 0, col = "black", lwd = 2)
hist(error_normal(median_value, samples), breaks = 30, col = "skyblue", main = "median-q", xlab = "Error", probability = TRUE)
abline(v = 0, col = "black", lwd = 2)
# hist(rel_error(median_value, samples), breaks = 30, col = "skyblue", main = "(m-q)/q", xlab = "Error", probability = TRUE)

