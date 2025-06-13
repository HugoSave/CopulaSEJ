# Log-Normal Distribution: Fixed Median with Varying Variance

library(ggplot2)
library(dplyr)

# Function to get log-normal parameters from mode and variance
lognormal_from_mode_var <- function(mode, variance) {

  # Define the function u^4 - u^3 - v/M^2 = 0
  # where u = exp(sigma^2)
  c_val <- variance / mode^2

  root_function <- function(u) {
    u^4 - u^3 - c_val
  }

  # Find the root using uniroot
  # We know u > 1 for positive variance, so search in (1, 10)
  # If that fails, expand the search range
  tryCatch({
    root_result <- uniroot(root_function, interval = c(1.001, 10))
    u <- root_result$root
  }, error = function(e) {
    # Try wider range if initial search fails
    print(c_val)
    root_result <- uniroot(root_function, interval = c(1.001, 100))
    u <- root_result$root
  })

  # Calculate sigma^2 and mu
  sigma_sq <- log(u)
  mu <- log(mode) + sigma_sq

  # Verify our solution by computing actual mode and variance
  actual_mode <- exp(mu - sigma_sq)
  actual_variance <- exp(2*mu + sigma_sq) * (exp(sigma_sq) - 1)

  return(list(
    mu = mu,
    sigma = sqrt(sigma_sq),
    sigma_sq = sigma_sq,
    u = u,
    actual_mode = actual_mode,
    actual_variance = actual_variance,
    mode_error = abs(actual_mode - mode),
    variance_error = abs(actual_variance - variance)
  ))
}

# Create data for multiple variance values
mode_val <- 1
variance_values <- seq(0.1, 5, by = 0.3)

# Generate data for plotting
plot_data <- data.frame()

for (var_val in variance_values) {
  params <- lognormal_from_mode_var(mode_val, var_val)

  # Create x values - need wider range for higher variances
  x_max <- mode_val + 6*sqrt(var_val)
  x <- seq(0.01, x_max, length.out = 500)

  # Calculate densities
  y <- dlnorm(x, meanlog = params$mu, sdlog = params$sigma)

  # Add to data frame
  temp_df <- data.frame(
    x = x,
    density = y,
    variance = var_val,
    variance_label = sprintf("Var = %.1f", var_val)
  )

  plot_data <- rbind(plot_data, temp_df)
}

# Create the plot
p1 <- ggplot(plot_data, aes(x = x, y = density, color = variance)) +
  geom_line(size = 0.8) +
  geom_vline(xintercept = mode_val, linetype = "dashed", color = "black", size = 1) +
  scale_color_viridis_c(name = "Variance", option = "plasma") +
  labs(
    title = "Log-Normal Distributions with Fixed Mode = 1",
    subtitle = "Variance ranging from 0.1 to 5.0",
    x = "x",
    y = "Density"
  ) +
  annotate("text", x = mode_val, y = max(plot_data$density) * 0.95,
           label = "Mode = 1", hjust = -0.1, size = 4) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  )

print(p1)

# Create a second plot showing selected variances more clearly
# EASILY CUSTOMIZABLE: Change these values as needed
selected_variances <- c(0.1, 0.5, 1.0, 2.0, 3.5, 5.0, 100)

# Adaptive color selection based on number of variances
n_colors <- length(selected_variances)
if (n_colors <= 2) {
  colors <- c("#440154", "#fde725")[1:n_colors]
} else if (n_colors <= 6) {
  colors <- c("#440154", "#31688e", "#35b779", "#fde725", "#ff6b35", "#d73027")[1:n_colors]
} else {
  # For more than 6 colors, use viridis palette
  colors <- viridis::viridis(n_colors, option = "plasma")
}

plot_data_selected <- data.frame()

# Sort selected_variances to ensure proper ordering
selected_variances <- sort(selected_variances)

for (i in 1:length(selected_variances)) {
  var_val <- selected_variances[i]
  params <- lognormal_from_mode_var(mode_val, var_val)

  x_max <- mode_val + 5*sqrt(var_val)
  x <- seq(0.01, x_max, length.out = 500)
  y <- dlnorm(x, meanlog = params$mu, sdlog = params$sigma)

  temp_df <- data.frame(
    x = x,
    density = y,
    variance = var_val,
    variance_label = paste0(var_val),
    color_group = as.factor(i)
  )

  plot_data_selected <- rbind(plot_data_selected, temp_df)
}

# Convert variance_label to ordered factor to ensure proper legend ordering
variance_labels <- paste0(sort(selected_variances))
plot_data_selected$variance_label <- factor(
  plot_data_selected$variance_label,
  levels = variance_labels
)

p2 <- ggplot(plot_data_selected, aes(x = x, y = density, color = variance_label)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = mode_val, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = colors, name = expression(sigma[prior]^2),
                     labels = parse(text = variance_labels)) +
  labs(
    title = "Log-Normal Distributions: Selected Variance Values",
    subtitle = "Mode = 1 (shown as dashed line)",
    x = "x",
    y = "Denssity"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  ) +
  xlim(0,10)

print(p2)


# Show the parameters for selected variances
cat("Parameters for Log-Normal with Mode = 1:\n")
cat(paste(rep("=", 50), collapse="") + "\n")

for (var_val in selected_variances) {
  params <- lognormal_from_mode_var(mode_val, var_val)
  cat(sprintf("Variance = %.1f: μ = %.4f, σ = %.4f\n",
              var_val, params$mu, params$sigma))
}

# Create a plot showing how parameters change with variance
param_data <- data.frame()
variance_range <- seq(0.1, 5, by = 0.1)

for (var_val in variance_range) {
  params <- lognormal_from_mode_var(mode_val, var_val)
  param_data <- rbind(param_data, data.frame(
    variance = var_val,
    mu = params$mu,
    sigma = params$sigma
  ))
}

# Plot both mu and sigma
library(tidyr)
param_long <- param_data %>%
  pivot_longer(cols = c(mu, sigma), names_to = "parameter", values_to = "value")

p3 <- ggplot(param_long, aes(x = variance, y = value, color = parameter)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("mu" = "blue", "sigma" = "red"),
                     labels = c("μ", "σ")) +
  labs(
    title = "Log-Normal Parameters vs Variance",
    subtitle = "With fixed mode = 1",
    x = "Variance",
    y = "Parameter Value",
    color = "Parameter"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

print(p3)

# Summary statistics table
summary_stats <- data.frame()

for (var_val in selected_variances) {
  params <- lognormal_from_mode_var(mode_val, var_val)

  # Calculate mean and median
  mean_val <- exp(params$mu + params$sigma_sq/2)
  median_val <- exp(params$mu)

  summary_stats <- rbind(summary_stats, data.frame(
    Variance = var_val,
    Mu = round(params$mu, 4),
    Sigma = round(params$sigma, 4),
    Mean = round(mean_val, 4),
    Median = round(median_val, 4),
    Mode = mode_val,
    CV = round(sqrt(var_val)/mode_val, 4)
  ))
}

cat("\n\nSummary Statistics:\n")
print(summary_stats)
