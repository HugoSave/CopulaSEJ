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
    sigma_label = sprintf("Var = %.1f", var_val)
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

selected_variances <- c(0.1^2, 0.5^2, 1.0^2, 5.0^2, 10^2)

# Adaptive color selection based on number of variances
n_colors <- length(selected_variances)
  # For more than 6 colors, use viridis palette
colors <- viridis::viridis(n_colors, option = "plasma")

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
    sigma_label = round(sqrt(var_val), 3),
    color_group = as.factor(i)
  )

  plot_data_selected <- rbind(plot_data_selected, temp_df)
}

# Convert sigma_label to ordered factor to ensure proper legend ordering
# sort string in numeric fashion
variance_labels <- plot_data_selected$sigma_label |> unique() |> sort()
plot_data_selected$sigma_label <- factor(
  plot_data_selected$sigma_label,
  levels = variance_labels
)

p2 <- ggplot(plot_data_selected, aes(x = x, y = density, color = sigma_label)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = mode_val, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = colors, name = expression(sigma[prior]),
                     labels = parse(text = variance_labels)) +
  labs(
    # title = "Log-Normal Distributions: Selected Variance Values",
    # subtitle = "Mode = 1 (shown as dashed line)",
    x = "x",
    y = "LogNormal Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  ) +
  xlim(0,10)

print(p2)
# Save the plot
ggsave("dev/output/log_normal_prior_example.pdf", plot = p2, width = 10, height = 6, dpi = 300)

