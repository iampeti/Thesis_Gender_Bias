# Load necessary library
library(ggplot2)

# Create sample data
data <- data.frame(
  genotype = rep(c("level1", "level2"), each = 5),
  expression = c(rnorm(5, mean = 2, sd = 0.1), rnorm(5, mean = 5, sd = 0.1))
)

# Define beta coefficients
beta_1 <- mean(data$expression[data$genotype == "level1"])
beta_2 <- mean(data$expression[data$genotype == "level2"])

# Create the plot
ggplot(data, aes(x = genotype, y = expression)) +
  geom_point(size = 3, color = "black") +  # Plot points
  geom_hline(yintercept = beta_1, color = "green", size = 1.5) +  # β₀ line
  geom_hline(yintercept = beta_2, color = "red", size = 1.5) +  # β₁ line
  labs(x = "factor", y = "expression") +
  theme_minimal(base_size = 14) +
  annotate("text", x = 0.5, y = beta_1, label = expression(beta[1]), color = "black", parse = TRUE, size = 5) +
  annotate("text", x = 0.5, y = beta_2, label = expression(beta[2]), color = "black", parse = TRUE, size = 5)
