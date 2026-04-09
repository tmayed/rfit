# Proof of concept: fit lognormal and pareto to lognormal data

source("app/dist/lognormal.R")
source("app/dist/pareto.R")

# Set seed for reproducibility
set.seed(42)

# Generate random sample data from lognormal distribution
cat("Generating sample data from lognormal(mu=1, sigma=0.5)...\n")
sample_data <- rlnorm(1000, meanlog = 1, sdlog = 0.5)
cat(sprintf("  Generated %d samples\n", length(sample_data)))
cat(sprintf("  Sample mean: %.4f\n", mean(sample_data)))
cat(sprintf("  Sample std: %.4f\n", sd(sample_data)))
cat("\n")

# Fit lognormal distribution
cat("=== Fitting Lognormal Distribution ===\n")
fit_lognormal <- lognormal_fit(sample_data)
cat(sprintf("  mu: %.4f\n", fit_lognormal$mu))
cat(sprintf("  sigma: %.4f\n", fit_lognormal$sigma))
cat(sprintf("  log-likelihood: %.4f\n", fit_lognormal$log_likelihood))
cat("\n")

# Fit pareto distribution
cat("=== Fitting Pareto Distribution ===\n")
fit_pareto <- pareto_fit(sample_data)
cat(sprintf("  shape: %.4f\n", fit_pareto$shape))
cat(sprintf("  scale: %.4f\n", fit_pareto$scale))
cat(sprintf("  log-likelihood: %.4f\n", fit_pareto$log_likelihood))
cat("\n")

# Compare log-likelihoods
cat("=== Comparison ===\n")
cat(sprintf("  Lognormal log-likelihood: %.4f\n", fit_lognormal$log_likelihood))
cat(sprintf("  Pareto log-likelihood: %.4f\n", fit_pareto$log_likelihood))
cat("\n")

if (fit_lognormal$log_likelihood > fit_pareto$log_likelihood) {
  cat("Lognormal provides a better fit (higher log-likelihood).\n")
} else {
  cat("Pareto provides a better fit (higher log-likelihood).\n")
}
