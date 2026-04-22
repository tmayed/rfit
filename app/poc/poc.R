# Proof of concept: fit lognormal and pareto to lognormal data

source("pkg/dist/lognormal_2p.R")
source("pkg/dist/pareto_2p.R")
source("pkg/plots/cdf_plot.R")

# Prevent Rplots.pdf from being created
if (!interactive()) {
  pdf(NULL)
}

# Set seed for reproducibility
set.seed(42)

sample_nl_fit = list(
  mu = 1,
  sigma = 0.5
)

# Generate random sample data from lognormal distribution
cat("Generating sample data from lognormal(mu=1, sigma=0.5)...\n")
# sample_data <- rlnorm(1000, meanlog = 1, sdlog = 0.5)

sample_data <- lognormal_rand(1000, fit=sample_nl_fit)

x_vals <- seq(
  min(sample_data),
  max(sample_data),
  length.out = 100
)

cat(sprintf("  Generated %d samples\n", length(sample_data)))
cat(sprintf("  Sample mean: %.4f\n", mean(sample_data)))
cat(sprintf("  Sample std: %.4f\n", sd(sample_data)))
cat("\n")

# Fit lognormal distribution
cat("=== Fitting Lognormal Distribution ===\n")
fit_lognormal <- lognormal_2p_fit(sample_data)
cat(sprintf("  mu: %.4f\n", fit_lognormal$mu))
cat(sprintf("  sigma: %.4f\n", fit_lognormal$sigma))
cat(sprintf("  log-likelihood: %.4f\n", fit_lognormal$log_likelihood))
cat("\n")

# Fit pareto distribution
cat("=== Fitting Pareto Distribution ===\n")
fit_pareto <- pareto_2p_fit(sample_data)
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

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
output_dir <- file.path(dirname(normalizePath(script_path)), "outputs")

plot_cdf_comparison(
  sample_data = sample_data,
  fit = fit_lognormal,
  dist_cdf = lognormal_cdf,
  output_dir = output_dir,
  output_file = "lognormal_cdf"
)

plot_cdf_comparison(
  sample_data = sample_data,
  fit = fit_pareto,
  dist_cdf = pareto_cdf,
  output_dir = output_dir,
  output_file = "pareto_cdf"
)
