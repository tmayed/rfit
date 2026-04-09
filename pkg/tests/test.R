# Test suite for rfit package

source("../app/dist/lognormal.R")
source("../app/dist/pareto.R")

# Set seed for reproducibility
set.seed(42)

cat("=== Testing Lognormal Distribution ===\n\n")

# Generate test data from lognormal
test_data <- rlnorm(1000, meanlog = 1, sdlog = 0.5)

# Test lognormal_log_likelihood
cat("Testing lognormal_log_likelihood...\n")
ll <- lognormal_log_likelihood(test_data)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# Test lognormal_fit
cat("Testing lognormal_fit...\n")
fit <- lognormal_fit(test_data)
cat(sprintf("  mu: %.4f (expected ~1.0)\n", fit$mu))
cat(sprintf("  sigma: %.4f (expected ~0.5)\n", fit$sigma))
cat(sprintf("  log_likelihood: %.4f\n", fit$log_likelihood))
stopifnot(abs(fit$mu - 1) < 0.1)
stopifnot(abs(fit$sigma - 0.5) < 0.1)

# Test lognormal_fit_truncated
cat("Testing lognormal_fit_truncated...\n")
trunc_data <- test_data[test_data > 1 & test_data < 10]
fit_trunc <- lognormal_fit_truncated(trunc_data, lower = 1, upper = 10)
cat(sprintf("  mu: %.4f\n", fit_trunc$mu))
cat(sprintf("  sigma: %.4f\n", fit_trunc$sigma))
stopifnot(is.finite(fit_trunc$log_likelihood))

# Test lognormal_pdf
cat("Testing lognormal_pdf...\n")
x <- seq(0.1, 10, length.out = 10)
dens <- lognormal_pdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", head(dens)), collapse = ", ")))
stopifnot(all(dens >= 0))

# Test lognormal_cdf
cat("Testing lognormal_cdf...\n")
cdf_vals <- lognormal_cdf(x, fit)
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", head(cdf_vals)), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(diff(cdf_vals) >= 0)  # Monotonically increasing

# Test lognormal_quantile
cat("Testing lognormal_quantile...\n")
p <- seq(0, 1, length.out = 10)
quantiles <- lognormal_quantile(p, fit)
cat(sprintf("  Quantiles: %s\n", paste(sprintf("%.4f", quantiles), collapse = ", ")))
stopifnot(diff(quantiles) >= 0)  # Monotonically increasing

# Test lognormal_rand
cat("Testing lognormal_rand...\n")
samples <- lognormal_rand(100, fit)
cat(sprintf("  Generated %d samples, mean: %.4f\n", length(samples), mean(samples)))
stopifnot(length(samples) == 100)

# Test lognormal_mean
cat("Testing lognormal_mean...\n")
mean_val <- lognormal_mean(fit)
cat(sprintf("  Mean: %.4f\n", mean_val))
stopifnot(is.finite(mean_val))

# Test lognormal_std
cat("Testing lognormal_std...\n")
std_val <- lognormal_std(fit)
cat(sprintf("  Std: %.4f\n", std_val))
stopifnot(is.finite(std_val))

cat("\n=== Testing Pareto Distribution ===\n\n")

# Generate test data from Pareto
pareto_data <- vapply(runif(1000), function(u) 1 / (1 - u)^(1/2.5), numeric(1))  # shape = 2.5, scale = 1

# Test pareto_log_likelihood
cat("Testing pareto_log_likelihood...\n")
ll <- pareto_log_likelihood(pareto_data)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# Test pareto_fit
cat("Testing pareto_fit...\n")
fit <- pareto_fit(pareto_data)
cat(sprintf("  shape: %.4f (expected ~2.5)\n", fit$shape))
cat(sprintf("  scale: %.4f (expected ~min(data))\n", fit$scale))
cat(sprintf("  log_likelihood: %.4f\n", fit$log_likelihood))
stopifnot(is.finite(fit$shape))
stopifnot(is.finite(fit$scale))

# Test pareto_fit_truncated
cat("Testing pareto_fit_truncated...\n")
trunc_data <- pareto_data[pareto_data > 1.5 & pareto_data < 10]
fit_trunc <- pareto_fit_truncated(trunc_data, lower = 1.5, upper = 10)
cat(sprintf("  shape: %.4f\n", fit_trunc$shape))
cat(sprintf("  scale: %.4f\n", fit_trunc$scale))
stopifnot(is.finite(fit_trunc$log_likelihood))

# Test pareto_pdf
cat("Testing pareto_pdf...\n")
x <- seq(fit$scale, fit$scale * 5, length.out = 10)
dens <- pareto_pdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", head(dens)), collapse = ", ")))
stopifnot(all(dens >= 0))

# Test pareto_cdf
cat("Testing pareto_cdf...\n")
cdf_vals <- pareto_cdf(x, fit)
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", head(cdf_vals)), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(diff(cdf_vals) >= 0)  # Monotonically increasing

# Test pareto_quantile
cat("Testing pareto_quantile...\n")
p <- seq(0, 0.99, length.out = 10)
quantiles <- pareto_quantile(p, fit)
cat(sprintf("  Quantiles: %s\n", paste(sprintf("%.4f", quantiles), collapse = ", ")))
stopifnot(diff(quantiles) >= 0)  # Monotonically increasing

# Test pareto_rand
cat("Testing pareto_rand...\n")
samples <- pareto_rand(100, fit)
cat(sprintf("  Generated %d samples, mean: %.4f\n", length(samples), mean(samples)))
stopifnot(length(samples) == 100)

# Test pareto_mean
cat("Testing pareto_mean...\n")
mean_val <- pareto_mean(fit)
cat(sprintf("  Mean: %.4f\n", mean_val))
stopifnot(is.finite(mean_val) || is.infinite(mean_val))

# Test pareto_std
cat("Testing pareto_std...\n")
std_val <- pareto_std(fit)
cat(sprintf("  Std: %.4f\n", std_val))
stopifnot(is.finite(std_val) || is.infinite(std_val))

cat("\n=== All tests passed! ===\n")
