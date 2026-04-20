# Test suite for Lognormal Distribution

source("../../pkg/dist/lognormal.R")

# Set seed for reproducibility
set.seed(42)

cat("=== Testing Lognormal Distribution ===\n\n")

# Generate test data from lognormal
test_data <- rlnorm(1000, meanlog = 1, sdlog = 0.5)

# Test lognormal_log_likelihood
cat("Testing lognormal_log_likelihood...\n")
ll <- lognormal_log_likelihood(test_data, 1, 0.5)
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

# Test lognormal_sf
cat("Testing lognormal_sf...\n")
sf <- lognormal_sf(x, fit)
cat(sprintf("  SF: %s\n", paste(sprintf("%.4f", head(sf)), collapse = ", ")))
stopifnot(all(abs(sf - (1 - cdf_vals)) < 1e-10))

# Test lognormal_quantile
cat("Testing lognormal_quantile...\n")
p <- seq(0.01, 0.99, length.out = 10)
quantiles <- lognormal_quantile(p, fit)
cat(sprintf("  Quantiles: %s\n", paste(sprintf("%.4f", quantiles), collapse = ", ")))
stopifnot(diff(quantiles) >= 0)  # Monotonically increasing

# Test lognormal_isf
cat("Testing lognormal_isf...\n")
isf <- lognormal_isf(p, fit)
cat(sprintf("  ISF: %s\n", paste(sprintf("%.4f", isf), collapse = ", ")))
stopifnot(all(abs(lognormal_cdf(isf, fit) - (1 - p)) < 1e-10))

# Test Logs
cat("Testing Log PDF/CDF/SF...\n")
lpdf <- lognormal_logpdf(x, fit)
lcdf <- lognormal_logcdf(x, fit)
lsf <- lognormal_logsf(x, fit)
stopifnot(all(abs(lpdf - log(lognormal_pdf(x, fit))) < 1e-10))
stopifnot(all(abs(lcdf - log(lognormal_cdf(x, fit))) < 1e-10))
stopifnot(all(abs(lsf - log(lognormal_sf(x, fit))) < 1e-10))

# Test lognormal_rand
cat("Testing lognormal_rand...\n")
samples <- lognormal_rand(100, fit)
cat(sprintf("  Generated %d samples, mean: %.4f\n", length(samples), mean(samples)))
stopifnot(length(samples) == 100)

# Test lognormal_mean/var/std
cat("Testing lognormal_mean/var/std...\n")
mean_val <- lognormal_mean(fit)
var_val <- lognormal_var(fit)
std_val <- lognormal_std(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", mean_val, var_val, std_val))
stopifnot(abs(std_val - sqrt(var_val)) < 1e-10)

# Test lognormal_moment
cat("Testing lognormal_moment...\n")
m1 <- lognormal_moment(1, fit)
m2 <- lognormal_moment(2, fit)
cat(sprintf("  1st Moment: %.4f, 2nd Moment: %.4f\n", m1, m2))
stopifnot(abs(m1 - mean_val) < 1e-10)

# Test lognormal_skew/kurtosis
cat("Testing lognormal_skew/kurtosis...\n")
sk <- lognormal_skew(fit)
ku <- lognormal_kurtosis(fit)
cat(sprintf("  Skew: %.4f, Kurtosis: %.4f\n", sk, ku))
stopifnot(is.finite(sk) && is.finite(ku))

# Test lognormal_median
cat("Testing lognormal_median...\n")
med <- lognormal_median(fit)
cat(sprintf("  Median: %.4f\n", med))
stopifnot(abs(med - exp(fit$mu)) < 1e-10)

# Test lognormal_interval
cat("Testing lognormal_interval...\n")
int <- lognormal_interval(0.95, fit)
cat(sprintf("  95%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(abs(lognormal_cdf(int[2], fit) - lognormal_cdf(int[1], fit) - 0.95) < 1e-10)

# Test lognormal_entropy
cat("Testing lognormal_entropy...\n")
ent <- lognormal_entropy(fit)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))

# Test lognormal_expect
cat("Testing lognormal_expect...\n")
e_x <- lognormal_expect(function(x) x, fit)
cat(sprintf("  E[x]: %.4f\n", e_x))
stopifnot(abs(e_x - mean_val) < 1e-5)

cat("\n=== Lognormal tests passed! ===\n")
