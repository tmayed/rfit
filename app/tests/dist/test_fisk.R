# Test suite for Fisk Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Fisk Distribution ===\n\n")

# Generate test data (Log-logistic)
# Using relationship with rlogis
test_scale <- 5
test_shape <- 5
test_data <- exp(rlogis(1000, location = log(test_scale), scale = 1/test_shape))

# 1. Test log-likelihood
cat("Testing fisk_2p_log_likelihood...\n")
ll <- fisk_2p_log_likelihood(test_data, test_scale, test_shape)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# 2. Test fit
cat("Testing fisk_2p_fit...\n")
fit <- fisk_2p_fit(test_data)
cat(sprintf("  scale: %.4f (expected ~%.4f)\n", fit$scale, test_scale))
cat(sprintf("  shape: %.4f (expected ~%.4f)\n", fit$shape, test_shape))
stopifnot(abs(fit$scale - test_scale) < 0.5)
stopifnot(abs(fit$shape - test_shape) < 0.5)

# 3. Test PDF / LogPDF
cat("Testing fisk_2p_pdf/logpdf...\n")
x <- seq(0.1, 10, length.out = 5)
dens <- fisk_2p_pdf(x, fit)
log_dens <- fisk_2p_logpdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", dens), collapse = ", ")))
stopifnot(all(dens >= 0))
stopifnot(all(abs(log(dens + 1e-100) - log_dens) < 1e-7))

# 4. Test CDF / SF / LogCDF / LogSF
cat("Testing fisk_2p_cdf/sf/logcdf/logsf...\n")
cdf_vals <- fisk_2p_cdf(x, fit)
sf_vals <- fisk_2p_sf(x, fit)
log_cdf_vals <- fisk_2p_logcdf(x, fit)
log_sf_vals <- fisk_2p_logsf(x, fit)
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", cdf_vals), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-7))
stopifnot(all(abs(log(cdf_vals + 1e-100) - log_cdf_vals) < 1e-7))
stopifnot(all(abs(log(sf_vals + 1e-100) - log_sf_vals) < 1e-7))

# 5. Test Quantile / ISF
cat("Testing fisk_2p_quantile/isf...\n")
p <- c(0.025, 0.5, 0.975)
qs <- fisk_2p_quantile(p, fit)
isfs <- fisk_2p_isf(1 - p, fit)
cat(sprintf("  Quantiles: %.4f, %.4f, %.4f\n", qs[1], qs[2], qs[3]))
stopifnot(all(diff(qs) > 0))
stopifnot(all(abs(qs - isfs) < 1e-7))

# 6. Test Mean / Var / Std / Skew / Kurtosis
cat("Testing moments/statistics...\n")
m <- fisk_2p_mean(fit)
v <- fisk_2p_var(fit)
s <- fisk_2p_std(fit)
sk <- fisk_2p_skew(fit)
kt <- fisk_2p_kurtosis(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
cat(sprintf("  Skewness: %.4f, Kurtosis: %.4f\n", sk, kt))
stopifnot(is.finite(m))
if (fit$shape > 2) {
  stopifnot(is.finite(v) && is.finite(s))
  stopifnot(all.equal(sqrt(v), s))
}

# 7. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- fisk_2p_median(fit)
int <- fisk_2p_interval(0.95, fit)
m2 <- fisk_2p_moment(2, fit)
cat(sprintf("  Median: %.4f\n", med))
cat(sprintf("  95%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(abs(med - fisk_2p_quantile(0.5, fit)) < 1e-7)
if (fit$shape > 2) {
  stopifnot(abs(m2 - (v + m^2)) < 1e-7)
}

# 8. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- fisk_2p_entropy(fit)
exp_x <- fisk_2p_expect(function(x) x, fit)
cat(sprintf("  Entropy: %.4f\n", ent))
cat(sprintf("  Expect(x): %.4f\n", exp_x))
stopifnot(is.finite(ent))
if (fit$shape > 1) {
    stopifnot(abs(exp_x - m) < 1e-4)
}

# 9. Test Truncated Fit
cat("Testing truncated fit...\n")
lower_bound <- 2
upper_bound <- 15
trunc_data <- test_data[test_data >= lower_bound & test_data <= upper_bound]
trunc_fit <- fisk_2p_fit_truncated(trunc_data, lower = lower_bound, upper = upper_bound)
cat(sprintf("  Truncated scale: %.4f, shape: %.4f\n", trunc_fit$scale, trunc_fit$shape))
stopifnot(is.finite(trunc_fit$log_likelihood))
stopifnot(trunc_fit$scale > 0 && trunc_fit$shape > 0)

cat("\n=== Fisk tests passed! ===\n")
