# Test suite for Weibull Distribution
source("../../app/rfit.R")

set.seed(42)

cat("=== Testing Weibull Distribution ===\n\n")

# Generate test data
test_data <- rweibull(1000, shape = 2, scale = 5)

# 1. Test log-likelihood
cat("Testing weibull_log_likelihood...\n")
ll <- weibull_log_likelihood(test_data, 2, 5)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# 2. Test fit
cat("Testing weibull_fit...\n")
fit <- weibull_fit(test_data)
cat(sprintf("  shape: %.4f (expected ~2.0)\n", fit$shape))
cat(sprintf("  scale: %.4f (expected ~5.0)\n", fit$scale))
stopifnot(abs(fit$shape - 2) < 0.2)
stopifnot(abs(fit$scale - 5) < 0.2)

# 3. Test PDF / LogPDF
cat("Testing weibull_pdf/logpdf...\n")
x <- seq(0.1, 10, length.out = 5)
dens <- weibull_pdf(x, fit)
log_dens <- weibull_logpdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", dens), collapse = ", ")))
stopifnot(all(dens >= 0))
stopifnot(all.equal(log(dens), log_dens))

# 4. Test CDF / SF / LogCDF / LogSF
cat("Testing weibull_cdf/sf/logcdf/logsf...\n")
cdf_vals <- weibull_cdf(x, fit)
sf_vals <- weibull_sf(x, fit)
log_cdf_vals <- weibull_logcdf(x, fit)
log_sf_vals <- weibull_logsf(x, fit)
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", cdf_vals), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all.equal(cdf_vals + sf_vals, rep(1, length(x))))
stopifnot(all.equal(log(cdf_vals), log_cdf_vals))
stopifnot(all.equal(log(sf_vals), log_sf_vals))

# 5. Test Quantile / ISF
cat("Testing weibull_quantile/isf...\n")
p <- c(0.025, 0.5, 0.975)
qs <- weibull_quantile(p, fit)
isfs <- weibull_isf(1 - p, fit)
cat(sprintf("  Quantiles: %.4f, %.4f, %.4f\n", qs[1], qs[2], qs[3]))
stopifnot(all(diff(qs) > 0))
stopifnot(all.equal(qs, isfs))

# 6. Test Mean / Var / Std / Skew / Kurtosis
cat("Testing moments/statistics...\n")
m <- weibull_mean(fit)
v <- weibull_var(fit)
s <- weibull_std(fit)
sk <- weibull_skew(fit)
kt <- weibull_kurtosis(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
cat(sprintf("  Skewness: %.4f, Kurtosis: %.4f\n", sk, kt))
stopifnot(is.finite(m) && is.finite(v) && is.finite(s))
stopifnot(all.equal(sqrt(v), s))

# 7. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- weibull_median(fit)
int <- weibull_interval(0.95, fit)
m2 <- weibull_moment(2, fit)
cat(sprintf("  Median: %.4f\n", med))
cat(sprintf("  95%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(all.equal(med, weibull_quantile(0.5, fit)))
stopifnot(all.equal(m2, v + m^2))

# 8. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- weibull_entropy(fit)
exp_x <- weibull_expect(function(x) x, fit)
cat(sprintf("  Entropy: %.4f\n", ent))
cat(sprintf("  Expect(x): %.4f\n", exp_x))
stopifnot(is.finite(ent))
stopifnot(all.equal(exp_x, m, tolerance = 1e-5))

cat("\n=== Weibull tests passed! ===\n")
