# Test suite for Kappa 4 Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Kappa 4 Distribution ===\n\n")

# Generate test data
# xi=10, alpha=5, k=0.1, h=0.2
params <- list(xi = 10, alpha = 5, k = 0.1, h = 0.2)
test_data <- kappa4_4p_rand(1000, params)

# 1. Test log-likelihood
cat("Testing kappa4_4p_log_likelihood...\n")
ll <- kappa4_4p_log_likelihood(test_data, 10, 5, 0.1, 0.2)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# 2. Test fit
cat("Testing kappa4_4p_fit...\n")
fit <- kappa4_4p_fit(test_data)
cat(sprintf("  xi: %.4f\n", fit$xi))
cat(sprintf("  alpha: %.4f\n", fit$alpha))
cat(sprintf("  k: %.4f\n", fit$k))
cat(sprintf("  h: %.4f\n", fit$h))
stopifnot(fit$log_likelihood > -4000) 

# 3. Test PDF / LogPDF
cat("Testing kappa4_4p_pdf/logpdf...\n")
x <- seq(min(test_data), max(test_data), length.out = 10)
dens <- kappa4_4p_pdf(x, fit)
log_dens <- kappa4_4p_logpdf(x, fit)
cat(sprintf("  Densities (head): %s\n", paste(sprintf("%.4f", head(dens, 3)), collapse = ", ")))
stopifnot(all(dens >= 0))
stopifnot(all.equal(log(dens), log_dens, tolerance = 1e-7))

# 4. Test CDF / SF / LogCDF / LogSF
cat("Testing kappa4_4p_cdf/sf/logcdf/logsf...\n")
cdf_vals <- kappa4_4p_cdf(x, fit)
sf_vals <- kappa4_4p_sf(x, fit)
log_cdf_vals <- kappa4_4p_logcdf(x, fit)
log_sf_vals <- kappa4_4p_logsf(x, fit)
cat(sprintf("  CDF values (head): %s\n", paste(sprintf("%.4f", head(cdf_vals, 3)), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all.equal(cdf_vals + sf_vals, rep(1, length(x))))
stopifnot(all.equal(log(cdf_vals), log_cdf_vals, tolerance = 1e-7))
stopifnot(all.equal(log(sf_vals), log_sf_vals, tolerance = 1e-7))

# 5. Test Quantile / ISF
cat("Testing kappa4_4p_quantile/isf...\n")
p <- seq(0.1, 0.9, length.out = 5)
qs <- kappa4_4p_quantile(p, fit)
isfs <- kappa4_4p_isf(1 - p, fit)
cat(sprintf("  Quantiles: %s\n", paste(sprintf("%.4f", qs), collapse = ", ")))
stopifnot(all(diff(qs) >= 0))
stopifnot(all.equal(qs, isfs))

# 6. Test Mean / Var / Std / Skew / Kurtosis
cat("Testing moments/statistics...\n")
m <- kappa4_4p_mean(fit)
v <- kappa4_4p_var(fit)
s <- kappa4_4p_std(fit)
sk <- kappa4_4p_skew(fit)
kt <- kappa4_4p_kurtosis(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
cat(sprintf("  Skewness: %.4f, Kurtosis: %.4f\n", sk, kt))
stopifnot(is.finite(m) && is.finite(v) && is.finite(s))
stopifnot(all.equal(sqrt(v), s, tolerance = 1e-10))

# 7. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- kappa4_4p_median(fit)
int <- kappa4_4p_interval(0.95, fit)
m2 <- kappa4_4p_moment(2, fit)
cat(sprintf("  Median: %.4f\n", med))
cat(sprintf("  95%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(all.equal(med, kappa4_4p_quantile(0.5, fit)))
stopifnot(m2 > 0)

# 8. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- kappa4_4p_entropy(fit)
exp_x <- kappa4_4p_expect(function(x) x, fit)
cat(sprintf("  Entropy: %.4f\n", ent))
cat(sprintf("  Expect(x): %.4f\n", exp_x))
stopifnot(is.finite(ent))
# Expect(x) should be close to Mean
stopifnot(abs(exp_x - m) / m < 0.05)

# 9. Test k=0 stability
cat("Testing k=0 stability...\n")
fit0 <- list(xi = 10, alpha = 5, k = 0, h = 0.2)
q0 <- kappa4_4p_quantile(0.5, fit0)
c0 <- kappa4_4p_cdf(10, fit0)
cat(sprintf("  k=0 Quantile(0.5): %.4f\n", q0))
stopifnot(is.finite(q0) && is.finite(c0))

cat("\n=== All Kappa 4 tests passed! ===\n")
