# Exhaustive Test suite for Johnson SU Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Johnson SU Distribution ===\n\n")

# Parameters: gamma=0.5, delta=1.5, xi=10, lambda=2
fit_true <- list(gamma = 0.5, delta = 1.5, xi = 10, lambda = 2)

# 1. Test Random Generation
cat("Testing johnsonsu_4p_rand...\n")
test_data <- johnsonsu_4p_rand(1000, fit_true)
stopifnot(length(test_data) == 1000)
cat(sprintf("  Mean of sample: %.4f\n", mean(test_data)))

# 2. Test Log-likelihood
cat("Testing johnsonsu_4p_log_likelihood...\n")
ll <- johnsonsu_4p_log_likelihood(test_data, 0.5, 1.5, 10, 2)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# 3. Test Fit
cat("Testing johnsonsu_4p_fit...\n")
fit_hat <- johnsonsu_4p_fit(test_data)
cat(sprintf("  Fitted gamma:  %.4f\n", fit_hat$gamma))
cat(sprintf("  Fitted delta:  %.4f\n", fit_hat$delta))
cat(sprintf("  Fitted xi:     %.4f\n", fit_hat$xi))
cat(sprintf("  Fitted lambda: %.4f\n", fit_hat$lambda))
stopifnot(abs(fit_hat$gamma - 0.5) < 0.3)
stopifnot(abs(fit_hat$delta - 1.5) < 0.3)

# 4. Test Truncated Fit
cat("Testing johnsonsu_4p_fit_truncated...\n")
# Truncate data to [5, 15]
fit_t <- johnsonsu_4p_fit_truncated(test_data, lower = 5, upper = 15)
cat(sprintf("  Truncated Log-likelihood: %.4f\n", fit_t$log_likelihood))
stopifnot(is.finite(fit_t$log_likelihood))

# 5. Test PDF / LogPDF
cat("Testing johnsonsu_4p_pdf/logpdf...\n")
x <- seq(5, 15, length.out = 5)
dens <- johnsonsu_4p_pdf(x, fit_hat)
log_dens <- johnsonsu_4p_logpdf(x, fit_hat)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", dens), collapse = ", ")))
stopifnot(all(dens > 0))
stopifnot(all.equal(log(dens), log_dens))

# 6. Test CDF / SF / LogCDF / LogSF
cat("Testing johnsonsu_4p_cdf/sf/logcdf/logsf...\n")
cdf_vals <- johnsonsu_4p_cdf(x, fit_hat)
sf_vals <- johnsonsu_4p_sf(x, fit_hat)
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", cdf_vals), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all.equal(cdf_vals + sf_vals, rep(1, length(x))))
stopifnot(all.equal(log(cdf_vals), johnsonsu_4p_logcdf(x, fit_hat)))
stopifnot(all.equal(log(sf_vals), johnsonsu_4p_logsf(x, fit_hat)))

# 7. Test Quantile / ISF
cat("Testing johnsonsu_4p_quantile/isf...\n")
p <- c(0.25, 0.5, 0.75)
qs <- johnsonsu_4p_quantile(p, fit_hat)
isfs <- johnsonsu_4p_isf(1 - p, fit_hat)
cat(sprintf("  Quantiles: %s\n", paste(sprintf("%.4f", qs), collapse = ", ")))
stopifnot(all.equal(qs, isfs))
stopifnot(all.equal(johnsonsu_4p_cdf(qs, fit_hat), p))

# 8. Test Mean / Var / Std / Skew / Kurtosis
cat("Testing moments/statistics...\n")
m <- johnsonsu_4p_mean(fit_hat)
v <- johnsonsu_4p_var(fit_hat)
s <- johnsonsu_4p_std(fit_hat)
sk <- johnsonsu_4p_skew(fit_hat)
kt <- johnsonsu_4p_kurtosis(fit_hat)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
cat(sprintf("  Skewness: %.4f, Kurtosis: %.4f\n", sk, kt))
stopifnot(is.finite(m) && is.finite(v) && is.finite(s))
stopifnot(all.equal(sqrt(v), s))

# 9. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- johnsonsu_4p_median(fit_hat)
int <- johnsonsu_4p_interval(0.95, fit_hat)
m1 <- johnsonsu_4p_moment(1, fit_hat)
cat(sprintf("  Median: %.4f\n", med))
stopifnot(all.equal(med, johnsonsu_4p_quantile(0.5, fit_hat)))
stopifnot(all.equal(m1, m, tolerance = 1e-5))

# 10. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- johnsonsu_4p_entropy(fit_hat)
exp_x <- johnsonsu_4p_expect(function(x) x, fit_hat)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))
stopifnot(all.equal(exp_x, m, tolerance = 1e-5))

cat("\n=== All Johnson SU tests passed! ===\n")
