# Test suite for Johnson SL Distribution
source("../../app/rfit.R")

set.seed(42)

cat("=== Testing Johnson SL Distribution ===\n\n")

# Parameters: gamma=0.5, delta=1.5, xi=10
fit_true <- list(gamma = 0.5, delta = 1.5, xi = 10)

# 1. Test Random Generation
cat("Testing johnsonsl_rand...\n")
test_data <- johnsonsl_rand(1000, fit_true)
stopifnot(all(test_data > 10))
cat(sprintf("  Mean of sample: %.4f\n", mean(test_data)))

# 2. Test Fit
cat("Testing johnsonsl_fit...\n")
fit_hat <- johnsonsl_fit(test_data)
cat(sprintf("  Fitted gamma: %.4f\n", fit_hat$gamma))
cat(sprintf("  Fitted delta: %.4f\n", fit_hat$delta))
cat(sprintf("  Fitted xi:    %.4f\n", fit_hat$xi))
stopifnot(abs(fit_hat$gamma - 0.5) < 0.5)
stopifnot(abs(fit_hat$delta - 1.5) < 0.5)

# 3. Test PDF / LogPDF
cat("Testing johnsonsl_pdf/logpdf...\n")
x <- seq(11, 20, length.out = 5)
dens <- johnsonsl_pdf(x, fit_hat)
log_dens <- johnsonsl_logpdf(x, fit_hat)
stopifnot(all(dens > 0))
stopifnot(all.equal(log(dens), log_dens))

# 4. Test CDF / SF / LogCDF / LogSF
cat("Testing johnsonsl_cdf/sf/logcdf/logsf...\n")
cdf_vals <- johnsonsl_cdf(x, fit_hat)
sf_vals <- johnsonsl_sf(x, fit_hat)
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all.equal(cdf_vals + sf_vals, rep(1, length(x))))
stopifnot(all.equal(log(cdf_vals), johnsonsl_logcdf(x, fit_hat)))
stopifnot(all.equal(log(sf_vals), johnsonsl_logsf(x, fit_hat)))

# 5. Test Quantile / ISF
cat("Testing johnsonsl_quantile/isf...\n")
p <- c(0.25, 0.5, 0.75)
qs <- johnsonsl_quantile(p, fit_hat)
isfs <- johnsonsl_isf(1 - p, fit_hat)
stopifnot(all.equal(qs, isfs))
stopifnot(all.equal(johnsonsl_cdf(qs, fit_hat), p))

# 6. Test Mean / Var / Std
cat("Testing moments/statistics...\n")
m <- johnsonsl_mean(fit_hat)
v <- johnsonsl_var(fit_hat)
s <- johnsonsl_std(fit_hat)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
stopifnot(is.finite(m) && is.finite(v) && is.finite(s))
stopifnot(all.equal(sqrt(v), s))

# 7. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- johnsonsl_median(fit_hat)
int <- johnsonsl_interval(0.95, fit_hat)
m1 <- johnsonsl_moment(1, fit_hat)
stopifnot(all.equal(med, johnsonsl_quantile(0.5, fit_hat)))
stopifnot(all.equal(m1, m, tolerance = 1e-3))

# 8. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- johnsonsl_entropy(fit_hat)
exp_x <- johnsonsl_expect(function(x) x, fit_hat)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))
stopifnot(all.equal(exp_x, m, tolerance = 1e-4))

cat("\n=== Johnson SL tests passed! ===\n")
