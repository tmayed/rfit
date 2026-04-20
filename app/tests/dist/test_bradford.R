# Test suite for Bradford Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Bradford Distribution ===\n\n")

# Parameters: shape=5, lower=0, upper=10
fit_true <- list(shape = 5, lower = 0, upper = 10)

# 1. Test Random Generation
cat("Testing bradford_rand...\n")
test_data <- bradford_rand(1000, fit_true)
stopifnot(all(test_data >= 0 & test_data <= 10))
cat(sprintf("  Mean of sample: %.4f\n", mean(test_data)))

# 2. Test Fit
cat("Testing bradford_fit...\n")
fit_hat <- bradford_fit(test_data)
cat(sprintf("  Fitted shape: %.4f (expected ~5.0)\n", fit_hat$shape))
cat(sprintf("  Fitted lower: %.4f (expected ~0.0)\n", fit_hat$lower))
cat(sprintf("  Fitted upper: %.4f (expected ~10.0)\n", fit_hat$upper))
stopifnot(abs(fit_hat$shape - 5) < 1.0) # Bradford shape fitting can be noisy

# 3. Test PDF / LogPDF
cat("Testing bradford_pdf/logpdf...\n")
x <- seq(1, 9, length.out = 5)
dens <- bradford_pdf(x, fit_hat)
log_dens <- bradford_logpdf(x, fit_hat)
stopifnot(all(dens > 0))
stopifnot(all.equal(log(dens), log_dens))

# 4. Test CDF / SF / LogCDF / LogSF
cat("Testing bradford_cdf/sf/logcdf/logsf...\n")
cdf_vals <- bradford_cdf(x, fit_hat)
sf_vals <- bradford_sf(x, fit_hat)
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all.equal(cdf_vals + sf_vals, rep(1, length(x))))
stopifnot(all.equal(log(cdf_vals), bradford_logcdf(x, fit_hat)))
stopifnot(all.equal(log(sf_vals), bradford_logsf(x, fit_hat)))

# 5. Test Quantile / ISF
cat("Testing bradford_quantile/isf...\n")
p <- c(0.25, 0.5, 0.75)
qs <- bradford_quantile(p, fit_hat)
isfs <- bradford_isf(1 - p, fit_hat)
stopifnot(all.equal(qs, isfs))
stopifnot(all.equal(bradford_cdf(qs, fit_hat), p))

# 6. Test Mean / Var / Std / Skew / Kurtosis
cat("Testing moments/statistics...\n")
m <- bradford_mean(fit_hat)
v <- bradford_var(fit_hat)
s <- bradford_std(fit_hat)
sk <- bradford_skew(fit_hat)
kt <- bradford_kurtosis(fit_hat)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
stopifnot(is.finite(m) && is.finite(v) && is.finite(s))
stopifnot(all.equal(sqrt(v), s))

# 7. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- bradford_median(fit_hat)
int <- bradford_interval(0.95, fit_hat)
m1 <- bradford_moment(1, fit_hat)
stopifnot(all.equal(med, bradford_quantile(0.5, fit_hat)))
stopifnot(all.equal(m1, m))

# 8. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- bradford_entropy(fit_hat)
exp_x <- bradford_expect(function(x) x, fit_hat)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))
stopifnot(all.equal(exp_x, m))

cat("\n=== Bradford tests passed! ===\n")
