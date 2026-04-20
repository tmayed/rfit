# Test suite for Johnson SB Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Johnson SB Distribution ===\n\n")

# Parameters: gamma=0.5, delta=1.5, xi=0, lambda=10
fit_true <- list(gamma = 0.5, delta = 1.5, xi = 0, lambda = 10)

# 1. Test Random Generation
cat("Testing johnsonsb_rand...\n")
test_data <- johnsonsb_rand(1000, fit_true)
stopifnot(all(test_data >= 0 & test_data <= 10))
cat(sprintf("  Mean of sample: %.4f\n", mean(test_data)))

# 2. Test Fit
cat("Testing johnsonsb_fit...\n")
fit_hat <- johnsonsb_fit(test_data)
cat(sprintf("  Fitted gamma:  %.4f\n", fit_hat$gamma))
cat(sprintf("  Fitted delta:  %.4f\n", fit_hat$delta))
cat(sprintf("  Fitted xi:     %.4f\n", fit_hat$xi))
cat(sprintf("  Fitted lambda: %.4f\n", fit_hat$lambda))
# Tolerance might need to be loose for 4 parameters
stopifnot(abs(fit_hat$gamma - 0.5) < 0.5)
stopifnot(abs(fit_hat$delta - 1.5) < 0.5)

# 3. Test PDF / LogPDF
cat("Testing johnsonsb_pdf/logpdf...\n")
x <- seq(1, 9, length.out = 5)
dens <- johnsonsb_pdf(x, fit_hat)
log_dens <- johnsonsb_logpdf(x, fit_hat)
stopifnot(all(dens > 0))
stopifnot(all.equal(log(dens), log_dens))

# 4. Test CDF / SF / LogCDF / LogSF
cat("Testing johnsonsb_cdf/sf/logcdf/logsf...\n")
cdf_vals <- johnsonsb_cdf(x, fit_hat)
sf_vals <- johnsonsb_sf(x, fit_hat)
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all.equal(cdf_vals + sf_vals, rep(1, length(x))))
stopifnot(all.equal(log(cdf_vals), johnsonsb_logcdf(x, fit_hat)))
stopifnot(all.equal(log(sf_vals), johnsonsb_logsf(x, fit_hat)))

# 5. Test Quantile / ISF
cat("Testing johnsonsb_quantile/isf...\n")
p <- c(0.25, 0.5, 0.75)
qs <- johnsonsb_quantile(p, fit_hat)
isfs <- johnsonsb_isf(1 - p, fit_hat)
stopifnot(all.equal(qs, isfs))
stopifnot(all.equal(johnsonsb_cdf(qs, fit_hat), p))

# 6. Test Mean / Var / Std
cat("Testing moments/statistics...\n")
m <- johnsonsb_mean(fit_hat)
v <- johnsonsb_var(fit_hat)
s <- johnsonsb_std(fit_hat)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))
stopifnot(is.finite(m) && is.finite(v) && is.finite(s))
stopifnot(all.equal(sqrt(v), s))

# 7. Test Median / Interval / Moment
cat("Testing median/interval/moment...\n")
med <- johnsonsb_median(fit_hat)
int <- johnsonsb_interval(0.95, fit_hat)
m1 <- johnsonsb_moment(1, fit_hat)
stopifnot(all.equal(med, johnsonsb_quantile(0.5, fit_hat)))
stopifnot(all.equal(m1, m, tolerance = 1e-5))

# 8. Test Entropy / Expect
cat("Testing entropy/expect...\n")
ent <- johnsonsb_entropy(fit_hat)
exp_x <- johnsonsb_expect(function(x) x, fit_hat)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))
stopifnot(all.equal(exp_x, m, tolerance = 1e-5))

cat("\n=== Johnson SB tests passed! ===\n")
