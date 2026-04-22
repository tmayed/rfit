# Test suite for GB2 Distribution
source("../../pkg/rfit.R")

# Set seed for reproducibility
set.seed(42)

cat("=== Testing GB2 Distribution ===\n\n")

# Generate test data from GB2
# Parameters: a=2, b=10, p=1.5, q=2.5
true_params <- list(a = 2, b = 10, p = 1.5, q = 2.5)
test_data <- gb2_4p_rand(1000, true_params)

# Test gb2_4p_log_likelihood
cat("Testing gb2_4p_log_likelihood...\n")
ll <- gb2_4p_log_likelihood(test_data, 2, 10, 1.5, 2.5)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# Test gb2_4p_fit
cat("Testing gb2_4p_fit...\n")
fit <- gb2_4p_fit(test_data)
cat(sprintf("  a: %.4f (expected ~2.0)\n", fit$a))
cat(sprintf("  b: %.4f (expected ~10.0)\n", fit$b))
cat(sprintf("  p: %.4f (expected ~1.5)\n", fit$p))
cat(sprintf("  q: %.4f (expected ~2.5)\n", fit$q))
cat(sprintf("  log_likelihood: %.4f\n", fit$log_likelihood))
cat(sprintf("  convergence: %d\n", fit$convergence))
# Tolerances are wider for 4-param distributions
stopifnot(abs(fit$a - 2) < 0.5)
stopifnot(abs(fit$p - 1.5) < 0.5)

# Test gb2_4p_fit_truncated
cat("Testing gb2_4p_fit_truncated...\n")
trunc_data <- test_data[test_data > 2 & test_data < 50]
fit_trunc <- gb2_4p_fit_truncated(trunc_data, lower = 2, upper = 50)
cat(sprintf("  a: %.4f\n", fit_trunc$a))
cat(sprintf("  b: %.4f\n", fit_trunc$b))
cat(sprintf("  convergence: %d\n", fit_trunc$convergence))
stopifnot(is.finite(fit_trunc$log_likelihood))

# Test PDF/CDF
cat("Testing gb2_4p_pdf/cdf...\n")
x <- seq(1, 50, length.out = 6)
dens <- gb2_4p_pdf(x, fit)
cdf_vals <- gb2_4p_cdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", dens), collapse = ", ")))
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", cdf_vals), collapse = ", ")))
stopifnot(all(dens >= 0))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))

# Test SF
cat("Testing gb2_4p_sf...\n")
sf <- gb2_4p_sf(x, fit)
cat(sprintf("  SF: %s\n", paste(sprintf("%.4f", sf), collapse = ", ")))
stopifnot(all(abs(sf - (1 - cdf_vals)) < 1e-10))

# Test Quantile
cat("Testing gb2_4p_quantile...\n")
p <- seq(0.1, 0.9, length.out = 10)
quantiles <- gb2_4p_quantile(p, fit)
cat(sprintf("  Quantiles: %s\n", paste(sprintf("%.4f", quantiles), collapse = ", ")))
stopifnot(diff(quantiles) >= 0)

# Test ISF
cat("Testing gb2_4p_isf...\n")
isf <- gb2_4p_isf(p, fit)
cat(sprintf("  ISF: %s\n", paste(sprintf("%.4f", isf), collapse = ", ")))
stopifnot(all(abs(gb2_4p_cdf(isf, fit) - (1 - p)) < 1e-10))

# Test Logs
cat("Testing Log PDF/CDF/SF...\n")
lpdf <- gb2_4p_logpdf(x, fit)
lcdf <- gb2_4p_logcdf(x, fit)
lsf <- gb2_4p_logsf(x, fit)
stopifnot(all(abs(lpdf - log(gb2_4p_pdf(x, fit))) < 1e-10))
stopifnot(all(abs(lcdf - log(gb2_4p_cdf(x, fit))) < 1e-10))
stopifnot(all(abs(lsf - log(gb2_4p_sf(x, fit))) < 1e-10))

# Test Stats
cat("Testing Stats (Mean, Var, Std, Skew, Kurtosis, Median)...\n")
mean_val <- gb2_4p_mean(fit)
var_val <- gb2_4p_var(fit)
std_val <- gb2_4p_std(fit)
skew_val <- gb2_4p_skew(fit)
kurt_val <- gb2_4p_kurtosis(fit)
med_val <- gb2_4p_median(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f, Skew: %.4f, Kurtosis: %.4f, Median: %.4f\n", 
            mean_val, var_val, std_val, skew_val, kurt_val, med_val))
stopifnot(is.finite(mean_val) && is.finite(var_val))
stopifnot(abs(std_val - sqrt(var_val)) < 1e-10)
stopifnot(abs(med_val - gb2_4p_quantile(0.5, fit)) < 1e-10)

# Test Moment
cat("Testing gb2_4p_moment...\n")
m1 <- gb2_4p_moment(1, fit)
m2 <- gb2_4p_moment(2, fit)
cat(sprintf("  1st Moment: %.4f, 2nd Moment: %.4f\n", m1, m2))
stopifnot(abs(m1 - mean_val) < 1e-10)

# Test Interval
cat("Testing gb2_4p_interval...\n")
int <- gb2_4p_interval(0.90, fit)
cat(sprintf("  90%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(abs(gb2_4p_cdf(int[2], fit) - gb2_4p_cdf(int[1], fit) - 0.90) < 1e-10)

# Test Entropy/Expect
cat("Testing gb2_4p_entropy/expect...\n")
ent <- gb2_4p_entropy(fit)
e_x <- gb2_4p_expect(function(x) x, fit)
cat(sprintf("  Entropy: %.4f, E[x]: %.4f\n", ent, e_x))
stopifnot(is.finite(ent))
stopifnot(abs(e_x - mean_val) < 1e-4)

# Test AIC/BIC compatibility
cat("Testing AIC/BIC compatibility...\n")
log_lik_val <- logLik(fit)
aic_val <- AIC(log_lik_val)
bic_val <- BIC(log_lik_val)
cat(sprintf("  AIC: %.4f\n", aic_val))
cat(sprintf("  BIC: %.4f\n", bic_val))
stopifnot(is.finite(aic_val) && is.finite(bic_val))

cat("\n=== GB2 tests passed! ===\n")
