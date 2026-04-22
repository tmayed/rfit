# Test suite for Pareto Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Pareto Distribution ===\n\n")

# Generate test data: shape = 3, scale = 2
# Pareto quantile: scale / (1-u)^(1/shape)
shape_true <- 3.0
scale_true <- 2.0
test_data <- scale_true / (1 - runif(1000))^(1 / shape_true)

# Test log-likelihood
cat("Testing pareto_2p_log_likelihood...\n")
ll <- pareto_2p_log_likelihood(test_data, shape_true, scale_true)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# Test fit
cat("Testing pareto_2p_fit...\n")
fit <- pareto_2p_fit(test_data)
cat(sprintf("  shape: %.4f (expected ~3.0)\n", fit$shape))
cat(sprintf("  scale: %.4f (expected ~2.0)\n", fit$scale))
stopifnot(abs(fit$shape - shape_true) < 0.5)
stopifnot(abs(fit$scale - min(test_data)) < 1e-10)

# Test fit_truncated
cat("Testing pareto_2p_fit_truncated...\n")
lower_bound <- 2.5
upper_bound <- 10.0
trunc_data <- test_data[test_data >= lower_bound & test_data <= upper_bound]
fit_t <- pareto_2p_fit_truncated(trunc_data, lower = lower_bound, upper = upper_bound)
cat(sprintf("  shape: %.4f\n", fit_t$shape))
cat(sprintf("  scale: %.4f\n", fit_t$scale))
stopifnot(is.finite(fit_t$log_likelihood))
stopifnot(fit_t$scale <= lower_bound)

# Test PDF/CDF
cat("Testing pareto_2p_pdf/cdf...\n")
x <- seq(fit$scale, fit$scale * 5, length.out = 10)
dens <- pareto_2p_pdf(x, fit)
cdf_vals <- pareto_2p_cdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", head(dens)), collapse = ", ")))
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", head(cdf_vals)), collapse = ", ")))
stopifnot(all(dens >= 0))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))
stopifnot(all(diff(cdf_vals) >= 0))

# Test SF
cat("Testing pareto_2p_sf...\n")
sf <- pareto_2p_sf(x, fit)
cat(sprintf("  SF: %s\n", paste(sprintf("%.4f", head(sf)), collapse = ", ")))
stopifnot(all(abs(sf - (1 - cdf_vals)) < 1e-10))

# Test Quantile
cat("Testing pareto_2p_quantile...\n")
p <- c(0.025, 0.5, 0.975)
qs <- pareto_2p_quantile(p, fit)
cat(sprintf("  Quantiles: %.4f, %.4f, %.4f\n", qs[1], qs[2], qs[3]))
stopifnot(all(diff(qs) > 0))

# Test ISF
cat("Testing pareto_2p_isf...\n")
isf <- pareto_2p_isf(p, fit)
cat(sprintf("  ISF: %s\n", paste(sprintf("%.4f", isf), collapse = ", ")))
stopifnot(all(abs(pareto_2p_cdf(isf, fit) - (1 - p)) < 1e-10))

# Test Logs
cat("Testing Log PDF/CDF/SF...\n")
lpdf <- pareto_2p_logpdf(x, fit)
lcdf <- pareto_2p_logcdf(x, fit)
lsf <- pareto_2p_logsf(x, fit)
stopifnot(all(abs(lpdf - log(pareto_2p_pdf(x, fit))) < 1e-10))
# CDF can be 0 at scale, so use is.finite for log
stopifnot(all(abs(lcdf[cdf_vals > 0] - log(cdf_vals[cdf_vals > 0])) < 1e-10))
stopifnot(all(abs(lsf - log(pareto_2p_sf(x, fit))) < 1e-10))

# Test pareto_2p_rand
cat("Testing pareto_2p_rand...\n")
samples <- pareto_2p_rand(100, fit)
cat(sprintf("  Generated %d samples, mean: %.4f\n", length(samples), mean(samples)))
stopifnot(length(samples) == 100)

# Test Stats (Mean, Var, Std, Skew, Kurtosis, Median)
cat("Testing Stats (Mean, Var, Std, Skew, Kurtosis, Median)...\n")
mean_val <- pareto_2p_mean(fit)
var_val <- pareto_2p_var(fit)
std_val <- pareto_2p_std(fit)
skew_val <- pareto_2p_skew(fit)
kurt_val <- pareto_2p_kurtosis(fit)
med_val <- pareto_2p_median(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f, Skew: %.4f, Kurtosis: %.4f, Median: %.4f\n", 
            mean_val, var_val, std_val, skew_val, kurt_val, med_val))
stopifnot(abs(std_val - sqrt(var_val)) < 1e-10)
# Theoretical skew for Pareto: (2(1+a)/(a-3)) * sqrt((a-2)/a) for a > 3
# Theoretical kurtosis for Pareto: 3(a-2)(3a^2+a+2)/(a(a-3)(a-4)) for a > 4
# For shape ~ 3, skew might be very large or Inf, kurtosis NA/Inf

# Test Moment
cat("Testing pareto_2p_moment...\n")
m1 <- pareto_2p_moment(1, fit)
m2 <- pareto_2p_moment(2, fit)
cat(sprintf("  1st Moment: %.4f, 2nd Moment: %.4f\n", m1, m2))
stopifnot(abs(m1 - mean_val) < 1e-10)
stopifnot(abs(m2 - (var_val + mean_val^2)) < 1e-10)

# Test Interval
cat("Testing pareto_2p_interval...\n")
int <- pareto_2p_interval(0.95, fit)
cat(sprintf("  95%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(abs(pareto_2p_cdf(int[2], fit) - pareto_2p_cdf(int[1], fit) - 0.95) < 1e-10)

# Test Entropy
cat("Testing pareto_2p_entropy...\n")
ent <- pareto_2p_entropy(fit)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))

# Test Expect
cat("Testing pareto_2p_expect...\n")
e_x <- pareto_2p_expect(function(x) x, fit)
cat(sprintf("  E[x]: %.4f\n", e_x))
stopifnot(abs(e_x - mean_val) < 1e-5)

cat("\n=== Pareto tests passed! ===\n")
