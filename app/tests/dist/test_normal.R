# Test suite for Normal Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Normal Distribution ===\n\n")

# Generate test data
test_data <- rnorm(1000, mean = 5, sd = 2)

# Test log-likelihood
cat("Testing normal_2p_log_likelihood...\n")
ll <- normal_2p_log_likelihood(test_data, 5, 2)
cat(sprintf("  Log-likelihood: %.4f\n", ll))
stopifnot(is.finite(ll))

# Test fit
cat("Testing normal_2p_fit...\n")
fit <- normal_2p_fit(test_data)
cat(sprintf("  mu: %.4f (expected ~5.0)\n", fit$mu))
cat(sprintf("  sigma: %.4f (expected ~2.0)\n", fit$sigma))
stopifnot(abs(fit$mu - 5) < 0.2)
stopifnot(abs(fit$sigma - 2) < 0.2)

# Test fit_truncated
cat("Testing normal_2p_fit_truncated...\n")
fit_t <- normal_2p_fit_truncated(test_data, lower = 0, upper = 10)
cat(sprintf("  mu: %.4f\n", fit_t$mean))
cat(sprintf("  sigma: %.4f\n", fit_t$sd))
stopifnot(is.finite(fit_t$log_likelihood))

# Test PDF/CDF
cat("Testing normal_2p_pdf/cdf...\n")
x <- seq(0, 10, length.out = 5)
dens <- normal_2p_pdf(x, fit)
cdf_vals <- normal_2p_cdf(x, fit)
cat(sprintf("  Densities: %s\n", paste(sprintf("%.4f", dens), collapse = ", ")))
cat(sprintf("  CDF values: %s\n", paste(sprintf("%.4f", cdf_vals), collapse = ", ")))
stopifnot(all(dens >= 0))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))

# Test SF
cat("Testing normal_2p_sf...\n")
sf <- normal_2p_sf(x, fit)
cat(sprintf("  SF: %s\n", paste(sprintf("%.4f", sf), collapse = ", ")))
stopifnot(all(abs(sf - (1 - cdf_vals)) < 1e-10))

# Test Quantile
cat("Testing normal_2p_quantile...\n")
p <- c(0.025, 0.5, 0.975)
qs <- normal_2p_quantile(p, fit)
cat(sprintf("  Quantiles (95%% CI): %.4f, %.4f, %.4f\n", qs[1], qs[2], qs[3]))
stopifnot(all(diff(qs) > 0))

# Test ISF
cat("Testing normal_2p_isf...\n")
isf <- normal_2p_isf(p, fit)
cat(sprintf("  ISF: %s\n", paste(sprintf("%.4f", isf), collapse = ", ")))
stopifnot(all(abs(normal_2p_cdf(isf, fit) - (1 - p)) < 1e-10))

# Test Logs
cat("Testing Log PDF/CDF/SF...\n")
lpdf <- normal_2p_logpdf(x, fit)
lcdf <- normal_2p_logcdf(x, fit)
lsf <- normal_2p_logsf(x, fit)
stopifnot(all(abs(lpdf - log(normal_2p_pdf(x, fit))) < 1e-10))
stopifnot(all(abs(lcdf - log(normal_2p_cdf(x, fit))) < 1e-10))
stopifnot(all(abs(lsf - log(normal_2p_sf(x, fit))) < 1e-10))

# Test normal_2p_rand
cat("Testing normal_2p_rand...\n")
samples <- normal_2p_rand(100, fit)
cat(sprintf("  Generated %d samples, mu: %.4f\n", length(samples), mean(samples)))
stopifnot(length(samples) == 100)

# Test Stats
cat("Testing Stats (Mean, Var, Std, Skew, Kurtosis, Median)...\n")
mean_val <- normal_2p_mean(fit)
var_val <- normal_2p_var(fit)
std_val <- normal_2p_std(fit)
skew_val <- normal_2p_skew(fit)
kurt_val <- normal_2p_kurtosis(fit)
med_val <- normal_2p_median(fit)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f, Skew: %.4f, Kurtosis: %.4f, Median: %.4f\n", 
            mean_val, var_val, std_val, skew_val, kurt_val, med_val))
stopifnot(abs(std_val - sqrt(var_val)) < 1e-10)
stopifnot(skew_val == 0)
stopifnot(kurt_val == 3)
stopifnot(abs(med_val - mean_val) < 1e-10)

# Test Moment
cat("Testing normal_2p_moment...\n")
m1 <- normal_2p_moment(1, fit)
m2 <- normal_2p_moment(2, fit)
cat(sprintf("  1st Moment: %.4f, 2nd Moment: %.4f\n", m1, m2))
# E[X^2] = Var[X] + E[X]^2
stopifnot(abs(m1 - mean_val) < 0.1)
stopifnot(abs(m2 - (var_val + mean_val^2)) < 0.5)

# Test Interval
cat("Testing normal_2p_interval...\n")
int <- normal_2p_interval(0.95, fit)
cat(sprintf("  95%% Interval: [%.4f, %.4f]\n", int[1], int[2]))
stopifnot(abs(normal_2p_cdf(int[2], fit) - normal_2p_cdf(int[1], fit) - 0.95) < 1e-10)

# Test Entropy
cat("Testing normal_2p_entropy...\n")
ent <- normal_2p_entropy(fit)
cat(sprintf("  Entropy: %.4f\n", ent))
stopifnot(is.finite(ent))

# Test Expect
cat("Testing normal_2p_expect...\n")
e_x <- normal_2p_expect(function(x) x, fit)
cat(sprintf("  E[x]: %.4f\n", e_x))
stopifnot(abs(e_x - mean_val) < 1e-5)

cat("\n=== Normal tests passed! ===\n")
