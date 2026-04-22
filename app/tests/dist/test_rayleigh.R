# Test suite for Rayleigh Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Rayleigh Distribution ===\n\n")

# Parameters: mu=5, sigma=2
fit_true <- list(mu = 5, sigma = 2, distribution = "rayleigh_2p")

# 1. Test Rand
cat("Testing rayleigh_2p_rand...\n")
data <- rayleigh_2p_rand(1000, fit_true)
cat(sprintf("  Mean of sample: %.4f (expected ~%.4f)\n", 
            mean(data), rayleigh_2p_mean(fit_true)))

# 2. Test Fit
cat("Testing rayleigh_2p_fit...\n")
fit_hat <- rayleigh_2p_fit(data)
cat(sprintf("  Fitted mu:    %.4f (expected ~5.0)\n", fit_hat$mu))
cat(sprintf("  Fitted sigma: %.4f (expected ~2.0)\n", fit_hat$sigma))

# 3. Test PDF/LogPDF
cat("Testing rayleigh_2p_pdf/logpdf...\n")
x <- c(6, 7, 8)
pdf_vals <- rayleigh_2p_pdf(x, fit_true)
logpdf_vals <- rayleigh_2p_logpdf(x, fit_true)
stopifnot(all(abs(log(pdf_vals) - logpdf_vals) < 1e-10))

# 4. Test CDF/SF/LogCDF/LogSF
cat("Testing rayleigh_2p_cdf/sf/logcdf/logsf...\n")
cdf_vals <- rayleigh_2p_cdf(x, fit_true)
sf_vals <- rayleigh_2p_sf(x, fit_true)
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-10))

# 5. Test Quantile/ISF
cat("Testing rayleigh_2p_quantile/isf...\n")
p <- c(0.1, 0.5, 0.9)
q_vals <- rayleigh_2p_quantile(p, fit_true)
stopifnot(all(abs(rayleigh_2p_cdf(q_vals, fit_true) - p) < 1e-10))

# 6. Test Moments
cat("Testing moments/statistics...\n")
m <- rayleigh_2p_mean(fit_true)
v <- rayleigh_2p_var(fit_true)
s <- rayleigh_2p_std(fit_true)
cat(sprintf("  Mean: %.4f, Var: %.4f, Std: %.4f\n", m, v, s))

# 7. Test Truncated Fit
cat("Testing truncated fit...\n")
trunc_data <- data[data >= 6 & data <= 10]
fit_trunc <- rayleigh_2p_fit_truncated(trunc_data, lower = 6, upper = 10)
cat(sprintf("  Truncated mu: %.4f, sigma: %.4f\n", fit_trunc$mu, fit_trunc$sigma))

cat("\n=== Rayleigh tests passed! ===\n")
