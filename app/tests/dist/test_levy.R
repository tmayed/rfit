# Test suite for Lévy Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Lévy Distribution ===\n\n")

# Parameters: mu=2, sigma=1
fit_true <- list(mu = 2, sigma = 1, distribution = "levy_2p")

# 1. Test Rand
cat("Testing levy_2p_rand...\n")
data <- levy_2p_rand(500, fit_true)
cat(sprintf("  Median of sample: %.4f (expected ~%.4f)\n", 
            median(data), levy_2p_median(fit_true)))

# 2. Test Fit
cat("Testing levy_2p_fit...\n")
fit_hat <- levy_2p_fit(data)
cat(sprintf("  Fitted mu:    %.4f (expected ~2.0)\n", fit_hat$mu))
cat(sprintf("  Fitted sigma: %.4f (expected ~1.0)\n", fit_hat$sigma))

# 3. Test PDF/LogPDF
cat("Testing levy_2p_pdf/logpdf...\n")
x <- c(3, 4, 10)
pdf_vals <- levy_2p_pdf(x, fit_true)
logpdf_vals <- levy_2p_logpdf(x, fit_true)
stopifnot(all(abs(log(pdf_vals) - logpdf_vals) < 1e-10))

# 4. Test CDF/SF
cat("Testing levy_2p_cdf/sf...\n")
cdf_vals <- levy_2p_cdf(x, fit_true)
sf_vals <- levy_2p_sf(x, fit_true)
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-10))

# 5. Test Quantile
cat("Testing levy_2p_quantile...\n")
p <- c(0.1, 0.5, 0.9)
q_vals <- levy_2p_quantile(p, fit_true)
stopifnot(all(abs(levy_2p_cdf(q_vals, fit_true) - p) < 1e-10))

# 6. Statistics
cat("Testing statistics...\n")
cat(sprintf("  Median: %.4f\n", levy_2p_median(fit_true)))
cat(sprintf("  Mean:   %s\n", as.character(levy_2p_mean(fit_true))))

# 7. Test Truncated Fit
cat("Testing truncated fit...\n")
trunc_data <- data[data >= 3 & data <= 20]
fit_trunc <- levy_2p_fit_truncated(trunc_data, lower = 3, upper = 20)
cat(sprintf("  Truncated mu: %.4f, sigma: %.4f\n", fit_trunc$mu, fit_trunc$sigma))

cat("\n=== Lévy tests passed! ===\n")
