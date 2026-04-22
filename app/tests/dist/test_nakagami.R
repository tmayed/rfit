# Test suite for Nakagami Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Nakagami Distribution ===\n\n")

# Parameters: m=1.5, omega=10
fit_true <- list(m = 1.5, omega = 10, distribution = "nakagami_2p")

# 1. Test Rand
cat("Testing nakagami_2p_rand...\n")
data <- nakagami_2p_rand(1000, fit_true)
cat(sprintf("  Mean of X^2 in sample: %.4f (expected ~%.4f)\n", 
            mean(data^2), fit_true$omega))

# 2. Test Fit
cat("Testing nakagami_2p_fit...\n")
fit_hat <- nakagami_2p_fit(data)
cat(sprintf("  Fitted m:     %.4f (expected ~1.5)\n", fit_hat$m))
cat(sprintf("  Fitted omega: %.4f (expected ~10.0)\n", fit_hat$omega))

# 3. Test PDF/LogPDF
cat("Testing nakagami_2p_pdf/logpdf...\n")
x <- c(1, 3, 5)
pdf_vals <- nakagami_2p_pdf(x, fit_true)
logpdf_vals <- nakagami_2p_logpdf(x, fit_true)
stopifnot(all(abs(log(pdf_vals) - logpdf_vals) < 1e-10))

# 4. Test CDF/SF
cat("Testing nakagami_2p_cdf/sf...\n")
cdf_vals <- nakagami_2p_cdf(x, fit_true)
sf_vals <- nakagami_2p_sf(x, fit_true)
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-10))

# 5. Test Quantile
cat("Testing nakagami_2p_quantile...\n")
p <- c(0.1, 0.5, 0.9)
q_vals <- nakagami_2p_quantile(p, fit_true)
stopifnot(all(abs(nakagami_2p_cdf(q_vals, fit_true) - p) < 1e-10))

# 6. Moments
cat("Testing moments/statistics...\n")
m_theoretical <- nakagami_2p_mean(fit_true)
cat(sprintf("  Mean: %.4f, Std: %.4f\n", m_theoretical, nakagami_2p_std(fit_true)))

# 7. Test Truncated Fit
cat("Testing truncated fit...\n")
trunc_data <- data[data >= 2 & data <= 6]
fit_trunc <- nakagami_2p_fit_truncated(trunc_data, lower = 2, upper = 6)
cat(sprintf("  Truncated m: %.4f, omega: %.4f\n", fit_trunc$m, fit_trunc$omega))

cat("\n=== Nakagami tests passed! ===\n")
