# Test suite for Beta-Rayleigh Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Beta-Rayleigh Distribution ===\n\n")

# Parameters: a=2, b=1.5, mu=10, sigma=3
fit_true <- list(a = 2, b = 1.5, mu = 10, sigma = 3, distribution = "betarayleigh_4p")

# 1. Test Rand
cat("Testing betarayleigh_4p_rand...\n")
data <- betarayleigh_4p_rand(1000, fit_true)
cat(sprintf("  Mean of sample: %.4f (expected ~%.4f)\n", 
            mean(data), betarayleigh_4p_mean(fit_true)))

# 2. Test Fit
cat("Testing betarayleigh_4p_fit...\n")
fit_hat <- betarayleigh_4p_fit(data)
cat(sprintf("  Fitted a:     %.4f (expected ~2.0)\n", fit_hat$a))
cat(sprintf("  Fitted b:     %.4f (expected ~1.5)\n", fit_hat$b))
cat(sprintf("  Fitted mu:    %.4f (expected ~10.0)\n", fit_hat$mu))
cat(sprintf("  Fitted sigma: %.4f (expected ~3.0)\n", fit_hat$sigma))

# 3. Test PDF/LogPDF
cat("Testing betarayleigh_4p_pdf/logpdf...\n")
x <- c(12, 15, 20)
pdf_vals <- betarayleigh_4p_pdf(x, fit_true)
logpdf_vals <- betarayleigh_4p_logpdf(x, fit_true)
stopifnot(all(abs(log(pdf_vals) - logpdf_vals) < 1e-10))

# 4. Test CDF/SF
cat("Testing betarayleigh_4p_cdf/sf...\n")
cdf_vals <- betarayleigh_4p_cdf(x, fit_true)
sf_vals <- betarayleigh_4p_sf(x, fit_true)
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-10))

# 5. Test Quantile
cat("Testing betarayleigh_4p_quantile...\n")
p <- c(0.1, 0.5, 0.9)
q_vals <- betarayleigh_4p_quantile(p, fit_true)
stopifnot(all(abs(betarayleigh_4p_cdf(q_vals, fit_true) - p) < 1e-10))

# 6. Test Moments
cat("Testing moments/statistics...\n")
m <- betarayleigh_4p_mean(fit_true)
v <- betarayleigh_4p_var(fit_true)
cat(sprintf("  Mean: %.4f, Var: %.4f\n", m, v))

# 7. Test Truncated Fit
cat("Testing truncated fit...\n")
trunc_data <- data[data >= 11 & data <= 25]
fit_trunc <- betarayleigh_4p_fit_truncated(trunc_data, lower = 11, upper = 25)
cat(sprintf("  Truncated mu: %.4f, sigma: %.4f\n", fit_trunc$mu, fit_trunc$sigma))

cat("\n=== Beta-Rayleigh tests passed! ===\n")
