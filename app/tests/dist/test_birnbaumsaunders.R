# Test suite for Birnbaum-Saunders Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing Birnbaum-Saunders Distribution ===\n\n")

# Parameters: alpha=0.5, beta=10, mu=5
fit_true <- list(alpha = 0.5, beta = 10, mu = 5, distribution = "birnbaumsaunders_3p")

# 1. Test Rand
cat("Testing birnbaumsaunders_3p_rand...\n")
data <- birnbaumsaunders_3p_rand(1000, fit_true)
cat(sprintf("  Median of sample: %.4f (expected ~%.4f)\n", 
            median(data), birnbaumsaunders_3p_median(fit_true)))

# 2. Test Fit
cat("Testing birnbaumsaunders_3p_fit...\n")
fit_hat <- birnbaumsaunders_3p_fit(data)
cat(sprintf("  Fitted alpha: %.4f (expected ~0.5)\n", fit_hat$alpha))
cat(sprintf("  Fitted beta:  %.4f (expected ~10.0)\n", fit_hat$beta))
cat(sprintf("  Fitted mu:    %.4f (expected ~5.0)\n", fit_hat$mu))

# 3. Test PDF/LogPDF
cat("Testing birnbaumsaunders_3p_pdf/logpdf...\n")
x <- c(12, 15, 20)
pdf_vals <- birnbaumsaunders_3p_pdf(x, fit_true)
logpdf_vals <- birnbaumsaunders_3p_logpdf(x, fit_true)
stopifnot(all(abs(log(pdf_vals) - logpdf_vals) < 1e-10))

# 4. Test CDF/SF
cat("Testing birnbaumsaunders_3p_cdf/sf...\n")
cdf_vals <- birnbaumsaunders_3p_cdf(x, fit_true)
sf_vals <- birnbaumsaunders_3p_sf(x, fit_true)
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-10))

# 5. Test Quantile
cat("Testing birnbaumsaunders_3p_quantile...\n")
p <- c(0.1, 0.5, 0.9)
q_vals <- birnbaumsaunders_3p_quantile(p, fit_true)
stopifnot(all(abs(birnbaumsaunders_3p_cdf(q_vals, fit_true) - p) < 1e-10))

# 6. Moments
cat("Testing moments/statistics...\n")
m <- birnbaumsaunders_3p_mean(fit_true)
v <- birnbaumsaunders_3p_var(fit_true)
cat(sprintf("  Mean: %.4f, Var: %.4f\n", m, v))

# 7. Test Truncated Fit
cat("Testing truncated fit...\n")
trunc_data <- data[data >= 10 & data <= 30]
fit_trunc <- birnbaumsaunders_3p_fit_truncated(trunc_data, lower = 10, upper = 30)
cat(sprintf("  Truncated mu: %.4f, beta: %.4f\n", fit_trunc$mu, fit_trunc$beta))

cat("\n=== Birnbaum-Saunders tests passed! ===\n")
