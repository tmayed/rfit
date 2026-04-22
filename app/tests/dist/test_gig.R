# Test suite for Generalized Inverse Gaussian (GIG) Distribution
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing GIG Distribution ===\n\n")

# Parameters: lambda=0.5, chi=2.0, psi=2.0
fit_true <- list(lambda = 0.5, chi = 2.0, psi = 2.0, distribution = "gig_3p")

# 1. Test Rand
cat("Testing gig_3p_rand...\n")
data <- gig_3p_rand(500, fit_true)
cat(sprintf("  Median of sample: %.4f (expected ~%.4f)\n", 
            median(data), gig_3p_median(fit_true)))

# 2. Test Fit
cat("Testing gig_3p_fit...\n")
fit_hat <- gig_3p_fit(data)
cat(sprintf("  Fitted lambda: %.4f (expected ~0.5)\n", fit_hat$lambda))
cat(sprintf("  Fitted chi:    %.4f (expected ~2.0)\n", fit_hat$chi))
cat(sprintf("  Fitted psi:    %.4f (expected ~2.0)\n", fit_hat$psi))

# 3. Test PDF/LogPDF
cat("Testing gig_3p_pdf/logpdf...\n")
x <- c(0.5, 1, 2)
pdf_vals <- gig_3p_pdf(x, fit_true)
logpdf_vals <- gig_3p_logpdf(x, fit_true)
stopifnot(all(abs(log(pdf_vals) - logpdf_vals) < 1e-10))

# 4. Test CDF/SF
cat("Testing gig_3p_cdf/sf...\n")
cdf_vals <- gig_3p_cdf(x, fit_true)
sf_vals <- gig_3p_sf(x, fit_true)
stopifnot(all(abs(cdf_vals + sf_vals - 1) < 1e-10))

# 5. Test Quantile
cat("Testing gig_3p_quantile...\n")
p <- c(0.1, 0.5, 0.9)
q_vals <- gig_3p_quantile(p, fit_true)
# Loosen tolerance as numerical integration of PDF is used
stopifnot(all(abs(gig_3p_cdf(q_vals, fit_true) - p) < 1e-4))

# 6. Moments
cat("Testing moments/statistics...\n")
cat(sprintf("  Mean: %.4f, Var: %.4f\n", gig_3p_mean(fit_true), gig_3p_var(fit_true)))

# 7. Test Truncated Fit
cat("Testing truncated fit...\n")
trunc_data <- data[data >= 0.5 & data <= 5]
fit_trunc <- gig_3p_fit_truncated(trunc_data, lower = 0.5, upper = 5)
cat(sprintf("  Truncated chi: %.4f, psi: %.4f\n", fit_trunc$chi, fit_trunc$psi))

cat("\n=== GIG tests passed! ===\n")
