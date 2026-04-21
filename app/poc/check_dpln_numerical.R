# Numerical check for dPLN
source("pkg/dist/dpln_4p.R")

fit <- list(alpha = 1.4757, beta = 12.2656, nu = 0.6020, tau = 0.0003)

# 1. Check if PDF integrates to 1
pdf_integrand <- function(x) dpln_pdf(x, fit)
total_prob <- integrate(pdf_integrand, lower = 0, upper = Inf)$value
cat(sprintf("Total Probability (integral of PDF): %.6f\n", total_prob))

# 2. Check numerical mean
mean_integrand <- function(x) x * dpln_pdf(x, fit)
numerical_mean <- integrate(mean_integrand, lower = 0, upper = Inf)$value
cat(sprintf("Numerical Mean: %.6f\n", numerical_mean))

# 3. Theoretical mean from dpln_mean
theoretical_mean <- dpln_mean(fit)
cat(sprintf("Theoretical Mean (dpln_mean): %.6f\n", theoretical_mean))
cat(sprintf("Difference: %.6e\n", numerical_mean - theoretical_mean))

# 4. Try with larger tau to ensure it's not a fluke
fit2 <- list(alpha = 3, beta = 2, nu = 1, tau = 0.5)
numerical_mean2 <- integrate(function(x) x * dpln_pdf(x, fit2), lower = 0, upper = Inf)$value
theoretical_mean2 <- dpln_mean(fit2)
cat(sprintf("\nWith tau=0.5:\n"))
cat(sprintf("Numerical Mean: %.6f\n", numerical_mean2))
cat(sprintf("Theoretical Mean: %.6f\n", theoretical_mean2))
cat(sprintf("Difference: %.6e\n", numerical_mean2 - theoretical_mean2))
