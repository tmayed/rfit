# Test dPLN consistency
source("pkg/dist/dpln_4p.R")

set.seed(42)
true_alpha <- 3.0
true_beta <- 2.0
true_nu <- 1.0
true_tau <- 0.5

true_fit <- list(alpha = true_alpha, beta = true_beta, nu = true_nu, tau = true_tau)
theoretical_mean <- dpln_mean(true_fit)

cat(sprintf("True parameters: alpha=%.2f, beta=%.2f, nu=%.2f, tau=%.2f\n", 
            true_alpha, true_beta, true_nu, true_tau))
cat(sprintf("Theoretical Mean: %.4f\n", theoretical_mean))

# Generate large sample
n <- 100000
sample_data <- dpln_rand(n, true_fit)
empirical_mean <- mean(sample_data)

cat(sprintf("Sample Mean (n=%d): %.4f\n", n, empirical_mean))
cat(sprintf("Difference: %.4f\n", empirical_mean - theoretical_mean))

# Fit the data
cat("\nFitting dPLN to sample...\n")
fitted_res <- dpln_fit(sample_data)

cat("\nFitted parameters:\n")
cat(sprintf("  alpha: %.4f\n", fitted_res$alpha))
cat(sprintf("  beta:  %.4f\n", fitted_res$beta))
cat(sprintf("  nu:    %.4f\n", fitted_res$nu))
cat(sprintf("  tau:   %.4f\n", fitted_res$tau))

fitted_theoretical_mean <- dpln_mean(fitted_res)
cat(sprintf("\nFitted Theoretical Mean: %.4f\n", fitted_theoretical_mean))
cat(sprintf("Sample Mean:            %.4f\n", empirical_mean))
cat(sprintf("Difference:             %.4f\n", fitted_theoretical_mean - empirical_mean))
