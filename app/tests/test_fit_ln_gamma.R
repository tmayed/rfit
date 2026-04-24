# Test suite for fit_ln_gamma
source("../pkg/rfit.R")

set.seed(42)

cat("=== Testing fit_ln_gamma ===\n\n")

# Generate synthetic data from a known mixture
n <- 1000
w_true <- c(0.6, 0.4)
mu_true <- 1.2; sigma_true <- 0.3
shape_true <- 4.0; scale_true <- 0.8

data1 <- rlnorm(n * w_true[1], meanlog = mu_true, sdlog = sigma_true)
data2 <- rgamma(n * w_true[2], shape = shape_true, scale = scale_true)

data <- c(data1, data2)

# Fit
fit <- fit_ln_gamma(data, n_starts = 5, maxit = 5000)

cat(sprintf("Weights: %.4f, %.4f\n", fit$weights[1], fit$weights[2]))
cat(sprintf("Log-likelihood: %.4f\n", fit$log_likelihood))
cat(sprintf("Convergence: %d\n", fit$convergence))

# Assertions
stopifnot(is.list(fit))
stopifnot(fit$distribution == "ln_gamma_mixture")
stopifnot(is.finite(fit$log_likelihood))
stopifnot(abs(sum(fit$weights) - 1) < 1e-6)
stopifnot(fit$convergence == 0) # Strict convergence check

# Sane parameters check
stopifnot(fit$components$gamma$shape > 0)
stopifnot(fit$components$lognormal$sigma > 0)

# Accuracy check
stopifnot(abs(fit$components$lognormal$mu - mu_true) < 0.5)
stopifnot(abs(fit$components$lognormal$sigma - sigma_true) < 0.3)

# Weight check
stopifnot(all(fit$weights > 0.01))

cat("\n=== fit_ln_gamma tests passed! ===\n")
