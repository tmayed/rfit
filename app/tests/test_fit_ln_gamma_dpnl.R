# Test suite for fit_ln_gamma_dpnl
source("../pkg/rfit.R")

set.seed(42)

cat("=== Testing fit_ln_gamma_dpnl ===\n\n")

# Generate synthetic data from a known mixture
n <- 1000
w_true <- c(0.5, 0.3, 0.2)
mu_true <- 1.0; sigma_true <- 0.4
shape_true <- 5.0; scale_true <- 0.5
alpha_true <- 3.0; beta_true <- 3.0; nu_true <- 2.0; tau_true <- 0.3

data1 <- rlnorm(n * w_true[1], meanlog = mu_true, sdlog = sigma_true)
data2 <- rgamma(n * w_true[2], shape = shape_true, scale = scale_true)
# dPLN rand: exp(nu + tau*Z + E1 - E2)
z <- rnorm(n * w_true[3])
e1 <- rexp(n * w_true[3], rate = alpha_true)
e2 <- rexp(n * w_true[3], rate = beta_true)
data3 <- exp(nu_true + tau_true * z + e1 - e2)

data <- c(data1, data2, data3)

# Fit
fit <- fit_ln_gamma_dpnl(data, n_starts = 5, maxit = 5000)

cat(sprintf("Weights: %.4f, %.4f, %.4f\n", fit$weights[1], fit$weights[2], fit$weights[3]))
cat(sprintf("Log-likelihood: %.4f\n", fit$log_likelihood))
cat(sprintf("Convergence: %d\n", fit$convergence))

# Assertions
stopifnot(is.list(fit))
stopifnot(fit$distribution == "ln_gamma_dpnl_mixture")
stopifnot(is.finite(fit$log_likelihood))
stopifnot(abs(sum(fit$weights) - 1) < 1e-6)
stopifnot(fit$convergence == 0) # Strict convergence check

# Sane parameters check
stopifnot(fit$components$gamma$shape > 1)
stopifnot(fit$components$dpln$alpha > 1)
stopifnot(fit$components$dpln$beta > 1)
stopifnot(fit$components$lognormal$sigma > 0)

# Accuracy check (broad ranges for reliability across random seeds)
# Lognormal component
stopifnot(abs(fit$components$lognormal$mu - mu_true) < 0.5)
stopifnot(abs(fit$components$lognormal$sigma - sigma_true) < 0.3)

# Weight check
stopifnot(all(fit$weights > 0.01)) # No component should collapse to zero on this data

cat("\n=== fit_ln_gamma_dpnl tests passed! ===\n")
