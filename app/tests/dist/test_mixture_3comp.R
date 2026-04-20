# Test for 3-component mixture
source("../../pkg/rfit.R")

set.seed(42)

cat("=== Testing 3-Component Mixture ===\n\n")

# Generate 3-component data
n <- 1500
# 1/3 Lognormal, 1/3 Pareto, 1/3 Kappa4
d1 <- lognormal_rand(n/3, list(mu = 1, sigma = 0.2))
d2 <- pareto_rand(n/3, list(shape = 4, scale = 10))
d3 <- kappa4_rand(n/3, list(xi = 20, alpha = 2, k = 0.1, h = 0.5))
data <- c(d1, d2, d3)

dist_names <- c("lognormal", "pareto", "kappa4")
cat(sprintf("Fitting mixture of: %s...\n", paste(dist_names, collapse = ", ")))

fit <- mixture_fit(data, dist_names)

cat("\nFitted weights (expected ~0.33 each):\n")
print(round(fit$weights, 3))

cat("\nLog-likelihood:\n")
print(fit$log_likelihood)

# Basic validation
stopifnot(length(fit$weights) == 3)
stopifnot(abs(sum(fit$weights) - 1) < 1e-7)

cat("\n=== 3-Component Mixture test passed! ===\n")
