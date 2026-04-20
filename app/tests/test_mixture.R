# Test suite for Mixture Distribution
source("../pkg/rfit.R")

set.seed(42)

cat("=== Testing Mixture Distribution ===\n\n")

# Generate mixture data: 50% lognormal, 50% pareto
n <- 1000
cat(sprintf("Generating %d mixture samples (50%% lognormal, 50%% pareto)...\n", n))

data_ln <- lognormal_rand(n/2, list(mu = 1, sigma = 0.5))
data_pa <- pareto_rand(n/2, list(shape = 3, scale = 1))
data <- c(data_ln, data_pa)

# Fit mixture
dist_names <- c("lognormal", "pareto")
cat(sprintf("Fitting mixture of: %s...\n", paste(dist_names, collapse = ", ")))

fit <- mixture_fit(data, dist_names)

cat("Fitted weights:\n")
print(fit$weights)

cat("\nFitted components:\n")
print(fit$components)

# Test PDF
x <- seq(min(data), max(data), length.out = 10)
dens <- mixture_pdf(x, fit)
cat(sprintf("\nPDF values: %s\n", paste(sprintf("%.4f", head(dens)), collapse = ", ")))
stopifnot(all(dens >= 0))

# Test CDF
cdf_vals <- mixture_cdf(x, fit)
cat(sprintf("CDF values: %s\n", paste(sprintf("%.4f", head(cdf_vals)), collapse = ", ")))
stopifnot(all(cdf_vals >= 0 & cdf_vals <= 1))

# Test Random Generation
samples <- mixture_rand(100, fit)
cat(sprintf("Generated %d samples, mean: %.4f\n", length(samples), mean(samples)))
stopifnot(length(samples) == 100)

cat("\n=== Mixture tests passed! ===\n")
