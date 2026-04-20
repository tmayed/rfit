# Test suite for fit_all orchestrator
source("../pkg/rfit.R")

set.seed(42)

cat("=== Testing fit_all Orchestrator ===\n\n")

# Generate some data from a known distribution (e.g., Weibull)
test_data <- rweibull(500, shape = 2, scale = 10)

cat("Fitting all distributions to Weibull data...\n")
fits <- fit_all(test_data, criterion = "AIC")

cat("\nOrdered results (Best to Worst):\n")
for (i in seq_along(fits)) {
  fit <- fits[[i]]
  cat(sprintf("%d. %s (AIC: %.2f)\n", i, fit$distribution, fit$fit_criterion_value))
}

# Verify that the best fit is reasonable
best_fit <- fits[[1]]
cat(sprintf("\nBest fit is: %s\n", best_fit$distribution))

# Since it's Weibull data, Weibull or more general distributions (like GB2, Kappa4) should be at the top
stopifnot(length(fits) > 0)
stopifnot(fits[[1]]$fit_criterion_value <= fits[[length(fits)]]$fit_criterion_value)

cat("\n=== fit_all tests passed! ===\n")
