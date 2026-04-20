# Test suite for fit_all_3m orchestrator
source("../pkg/rfit.R")

set.seed(42)

cat("=== Testing fit_all_3m Orchestrator ===\n\n")

# Generate some data from a 3-component mixture (Lognormal + Weibull + Normal)
data1 <- rlnorm(100, meanlog = 0.5, sdlog = 0.3)
data2 <- rweibull(150, shape = 2, scale = 5)
data3 <- rnorm(100, mean = 15, sd = 2)
data3 <- data3[data3 > 0] # Ensure positive
test_data <- c(data1, data2, data3)

cat("Fitting all 3-component mixture combinations (subset for speed if needed, but let's try)...\n")
# Note: This might be slow if we do all 84 combinations in a test.
# For testing purposes, we might want to check if the function exists and runs on a few.
# But let's run it and see.

fits <- fit_all_3m(test_data, criterion = "AIC")

cat("\nTop 5 Results (Best to Worst):\n")
for (i in 1:min(5, length(fits))) {
  fit <- fits[[i]]
  cat(sprintf("%d. %s (AIC: %.2f)\n", i, fit$distribution_mixture, fit$fit_criterion_value))
}

stopifnot(length(fits) > 0)
stopifnot(fits[[1]]$fit_criterion_value <= fits[[length(fits)]]$fit_criterion_value)

cat("\n=== fit_all_3m tests passed! ===\n")
