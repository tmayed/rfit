# Test suite for fit_all_2m orchestrator
source("../pkg/rfit.R")

set.seed(42)

cat("=== Testing fit_all_2m Orchestrator ===\n\n")

# Generate some data from a 2-component mixture (Lognormal + Weibull)
w1 <- 0.4
data1 <- rlnorm(200, meanlog = 1, sdlog = 0.5)
data2 <- rweibull(300, shape = 2, scale = 10)
test_data <- c(data1, data2)

cat("Fitting all 2-component mixture combinations...\n")
# To keep test time reasonable, we might want to subset distributions
# but the request was for all combinations.
# Let's run it and see.

# For the test, I'll use a small subset to ensure it finishes quickly if it's too slow
# but fit_all_2m is already defined to use the main ones.

fits <- fit_all_2m(test_data, criterion = "AIC")

cat("\nTop 5 Results (Best to Worst):\n")
for (i in 1:min(5, length(fits))) {
  fit <- fits[[i]]
  cat(sprintf("%d. %s (AIC: %.2f)\n", i, fit$distribution_mixture, fit$fit_criterion_value))
}

stopifnot(length(fits) > 0)
stopifnot(fits[[1]]$fit_criterion_value <= fits[[length(fits)]]$fit_criterion_value)

cat("\n=== fit_all_2m tests passed! ===\n")
