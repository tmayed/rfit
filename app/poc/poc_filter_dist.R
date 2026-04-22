# POC: Testing the filter_distributions function

source("../pkg/definitions.R")

cat("--- rfit Distribution Filtering POC ---\n\n")

# Example 0:
cat("Example 0: 2-parameter distributions supported on [0, Inf)\n")
dist_list_1 <- filter_distributions(lb = 0, f_0_strict = TRUE, mean_def = TRUE)
print(names(dist_list_1))
cat("\n")

# Example 1: All 2-parameter distributions with lb = 0
cat("Example 1: 2-parameter distributions supported on [0, Inf)\n")
dist_list_1 <- filter_distributions(np = 2, lb = 0)
print(names(dist_list_1))
cat("\n")

# Example 2: Distributions with well-defined means and lb = 0
cat("Example 2: Distributions supported on [0, Inf) with well-defined means\n")
dist_list_2 <- filter_distributions(lb = 0, mean_def = TRUE)
print(names(dist_list_2))
cat("\n")

# Example 3: Distributions where f(0) = 0 is guaranteed
cat("Example 3: Distributions where f(0) = 0 is strictly guaranteed\n")
dist_list_3 <- filter_distributions(f_0_strict = TRUE)
print(names(dist_list_3))
cat("\n")

# Example 4: Complex filter - 4 parameters, positive support, mean defined
cat("Example 4: 4-parameter distributions, positive support, mean defined\n")
dist_list_4 <- filter_distributions(np = 4, lb = 0, mean_def = TRUE)
if (length(dist_list_4) > 0) {
  print(names(dist_list_4))
} else {
  cat("No distributions match these criteria.\n")
}
cat("\n")

cat("--- POC Complete ---\n")
