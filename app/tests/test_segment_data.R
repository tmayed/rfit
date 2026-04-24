# Test suite for segment_data
source("../pkg/segment_data.R")

set.seed(42)

cat("=== Testing segment_data ===\n\n")

data_three <- c(
  stats::runif(45, min = -3, max = 1),
  stats::runif(50, min = 5, max = 9),
  stats::runif(40, min = 15, max = 19)
)

# Test 1: Automatic selection (should find 3)
fit_auto <- segment_data(data_three, n_boot = 10)

cat("Three-regime data (Auto)\n")
cat(sprintf("x1: %.4f\n", fit_auto$x1))
cat(sprintf("x2: %.4f\n", fit_auto$x2))
cat(sprintf("Segments: %d, Splits: %d\n", fit_auto$n_segments, fit_auto$n_splits))
cat(sprintf("Counts: %s\n", paste(fit_auto$counts, collapse = ", ")))

stopifnot(fit_auto$n_segments == 3L)
stopifnot(fit_auto$n_splits == 2L)

# Test 2: Forced k=3 (previously HBT)
fit_k3 <- segment_data(data_three, k = 3, n_boot = 5)
cat("\nThree-regime data (Forced k=3)\n")
cat(sprintf("Segments: %d, Splits: %d\n", fit_k3$n_segments, fit_k3$n_splits))
stopifnot(fit_k3$n_segments == 3L)
stopifnot(fit_k3$method == "Forced segmentation into k=3 segments")

# Test 3: Forced k=2
fit_k2 <- segment_data(data_three, k = 2, n_boot = 5)
cat("\nThree-regime data (Forced k=2)\n")
cat(sprintf("Segments: %d, Splits: %d\n", fit_k2$n_segments, fit_k2$n_splits))
stopifnot(fit_k2$n_segments == 2L)
stopifnot(fit_k2$n_splits == 1L)

# Test 4: Single-regime data (Auto)
data_one <- stats::rnorm(120, mean = 0, sd = 1)
fit_one <- segment_data(data_one, n_boot = 10)

cat("\nSingle-regime data\n")
cat(sprintf("Segments: %d, Splits: %d\n", fit_one$n_segments, fit_one$n_splits))

stopifnot(fit_one$n_segments == 1L)
stopifnot(fit_one$n_splits == 0L)

# Test 5: Single-regime data forced k=4
fit_one_k4 <- segment_data(data_one, k = 4, n_boot = 5)
cat("\nSingle-regime data (Forced k=4)\n")
cat(sprintf("Segments: %d, Splits: %d\n", fit_one_k4$n_segments, fit_one_k4$n_splits))
stopifnot(fit_one_k4$n_segments == 4L)

cat("\n=== segment_data tests passed! ===\n")
