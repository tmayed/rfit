# POC: Testing the Elbow Method for reasonable regime detection
source("../pkg/segment_data.R")

# Load dataset
data_path <- "poc/inputs/router_traffic_2026_03_30.csv"
if (!file.exists(data_path)) {
  data_path <- "inputs/router_traffic_2026_03_30.csv"
}

df <- read.csv(data_path)
data <- df$traffic
data <- data[is.finite(data) & data > 0]

cat(sprintf("Loaded dataset with %d observations.\n", length(data)))

# Test Elbow Method with different thresholds
thresholds <- c(0.2, 0.1, 0.05, 0.02)

cat("\nElbow Method Tuning Results:\n")
cat(sprintf("%-12s | %-10s | %s\n", "Threshold", "Segments", "Splits"))
cat("--------------------------------------------\n")

for (t in thresholds) {
  fit <- segment_data(data, criterion = "elbow", threshold = t)
  splits_str <- paste(round(fit$splits, 2), collapse = ", ")
  cat(sprintf("%-12.2f | %-10d | %s\n", t, fit$n_segments, splits_str))
}

cat("\nSame test on LOG data:\n")
log_data <- log(data)
for (t in thresholds) {
  fit <- segment_data(log_data, criterion = "elbow", threshold = t)
  # Map splits back to original scale for display
  orig_splits_str <- paste(round(exp(fit$splits), 2), collapse = ", ")
  cat(sprintf("%-12.2f | %-10d | %s\n", t, fit$n_segments, orig_splits_str))
}
