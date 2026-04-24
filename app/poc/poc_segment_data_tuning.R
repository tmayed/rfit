# POC: Tuning the BIC penalty for reasonable regime detection
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

# Test different penalties
penalties <- c(1, 2, 5, 10, 20)

cat("\nBIC Penalty Tuning Results:\n")
cat(sprintf("%-10s | %-10s | %s\n", "Penalty", "Segments", "Splits"))
cat("------------------------------------------\n")

for (p in penalties) {
  fit <- segment_data(data, criterion = "BIC", penalty = p)
  splits_str <- paste(round(fit$splits, 2), collapse = ", ")
  cat(sprintf("%-10.1f | %-10d | %s\n", p, fit$n_segments, splits_str))
}

cat("\nSame test on LOG data:\n")
log_data <- log(data)
for (p in penalties) {
  fit <- segment_data(log_data, criterion = "BIC", penalty = p)
  # Map splits back to original scale for display
  orig_splits_str <- paste(round(exp(fit$splits), 2), collapse = ", ")
  cat(sprintf("%-10.1f | %-10d | %s\n", p, fit$n_segments, orig_splits_str))
}
