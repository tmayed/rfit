# POC for segment_data using router_traffic_2026_03_30.csv
source("../pkg/segment_data.R")

# Load dataset
data_path <- "inputs/router_traffic_2026_03_30.csv"
if (!file.exists(data_path)) {
  data_path <- "poc/inputs/router_traffic_2026_03_30.csv"
}
if (!file.exists(data_path)) {
  stop("Data file 'router_traffic_2026_03_30.csv' not found.")
}

df <- read.csv(data_path)
data <- df$traffic
data <- data[is.finite(data) & data > 0]

cat(sprintf("Loaded dataset with %d observations.\n", length(data)))

# Helper function to print results consistently
print_segment_results <- function(fit, label) {
  cat(sprintf("\n=== segment_data (%s) Results ===\n", label))
  cat(sprintf("Optimal number of segments: %d\n", fit$n_segments))
  cat(sprintf("Number of splits found: %d\n", fit$n_splits))
  
  if (fit$n_splits > 0) {
    for (i in seq_along(fit$splits)) {
      cat(sprintf("  Split %d (x%d): %.4f\n", i, i, fit$splits[i]))
    }
  }
  
  cat("Segment counts:", paste(fit$counts, collapse = ", "), "\n")
  cat(sprintf("Objective: %.4f\n", fit$objective))
}

# Run segment_data with k=2 to force 2 segments
cat("\nRunning segment_data with k=2 (Forced 2 segments)...\n")
fit_k2 <- segment_data(data, k = 2, n_boot = 5)
print_segment_results(fit_k2, "k=2")

# Run segment_data with k=3 to force 3 segments (HBT)
cat("\nRunning segment_data with k=3 (Forced 3 segments)...\n")
fit_k3 <- segment_data(data, k = 3, n_boot = 5)
print_segment_results(fit_k3, "k=3")

# Run segment_data with k=4 to force 4 segments
cat("\nRunning segment_data with k=4 (Forced 4 segments)...\n")
fit_k4 <- segment_data(data, k = 4, n_boot = 5)
print_segment_results(fit_k4, "k=4")

# Run segment_data with k=NULL (Automatic)
cat("\nRunning segment_data with k=NULL (Automatic selection)...\n")
fit_auto <- segment_data(data, k = NULL, n_boot = 5)
print_segment_results(fit_auto, "Automatic")
