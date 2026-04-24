# POC for segment_data using whole_dataset.csv
source("../pkg/segment_data.R")

# Load dataset
data_path <- "inputs/router_traffic_2026_03_30.csv"
if (!file.exists(data_path)) {
  stop(paste("Data file not found at:", data_path))
}

df <- read.csv(data_path)
data <- df$traffic

cat(sprintf("Loaded dataset with %d observations.\n", length(data)))

# Run segment_data with k=2 to force 2 segments
cat("\nRunning segment_data with k=2 (Forced 2 segments)...\n")
start_time <- Sys.time()
fit_k2 <- segment_data(data, k = 2, n_boot = 5)
end_time <- Sys.time()

cat("\n=== segment_data (k=2) Results ===\n")
cat(sprintf("Execution time: %.2f seconds\n", as.numeric(end_time - start_time, units = "secs")))
cat(sprintf("Segments: %d\n", fit_k2$n_segments))
cat(sprintf("Split 1 (x1): %.4f\n", fit_k2$x1))
cat("Segment counts:", paste(fit_k2$counts, collapse = ", "), "\n")

# Run segment_data with k=3 to force 3 segments (HBT)
cat("\nRunning segment_data with k=3 (Forced 3 segments)...\n")
start_time <- Sys.time()
fit_k3 <- segment_data(data, k = 3, n_boot = 5)
end_time <- Sys.time()

cat("\n=== segment_data (k=3) Results ===\n")
cat(sprintf("Execution time: %.2f seconds\n", as.numeric(end_time - start_time, units = "secs")))
cat(sprintf("Segments: %d\n", fit_k3$n_segments))
cat(sprintf("Split 1 (x1): %.4f\n", fit_k3$x1))
cat(sprintf("Split 2 (x2): %.4f\n", fit_k3$x2))
cat("Segment counts:", paste(fit_k3$counts, collapse = ", "), "\n")

# Run segment_data with k=NULL (Automatic)
cat("\nRunning segment_data with k=NULL (Automatic selection)...\n")
start_time <- Sys.time()
fit_auto <- segment_data(data, k = NULL, n_boot = 5)
end_time <- Sys.time()

cat("\n=== segment_data (Automatic) Results ===\n")
cat(sprintf("Execution time: %.2f seconds\n", as.numeric(end_time - start_time, units = "secs")))
cat(sprintf("Optimal number of segments: %d\n", fit_auto$n_segments))
cat("Segment counts:", paste(fit_auto$counts, collapse = ", "), "\n")
