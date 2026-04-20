# POC: Fitting distributions to Danish fire insurance data

# Initialize renv
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../pkg/rfit.R")
source("../pkg/plots/cdf_plot.R")

# Prevent Rplots.pdf from being created
if (!interactive()) {
  pdf(NULL)
}

# 1. Load data
data_path <- "inputs/danish.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found: %s", data_path))
}

danish_data <- read.csv(data_path)$dat
danish_data <- danish_data[!is.na(danish_data) & danish_data > 0]

cat(sprintf("Loaded %d observations from %s\n", length(danish_data), data_path))

# 2. Fit all available distributions
cat("\nFitting all distributions using fit_all()...\n")
all_fits <- fit_all(danish_data, criterion = "AIC")

# 3. Display results ordered by quality
cat("\nDistribution Fitting Results (ordered by AIC):\n")
cat(sprintf("%-15s | %-15s | %s\n", "Distribution", "AIC", "Log-Likelihood"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (i in seq_along(all_fits)) {
  fit <- all_fits[[i]]
  cat(sprintf("%-15s | %-15.2f | %-15.2f\n", 
              fit$distribution, 
              fit$fit_criterion_value, 
              fit$log_likelihood))
}

# 4. Get the best distribution
best_fit <- all_fits[[1]]
best_name <- best_fit$distribution
best_cdf_func <- get(paste0(best_name, "_cdf"))

cat(sprintf("\nBest fitting distribution: %s (%s: %.4f)\n", 
            best_name, best_fit$fit_criterion_name, best_fit$fit_criterion_value))

# 5. Create CDF plot for the best fit
output_dir <- "outputs"
output_file <- "danish_best_fit_cdf"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

plot_cdf_comparison(
  sample_data = danish_data,
  fit = best_fit,
  dist_cdf = best_cdf_func,
  output_dir = output_dir,
  output_file = output_file,
  title = sprintf("Empirical vs Fitted CDF (Best: %s)", best_name),
  x_label = "Loss Amount"
)

cat(sprintf("\nProcessing complete. Best fit plot saved to %s/%s.png\n", output_dir, output_file))
