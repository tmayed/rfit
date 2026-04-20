# POC: Fitting 3-component mixture distributions to Danish fire insurance data

# Initialize renv
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../app/rfit.R")
source("../app/plots/cdf_plot.R")

# 1. Load data
data_path <- "inputs/danish.csv"
if (!file.exists(data_path)) {
  stop(sprintf("Data file not found: %s", data_path))
}

danish_data <- read.csv(data_path)$dat
danish_data <- danish_data[!is.na(danish_data) & danish_data > 0]

cat(sprintf("Loaded %d observations from %s\n", length(danish_data), data_path))

# 2. Fit all 3-component mixture combinations
cat("\nFitting all 3-component mixtures using fit_all_3m()...\n")
# Note: This will take several minutes
all_fits <- fit_all_3m(danish_data, criterion = "AIC")

# 3. Display top 10 results ordered by quality
cat("\nTop 10 Mixture Fitting Results (ordered by AIC):\n")
cat(sprintf("%-30s | %-15s | %s\n", "Mixture", "AIC", "Log-Likelihood"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (i in 1:min(10, length(all_fits))) {
  fit <- all_fits[[i]]
  cat(sprintf("%-30s | %-15.2f | %-15.2f\n", 
              fit$distribution_mixture, 
              fit$fit_criterion_value, 
              fit$log_likelihood))
}

# 4. Get the best distribution
best_fit <- all_fits[[1]]
best_name <- best_fit$distribution_mixture

cat(sprintf("\nBest fitting mixture: %s (%s: %.4f)\n", 
            best_name, best_fit$fit_criterion_name, best_fit$fit_criterion_value))

# 5. Create CDF plot for the best fit
output_dir <- "outputs"
output_file <- "danish_best_fit_3m_cdf"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

plot_cdf_comparison(
  sample_data = danish_data,
  fit = best_fit,
  dist_cdf = mixture_cdf,
  output_dir = output_dir,
  output_file = output_file,
  title = sprintf("Empirical vs Fitted CDF (Best 3-Mixture: %s)", best_name),
  x_label = "Loss Amount"
)

cat(sprintf("\nProcessing complete. Best fit plot saved to %s/%s.png\n", output_dir, output_file))
