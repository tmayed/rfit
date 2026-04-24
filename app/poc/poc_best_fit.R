# POC: Finding the best single distribution fit for input data

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
input_file <- "inputs/router_traffic_2026_03_30.csv"
if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s", input_file))
}

data_raw <- read.csv(input_file)$traffic
data_clean <- data_raw[!is.na(data_raw) & data_raw > 0]

# 2. Take a subset of around 2000 data points
set.seed(42)
subset_size <- min(2000, length(data_clean))
sample_data <- sample(data_clean, subset_size)

cat(sprintf("Loaded %d observations from %s, using subset of %d\n", 
            length(data_clean), input_file, subset_size))

# 3. Fit all single distributions
cat("\nFitting all single distributions using fit_all()...\n")
# Note: Single distributions may have lower BIC than mixtures if the mixture
# does not improve the log-likelihood enough to justify the extra parameters.
all_fits <- fit_all(sample_data, criterion = "BIC")

# 4. Display top results
cat("\nTop Distribution Fitting Results (ordered by BIC):\n")
cat(sprintf("%-25s | %-12s | %-12s | %s\n", "Distribution", "AIC", "BIC", "Log-Lik"))
cat(paste(rep("-", 75), collapse = ""), "\n")

dist_names <- names(all_fits)
for (i in 1:min(25, length(all_fits))) {
  fit <- all_fits[[i]]
  dist_name <- dist_names[i]
  cat(sprintf("%-25s | %-12.2f | %-12.2f | %-15.2f\n", 
              dist_name, 
              AIC(fit), 
              BIC(fit),
              fit$log_likelihood))
}

# 5. Get the best distribution
best_fit <- all_fits[[1]]
best_name <- dist_names[1]

cat(sprintf("\nBest fitting distribution: %s (%s: %.4f)\n", 
            best_name, best_fit$fit_criterion_name, best_fit$fit_criterion_value))

# 6. Create CDF plot for the best fit
output_dir <- "outputs"
output_file <- "poc_best_fit_single"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get the CDF function for the best distribution
best_cdf_func <- get(paste0(best_name, "_cdf"))

plot_cdf_comparison(
  sample_data = sample_data,
  fit = best_fit,
  dist_cdf = best_cdf_func,
  output_dir = output_dir,
  output_file = output_file,
  title = sprintf("Best Single Distribution Fit (Data: %s, Best: %s)", 
                  basename(input_file), best_name),
  x_label = "Value"
)

cat(sprintf("\nProcessing complete. Best fit plot saved to %s/%s.png\n", output_dir, output_file))
