# POC: Fit a 2-component mixture with mean constraint

# Initialize renv if it exists
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../pkg/rfit.R")
source("../pkg/plots/cdf_plot.R")
source("../pkg/plots/pdf_plot.R")

# Prevent Rplots.pdf from being created
if (!interactive()) {
  pdf(NULL)
}

# 1. Load data
input_file <- "inputs/whole_dataset.csv"
cat(sprintf("Loading data from: %s\n", input_file))
raw_data <- read.csv(input_file)
data <- raw_data$traffic
data <- data[!is.na(data) & data > 0]
empirical_mean <- mean(data)

cat(sprintf("Loaded %d data points.\n", length(data)))
cat(sprintf("Empirical Mean: %.4f\n", empirical_mean))

# 2. Take sample
set.seed(42)
subset_size <- min(10000, length(data))
data <- sample(data, subset_size)
cat(sprintf("Sample %d data points.\n", length(data)))

# 3. Define mixture components
dist_names <- c("lognormal", "dpln")
cat(sprintf("\nFitting mixture WITH MEAN CONSTRAINT: %s + %s\n", dist_names[1], dist_names[2]))

# 4. Fit the mixture with mean constraint
fit <- mixture_fit_mean(data, dist_names)

# 5. Output results
cat("\n=== Fit Results (Constrained) ===\n")
cat(sprintf("Log-Likelihood (penalized): %.4f\n", fit$log_likelihood))

for (i in 1:length(fit$weights)) {
  dist_name <- fit$dist_names[i]
  cat(sprintf("Component %d (%s):\n", i, dist_name))
  cat(sprintf("  Weight: %.4f\n", fit$weights[i]))
  comp_params <- fit$components[[i]]
  
  # Calculate component mean
  mean_func <- get(paste0(dist_name, "_mean"))
  c_mean <- mean_func(comp_params)
  cat(sprintf("  Component Mean: %.4f\n", c_mean))
}

mixture_theoretical_mean <- mixture_mean(fit)
cat("\n=== Mean Comparison ===\n")
cat(sprintf("Mixture Theoretical Mean: %.4f\n", mixture_theoretical_mean))
cat(sprintf("Empirical Mean:          %.4f\n", empirical_mean))
cat(sprintf("Difference:              %.4e\n", mixture_theoretical_mean - empirical_mean))

# 6. Plots
output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_cdf_comparison(
  sample_data = data,
  fit = fit,
  dist_cdf = mixture_cdf,
  output_dir = output_dir,
  output_file = "fit_2m_mean_constrained_cdf",
  title = "2-Comp Mixture Fit (Mean Constrained) - CDF",
  empirical_mean = empirical_mean,
  fitted_mean = mixture_theoretical_mean
)

plot_pdf_comparison(
  sample_data = data,
  fit = fit,
  dist_pdf = mixture_pdf,
  output_dir = output_dir,
  output_file = "fit_2m_mean_constrained_pdf",
  title = "2-Comp Mixture Fit (Mean Constrained) - PDF",
  empirical_mean = empirical_mean,
  fitted_mean = mixture_theoretical_mean
)

cat("\nProcessing complete.\n")
