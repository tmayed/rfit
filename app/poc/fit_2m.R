# POC: Fit a hardcoded 2-component mixture to CSV data

# Initialize renv if it exists
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../pkg/rfit.R")
source("../pkg/plots/cdf_plot.R")
source("../pkg/plots/pdf_plot.R")
source("../pkg/plots/diag_plot.R")

# Prevent Rplots.pdf from being created
if (!interactive()) {
  pdf(NULL)
}

# 1. Load data
input_file <- "inputs/whole_dataset.csv"
if (!file.exists(input_file)) {
  # Try without ../ if called from poc directory
  if (!file.exists(input_file)) {
     stop(sprintf("Input file not found: %s", input_file))
  }
}

input_name <- tools::file_path_sans_ext(basename(input_file))
output_dir <- file.path("outputs", paste0(input_name, "_fit_2m"))

# Create or empty the output directory
if (dir.exists(output_dir)) {
  unlink(list.files(output_dir, full.names = TRUE), recursive = TRUE)
} else {
  dir.create(output_dir, recursive = TRUE)
}

cat(sprintf("Loading data from: %s\n", input_file))
raw_data <- read.csv(input_file)

# Use 'traffic' column as requested
if ("traffic" %in% colnames(raw_data)) {
  data <- raw_data$traffic
} else {
  stop("Column 'traffic' not found in dataset")
}
data <- data[!is.na(data) & data > 0]

set.seed(42)
subset_size <- min(10000, length(data))
data <- sample(data, subset_size)

cat(sprintf("Loaded %d data points.\n", length(data)))

# 2. Define hardcoded mixture components
dist_names <- c("lognormal", "dpln")
cat(sprintf("Fitting mixture: %s + %s\n", dist_names[1], dist_names[2]))

# 3. Fit the mixture
fit <- mixture_fit(data, dist_names)

# 4. Output results to terminal
cat("\n=== Fit Results ===\n")
cat(sprintf("Log-Likelihood: %.4f\n", fit$log_likelihood))
cat(sprintf("AIC: %.4f\n", AIC(fit)))
cat(sprintf("BIC: %.4f\n", BIC(fit)))
cat("\nComponents:\n")

component_means <- numeric(length(fit$weights))
fit_results_df <- data.frame(
  distribution = character(),
  param = character(),
  value = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(fit$weights)) {
  dist_name <- fit$dist_names[i]
  cat(sprintf("Component %d (%s):\n", i, dist_name))
  cat(sprintf("  Weight: %.4f\n", fit$weights[i]))
  
  fit_results_df <- rbind(fit_results_df, data.frame(
    distribution = dist_name, param = "weight", value = fit$weights[i]
  ))
  
  comp_params <- fit$components[[i]]
  for (pname in names(comp_params)) {
    if (pname != "distribution") {
      cat(sprintf("  %s: %.4f\n", pname, comp_params[[pname]]))
      fit_results_df <- rbind(fit_results_df, data.frame(
        distribution = dist_name, param = pname, value = comp_params[[pname]]
      ))
    }
  }
  
  # Calculate component mean
  mean_func <- get(paste0(dist_name, "_mean"))
  c_mean <- mean_func(comp_params)
  component_means[i] <- c_mean
  cat(sprintf("  Component Mean: %.4f\n", c_mean))
  fit_results_df <- rbind(fit_results_df, data.frame(
    distribution = dist_name, param = "component_mean", value = c_mean
  ))
}

# 5. Save results to CSV
csv_file <- file.path(output_dir, "results.csv")
write.csv(fit_results_df, csv_file, row.names = FALSE)
cat(sprintf("\nResults exported to: %s\n", csv_file))

# 6. Calculate Mixture Mean and Empirical Mean
mixture_theoretical_mean <- sum(fit$weights * component_means)
empirical_mean <- mean(data)

cat("\n=== Mean Comparison ===\n")
cat(sprintf("Mixture Theoretical Mean: %.4f\n", mixture_theoretical_mean))
cat(sprintf("Empirical Mean:          %.4f\n", empirical_mean))

# 7. Create Plots
cat("\nGenerating plots...\n")

plot_cdf_comparison(
  sample_data = data,
  fit = fit,
  dist_cdf = mixture_cdf,
  output_dir = output_dir,
  output_file = "cdf",
  title = sprintf("2-Component Mixture Fit (Traffic Data: %s+%s)", 
                  dist_names[1], dist_names[2]),
  x_label = "Traffic",
)

plot_cdf_comparison(
  sample_data = data,
  fit = fit,
  dist_cdf = mixture_cdf,
  output_dir = output_dir,
  output_file = "cdf_log",
  title = sprintf("2-Component Mixture Fit (Log Scale) (Traffic Data: %s+%s)", 
                  dist_names[1], dist_names[2]),
  x_label = "Traffic",
  log_x = TRUE
)

plot_pdf_comparison(
  sample_data = data,
  fit = fit,
  dist_pdf = mixture_pdf,
  output_dir = output_dir,
  output_file = "pdf",
  title = sprintf("2-Component Mixture PDF Fit (Traffic Data: %s+%s)", 
                  dist_names[1], dist_names[2]),
  x_label = "Traffic",
  empirical_mean = empirical_mean,
  fitted_mean = mixture_theoretical_mean
)

plot_diag(
  sample_data = data,
  fit = fit,
  output_dir = output_dir,
  output_file = "diag"
)

cat(sprintf("Processing complete. All outputs saved to: %s/\n", output_dir))
