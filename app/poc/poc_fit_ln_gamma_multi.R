# POC: Fit the dedicated lognormal-gamma mixture to multiple CSV data files

# Configuration
# Set to NULL to automatically detect all available dates in the inputs/ directory
date_str_list <- NULL 
ln_weight <- 0.25
mean_fit <- TRUE
sample_size <- 2000

# Initialize renv
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../pkg/rfit.R")

# Function to detect dates from filenames of the form router_traffic_yyyy_mm_dd.csv
get_available_dates <- function() {
  dirs <- c("inputs", "poc/inputs")
  all_files <- unlist(lapply(dirs, function(d) {
    if (dir.exists(d)) list.files(d, pattern = "^router_traffic_.*\\.csv$") else NULL
  }))
  
  if (length(all_files) == 0) return(character(0))
  
  # Extract yyyy_mm_dd using regex
  dates <- gsub("^router_traffic_(\\d{4}_\\d{2}_\\d{2})\\.csv$", "\\1", all_files)
  
  # Return sorted unique dates
  return(sort(unique(dates)))
}

# Auto-detect if list is NULL
if (is.null(date_str_list)) {
  date_str_list <- get_available_dates()
  if (length(date_str_list) == 0) {
    stop("No data files found and date_str_list is NULL.")
  }
  cat(sprintf("Auto-detected %d dates: %s\n", length(date_str_list), paste(date_str_list, collapse = ", ")))
}

# Dedicated CDF for the fit_ln_gamma() result object.
fit_ln_gamma_cdf <- function(x, fit) {
  ln_cdf <- plnorm(
    x,
    meanlog = fit$components$lognormal$mu,
    sdlog = fit$components$lognormal$sigma
  )

  gamma_cdf <- pgamma(
    x,
    shape = fit$components$gamma$shape,
    scale = fit$components$gamma$scale
  )

  fit$weights["lognormal"] * ln_cdf +
    fit$weights["gamma"] * gamma_cdf
}

# Prevent Rplots.pdf from being created
if (!interactive()) {
  pdf(NULL)
}

all_results_df <- data.frame()

for (date_str in date_str_list) {
  cat(sprintf("\n--- Processing Date: %s ---\n", date_str))
  
  # 1. Load data
  input_file <- sprintf("inputs/router_traffic_%s.csv", date_str)
  if (!file.exists(input_file)) {
    input_file <- sprintf("poc/inputs/router_traffic_%s.csv", date_str)
  }
  if (!file.exists(input_file)) {
    cat(sprintf("Warning: Input file not found for %s, skipping.\n", date_str))
    next
  }

  data_raw <- read.csv(input_file)$traffic
  data_clean <- data_raw[!is.na(data_raw) & data_raw > 0]

  # 2. Take a subset robustly
  set.seed(42)
  subset_size <- min(sample_size, length(data_clean))
  sample_data <- sample(data_clean, subset_size)

  cat(sprintf("Loaded %d observations, using subset of %d\n", 
              length(data_clean), subset_size))

  # 3. Fit the dedicated lognormal-gamma mixture
  cat("Fitting lognormal-gamma mixture...\n")
  best_fit <- fit_ln_gamma(sample_data, w=ln_weight, mean_fit = mean_fit)

  # 4. Calculate Tracking Metrics
  norm_loglik <- best_fit$log_likelihood / length(sample_data)
  
  sorted_data <- sort(sample_data)
  theoretical_cdf <- fit_ln_gamma_cdf(sorted_data, best_fit)
  empirical_cdf <- seq_along(sorted_data) / length(sorted_data)
  ks_stat <- max(abs(theoretical_cdf - empirical_cdf))

  # 5. Collect results
  fit_results_df <- data.frame(
    date_str = character(),
    distribution = character(),
    param = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )

  # Helper function to add rows
  add_res <- function(df, dist, param, val) {
    rbind(df, data.frame(date_str = date_str, distribution = dist, param = param, value = val))
  }

  # Metrics
  fit_results_df <- add_res(fit_results_df, "metrics", "bic", BIC(best_fit))
  fit_results_df <- add_res(fit_results_df, "metrics", "norm_loglik", norm_loglik)
  fit_results_df <- add_res(fit_results_df, "metrics", "ks_stat", ks_stat)
  fit_results_df <- add_res(fit_results_df, "metrics", "sample_size", length(sample_data))
  fit_results_df <- add_res(fit_results_df, "metrics", "convergence", best_fit$convergence)

  # Lognormal component
  fit_results_df <- add_res(fit_results_df, "lognormal", "weight", best_fit$weights["lognormal"])
  for (p in names(best_fit$components$lognormal)) {
    fit_results_df <- add_res(fit_results_df, "lognormal", p, best_fit$components$lognormal[[p]])
  }
  comp_mean_ln <- lognormal_2p_mean(best_fit$components$lognormal)
  ln_exp_term <- best_fit$components$lognormal$mu + (best_fit$components$lognormal$sigma^2) / 2
  fit_results_df <- add_res(fit_results_df, "lognormal", "component_mean", comp_mean_ln)
  fit_results_df <- add_res(fit_results_df, "lognormal", "log_mean_term", ln_exp_term)

  # Gamma component
  fit_results_df <- add_res(fit_results_df, "gamma", "weight", best_fit$weights["gamma"])
  for (p in names(best_fit$components$gamma)) {
    fit_results_df <- add_res(fit_results_df, "gamma", p, best_fit$components$gamma[[p]])
  }
  comp_mean_g <- gamma_2p_mean(best_fit$components$gamma)
  gamma_exp_term <- log(best_fit$components$gamma$shape * best_fit$components$gamma$scale)
  fit_results_df <- add_res(fit_results_df, "gamma", "component_mean", comp_mean_g)
  fit_results_df <- add_res(fit_results_df, "gamma", "log_mean_term", gamma_exp_term)

  # Overall Mean
  mixture_theoretical_mean <- sum(best_fit$weights * c(comp_mean_ln, comp_mean_g))
  empirical_mean <- mean(sample_data)
  fit_results_df <- add_res(fit_results_df, "mixture", "theoretical_mean", mixture_theoretical_mean)
  fit_results_df <- add_res(fit_results_df, "mixture", "empirical_mean", empirical_mean)

  # Append to master dataframe
  all_results_df <- rbind(all_results_df, fit_results_df)
}

# 6. Save consolidated results
output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

csv_file <- file.path(output_dir, "multi_fit_results.csv")
write.csv(all_results_df, csv_file, row.names = FALSE)

cat(sprintf("\nMulti-fit processing complete. Results saved to: %s\n", csv_file))
