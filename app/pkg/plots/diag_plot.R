# Diagnostic plotting utilities for distribution fitting

#' Plot diagnostic dashboard for fitted distribution
#' @param sample_data Numeric vector of observed sample data
#' @param fit A fit object returned by one of the {dist}_fit functions
#' @param output_dir Character string path for output directory
#' @param output_file Character string path for output file (without .png extension)
#' @return NULL (saves plot to file)
#' @export
plot_diag <- function(sample_data, fit, output_dir, output_file) {
  # Validate sample_data
  sample_data <- sample_data[!is.na(sample_data)]
  if (length(sample_data) == 0) {
    stop("sample_data contains no valid values")
  }

  # Build full output path
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_file_full <- file.path(output_dir, paste0(output_file, ".png"))

  # Identify distribution and functions
  dist_name <- fit$distribution
  if (is.null(dist_name)) stop("Fit object must have a 'distribution' field")

  # Get required functions using naming convention
  pdf_func <- get(paste0(dist_name, "_pdf"))
  cdf_func <- get(paste0(dist_name, "_cdf"))
  q_func   <- get(paste0(dist_name, "_quantile"))

  # Prepare data
  n <- length(sample_data)
  sorted_data <- sort(sample_data)
  
  # Probabilities for P-P and Q-Q
  # Using (i - 0.5) / n for plotting positions
  p_points <- (seq_len(n) - 0.5) / n
  
  # Open PNG device
  png(output_file_full, width = 1200, height = 1000, res = 120)
  
  # 2x2 layout
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2))

  # 1. PDF Plot (Histogram + Fitted Line)
  h <- hist(sample_data, breaks = "FD", plot = FALSE)
  x_range <- range(sample_data)
  x_grid <- seq(x_range[1], x_range[2], length.out = 200)
  y_fit <- pdf_func(x_grid, fit)
  
  plot(h, freq = FALSE, col = "lightgrey", border = "white",
       main = paste("PDF:", dist_name), xlab = "Value", ylab = "Density")
  lines(x_grid, y_fit, col = "red", lwd = 2)

  # 2. CDF Plot (Empirical vs Fitted)
  plot(sorted_data, p_points, type = "s", col = "blue", lwd = 1,
       main = paste("CDF:", dist_name), xlab = "Value", ylab = "Cumulative Probability")
  y_cdf <- cdf_func(x_grid, fit)
  lines(x_grid, y_cdf, col = "red", lwd = 2)
  legend("bottomright", legend = c("Empirical", "Fitted"), 
         col = c("blue", "red"), lwd = 2, bty = "n")

  # 3. Q-Q Plot (Quantile-Quantile)
  theo_quantiles <- q_func(p_points, fit)
  
  plot(theo_quantiles, sorted_data, pch = 20, col = rgb(0,0,1,0.5),
       main = "Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles")
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  
  # 4. P-P Plot (Probability-Probability)
  theo_probs <- cdf_func(sorted_data, fit)
  
  plot(p_points, theo_probs, pch = 20, col = rgb(0,0,1,0.5),
       main = "P-P Plot", xlab = "Empirical Probabilities", ylab = "Theoretical Probabilities")
  abline(0, 1, col = "red", lwd = 2, lty = 2)

  # Reset layout and close device
  dev.off()
  
  message(sprintf("Diagnostic plot saved to: %s", output_file_full))
}
