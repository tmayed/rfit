# Plotting utilities for distribution fitting (PDF)

#' Plot PDF comparison between sample data and fitted distribution
#' @param sample_data Numeric vector of observed sample data
#' @param fit Fitted distribution object
#' @param dist_pdf Function that takes x and fit, and returns PDF values
#' @param output_dir Character string path for output directory
#' @param output_file Character string path for output file (without .png extension)
#' @param title Plot title (default: "Empirical vs Fitted PDF")
#' @param x_label X-axis label (default: "Value")
#' @param y_label Y-axis label (default: "Density")
#' @param empirical_mean Optional pre-calculated empirical mean
#' @param fitted_mean Optional pre-calculated fitted mean
#' @return ggplot2 plot object
#' @export
plot_pdf_comparison <- function(sample_data, fit, dist_pdf,
                                output_dir, output_file,
                                title = "Empirical vs Fitted PDF",
                                x_label = "Value",
                                y_label = "Density",
                                empirical_mean = NULL,
                                fitted_mean = NULL) {

  # Validate sample_data
  sample_data <- sample_data[!is.na(sample_data)]
  if (length(sample_data) == 0) {
    stop("sample_data contains no valid values")
  }

  # Build full output path
  output_file_full <- file.path(output_dir, paste0(output_file, ".png"))

  # Range for plotting
  x_range <- range(sample_data)
  # Increase resolution for better curves
  x_points <- seq(x_range[1], x_range[2], length.out = 1000)
  
  # Calculate fitted PDF
  fitted_dens <- dist_pdf(x=x_points, fit=fit)
  
  fitted_df <- data.frame(
    x = x_points,
    y = fitted_dens,
    type = "Fitted"
  )
  # Filter out non-finite values
  fitted_df <- fitted_df[is.finite(fitted_df$x) & is.finite(fitted_df$y), ]

  # Calculate means if not provided
  if (is.null(empirical_mean)) {
    empirical_mean <- mean(sample_data)
  }
  
  if (is.null(fitted_mean)) {
    # Distribution mean using naming convention
    dist_name <- fit$distribution
    dist_mean <- NA
    
    if (!is.null(dist_name)) {
      mean_func_name <- paste0(dist_name, "_mean")
      
      # Try direct lookup first
      dist_mean <- tryCatch({
        if (dist_name == "mixture") {
          mixture_mean(fit)
        } else {
          func <- get(mean_func_name)
          func(fit)
        }
      }, error = function(e) {
        tryCatch({
          func <- get(mean_func_name, envir = .GlobalEnv)
          func(fit)
        }, error = function(e2) {
          NA
        })
      })
    }
    fitted_mean <- dist_mean
  }
  
  # Create a data frame for vertical lines
  means_df <- data.frame(
    xintercept = c(empirical_mean, fitted_mean),
    type = c("Empirical Mean", "Fitted Mean")
  )
  # Filter out non-finite means
  means_df <- means_df[is.finite(means_df$xintercept), ]

  # Now load libraries for plotting
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  # Create plot
  # We use geom_density for empirical and geom_line for fitted
  p <- ggplot() +
    geom_density(data = data.frame(x = sample_data), aes(x = x, color = "Empirical"), linewidth = 1) +
    geom_line(data = fitted_df, aes(x = x, y = y, color = "Fitted"), linewidth = 1) +
    geom_vline(data = means_df, aes(xintercept = xintercept, color = type), 
               linetype = "dashed", linewidth = 0.8) +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = "Legend"
    ) +
    coord_cartesian(xlim = x_range) + # Focus on data range
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_color_manual(values = c(
      "Empirical" = "blue", 
      "Fitted" = "red",
      "Empirical Mean" = "darkblue",
      "Fitted Mean" = "darkred"
    ))

  # Save plot
  ggsave(output_file_full, plot = p, width = 10, height = 6, dpi = 150)

  message(sprintf("Plot saved to: %s", output_file_full))

  return(p)
}
