# Plotting utilities for distribution fitting

#' Plot CDF comparison between sample data and fitted distribution
#' @param sample_data Numeric vector of observed sample data
#' @param fitted_cdf Function that takes quantiles and returns CDF values
#'                     (e.g., a function closure with fitted parameters)
#' @param output_dir Character string path for output directory
#' @param output_file Character string path for output file (without .png extension)
#' @param title Plot title (default: "Empirical vs Fitted CDF")
#' @param x_label X-axis label (default: "Value")
#' @param y_label Y-axis label (default: "Cumulative Probability")
#' @return ggplot2 plot object
#' @examples
#' # Example usage:
#' data <- c(1, 2, 3, 4, 5)
#' fit <- lognormal_fit(data)
#' cdf_func <- function(q) plnorm(q, meanlog = fit$mu, sdlog = fit$sigma)
#' plot_cdf_comparison(sample_data = data, fitted_cdf = cdf_func,
#'                     output_dir = "poc/outputs", output_file = "lognormal_cdf")
#' @export
plot_cdf_comparison <- function(sample_data, fit, dist_cdf,
                                output_dir, output_file,
                                title = "Empirical vs Fitted CDF",
                                x_label = "Value",
                                y_label = "Cumulative Probability",
                                empirical_mean = NULL,
                                fitted_mean = NULL,
                                log_x = FALSE,
                                log_y = FALSE,
                                complementary = FALSE) {

  # Validate sample_data
  sample_data <- sample_data[!is.na(sample_data)]
  if (length(sample_data) == 0) {
    stop("sample_data contains no valid values")
  }

  # Build full output path
  output_file_full <- file.path(output_dir, paste0(output_file, ".png"))

  # Create empirical CDF data
  n <- length(sample_data)
  sorted_data <- sort(sample_data)
  empirical_probs <- seq_len(n) / n
  
  if (complementary) {
    empirical_probs <- 1 - (seq_len(n) - 0.5) / n # Standard plotting position for CCDF
  }

  empirical_df <- data.frame(
    x = sorted_data,
    y = empirical_probs,
    type = "Empirical"
  )

  # Create fitted CDF data
  x_range <- range(sample_data)
  if (log_x) {
    x_range[1] <- max(x_range[1], 1e-6)
    x_points <- exp(seq(log(x_range[1]), log(x_range[2]), length.out = 1000))
  } else {
    x_points <- seq(x_range[1], x_range[2], length.out = 1000)
  }
  
  fitted_probs <- dist_cdf(x=x_points, fit=fit)
  if (complementary) {
    fitted_probs <- 1 - fitted_probs
  }

  fitted_df <- data.frame(
    x = x_points,
    y = fitted_probs,
    type = "Fitted"
  )
  # Filter out non-finite values (should be fewer now)
  fitted_df <- fitted_df[is.finite(fitted_df$x) & is.finite(fitted_df$y), ]
  if (complementary) {
    fitted_df <- fitted_df[fitted_df$y > 0, ] # Avoid log(0)
  }

  # Prepare data for vertical lines (means)
  means_df <- data.frame(xintercept = numeric(0), type = character(0))
  
  if (!is.null(empirical_mean)) {
    means_df <- rbind(means_df, data.frame(
      xintercept = empirical_mean,
      type = "Empirical Mean"
    ))
  }
  
  if (!is.null(fitted_mean)) {
    means_df <- rbind(means_df, data.frame(
      xintercept = fitted_mean,
      type = "Fitted Mean"
    ))
  }
  
  # Filter out non-finite means
  if (nrow(means_df) > 0) {
    means_df <- means_df[is.finite(means_df$xintercept), ]
  }

  # Now load libraries for plotting
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  # Define all possible levels and their aesthetics
  all_colors <- c(
    "Empirical" = "blue", 
    "Fitted" = "red",
    "Empirical Mean" = "darkblue",
    "Fitted Mean" = "darkred"
  )
  all_linetypes <- c(
    "Empirical" = "solid", 
    "Fitted" = "solid",
    "Empirical Mean" = "dashed",
    "Fitted Mean" = "dashed"
  )

  # Create plot
  p <- ggplot() +
    geom_step(data = empirical_df, aes(x = x, y = y, color = type, linetype = type), linewidth = 1) +
    geom_line(data = fitted_df, aes(x = x, y = y, color = type, linetype = type), linewidth = 1)
  
  if (nrow(means_df) > 0) {
    p <- p + geom_vline(data = means_df, aes(xintercept = xintercept, color = type, linetype = type), 
                        linewidth = 0.8)
  }

  p <- p +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = "Legend",
      linetype = "Legend"
    ) +
    coord_cartesian(xlim = x_range) + 
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_color_manual(values = all_colors) +
    scale_linetype_manual(values = all_linetypes)

  if (log_x) {
    p <- p + scale_x_log10()
  }
  if (log_y) {
    p <- p + scale_y_log10()
  }

  # Save plot
  ggsave(output_file_full, plot = p, width = 10, height = 6, dpi = 150)

  message(sprintf("Plot saved to: %s", output_file_full))

  return(p)
}