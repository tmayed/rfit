# Plotting utilities for Mixture Distributions - Diagnostic Version

#' Plot mixture and components separately
#' @param sample_data Numeric vector of observed sample data
#' @param fit Fitted mixture object
#' @param output_dir Character string path for output directory
#' @param prefix Prefix for the three output files
#' @param transform "none" for standard PDF, "log" for density of log-data
#' @export
plot_mixture_diagnostic <- function(sample_data, fit, 
                                    output_dir, prefix,
                                    transform = c("none", "log")) {
  
  transform <- match.arg(transform)
  library(dplyr)
  library(ggplot2)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Data transformation
  plot_x_data <- if (transform == "log") log(sample_data[sample_data > 0]) else sample_data
  x_range <- range(plot_x_data, na.rm = TRUE)
  x_points_transformed <- seq(x_range[1], x_range[2], length.out = 1000)
  x_points_original <- if (transform == "log") exp(x_points_transformed) else x_points_transformed

  # Helper to get transformed density
  # g(y) = f(exp(y)) * exp(y)
  get_density <- function(x_orig, x_trans, dist_name, comp_fit, weight = 1.0) {
    pdf_func <- get(paste0(dist_name, "_pdf"))
    d <- pdf_func(x_orig, comp_fit) * weight
    if (transform == "log") {
      d <- d * x_orig
    }
    d
  }

  # 1. Individual Components (Unweighted)
  for (i in seq_along(fit$components)) {
    dist_name <- fit$dist_names[i]
    y_vals <- get_density(x_points_original, x_points_transformed, dist_name, fit$components[[i]], weight = 1.0)
    
    df <- data.frame(x = x_points_transformed, y = y_vals)
    
    p <- ggplot() +
      geom_density(data = data.frame(x = plot_x_data), aes(x = x), fill = "grey90", color = "grey80") +
      geom_line(data = df, aes(x = x, y = y), color = "blue", linewidth = 1) +
      labs(
        title = sprintf("Component %d: %s (Unweighted Pure PDF)", i, dist_name),
        subtitle = if(transform == "log") "Density of Log-Data" else "Standard PDF",
        x = if(transform == "log") "Log(Value)" else "Value",
        y = "Density"
      ) +
      theme_minimal()
    
    ggsave(file.path(output_dir, sprintf("%s_comp%d.png", prefix, i)), plot = p, width = 8, height = 5)
  }

  # 2. Total Mixture
  mixture_y <- numeric(length(x_points_original))
  for (i in seq_along(fit$components)) {
    mixture_y <- mixture_y + get_density(x_points_original, x_points_transformed, 
                                         fit$dist_names[i], fit$components[[i]], 
                                         weight = fit$weights[i])
  }
  
  mix_df <- data.frame(x = x_points_transformed, y = mixture_y)
  
  p_mix <- ggplot() +
    geom_density(data = data.frame(x = plot_x_data), aes(x = x), fill = "grey90", color = "grey80") +
    geom_line(data = mix_df, aes(x = x, y = y), color = "black", linewidth = 1.2) +
    labs(
      title = "Total Mixture Fit",
      subtitle = if(transform == "log") "Density of Log-Data" else "Standard PDF",
      x = if(transform == "log") "Log(Value)" else "Value",
      y = "Density"
    ) +
    theme_minimal()
  
  ggsave(file.path(output_dir, sprintf("%s_mixture.png", prefix)), plot = p_mix, width = 8, height = 5)
  
  message(sprintf("Diagnostic plots saved to %s with prefix %s", output_dir, prefix))
}
