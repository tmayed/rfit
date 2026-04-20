# Distribution Fitting Orchestrator

#' @title fit_all
#' @description Fits all available distributions to the provided data and returns them ordered by quality of fit.
#' @param data A numeric vector of positive data points.
#' @param criterion The criterion to use for ordering fits. Options are "AIC" (default) or "BIC".
#' @return A list of fit objects, ordered from best (lowest criterion) to worst.
#' @export
fit_all <- function(data, criterion = c("AIC", "BIC")) {
  criterion <- match.arg(criterion)
  
  # Ensure data is positive for distributions that require it
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 2) {
    stop("Need at least 2 valid positive data points for fitting.")
  }

  # List of distribution names available in the package
  # These correspond to the {dist}_fit functions
  distributions <- c(
    "lognormal", "normal", "weibull", "pareto", "gb2", 
    "dpln", "kappa4", "bradford", "fisk",
    "johnsonsu", "johnsonsb", "johnsonsl"
  )

  fits <- list()
  
  for (dist in distributions) {
    fit_func_name <- paste0(dist, "_fit")
    
    # Check if the fit function exists
    if (!exists(fit_func_name, mode = "function")) {
      warning(sprintf("Fit function %s not found. Skipping.", fit_func_name))
      next
    }
    
    fit_func <- get(fit_func_name)
    
    # Attempt to fit and handle potential errors during optimization
    fit <- tryCatch({
      fit_func(data)
    }, error = function(e) {
      warning(sprintf("Failed to fit %s: %s", dist, e$message))
      NULL
    })
    
    if (!is.null(fit)) {
      # Calculate the chosen criterion
      val <- tryCatch({
        if (criterion == "AIC") {
          AIC(fit)
        } else {
          BIC(fit)
        }
      }, error = function(e) {
        # Fallback if logLik method is missing or fails
        # Assuming most have 2 parameters if logLik is missing
        k <- if (is.null(attr(fit, "df"))) 2 else attr(fit, "df")
        ll <- fit$log_likelihood
        if (criterion == "AIC") {
          2 * k - 2 * ll
        } else {
          k * log(length(data)) - 2 * ll
        }
      })
      
      # Store the criterion value in the fit object for sorting
      fit$fit_criterion_value <- val
      fit$fit_criterion_name <- criterion
      fits[[dist]] <- fit
    }
  }
  
  if (length(fits) == 0) {
    stop("No distributions were successfully fitted.")
  }
  
  # Order fits by the criterion value (lowest is best)
  ordered_indices <- order(sapply(fits, function(f) f$fit_criterion_value))
  ordered_fits <- fits[ordered_indices]
  
  return(ordered_fits)
}
