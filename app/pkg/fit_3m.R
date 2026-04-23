# 3-Component Mixture Fitting Orchestrator

#' @title fit_all_3m
#' @description Fits all combinations of 3-component mixtures of available distributions to the provided data 
#'              and returns them ordered by quality of fit.
#' @param data A numeric vector of positive data points.
#' @param dist_names Optional character vector of distribution names to use for mixtures. If NULL, all registered distributions are used.
#' @param criterion The criterion to use for ordering fits. Options are "AIC" (default) or "BIC".
#' @return A list of mixture fit objects, ordered from best (lowest criterion) to worst.
#' @export
fit_all_3m <- function(data, dist_names = NULL, criterion = c("AIC", "BIC")) {
  criterion <- match.arg(criterion)
  
  # Ensure data is positive
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 15) {
    stop("Need at least 15 valid positive data points for 3-component mixture fitting.")
  }

  # List of distribution names available for mixtures
  if (is.null(dist_names)) {
    distributions <- names(.DIST_REGISTRY)
  } else {
    distributions <- dist_names
  }

  # Generate all unique triples with NO repeats (distinct distributions)
  # combn returns a matrix where each column is a combination
  triples_mat <- combn(distributions, 3)
  triples <- as.data.frame(t(triples_mat), stringsAsFactors = FALSE)
  colnames(triples) <- c("dist1", "dist2", "dist3")

  fits <- list()
  
  for (i in 1:nrow(triples)) {
    dists <- as.character(triples[i, ])
    mix_name <- paste(dists, collapse = "+")
    
    cat(sprintf("[%d/%d] Fitting 3-component mixture: %s...\n", i, nrow(triples), mix_name))
    
    # Attempt to fit the mixture
    fit <- tryCatch({
      mixture_fit(data, dists)
    }, error = function(e) {
      warning(sprintf("Failed to fit mixture %s: %s", mix_name, e$message))
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
        # Manual calculation if AIC/BIC fails
        ll <- fit$log_likelihood
        df <- length(fit$weights) - 1
        for (comp in fit$components) {
          df <- df + (length(comp) - 1)
        }
        
        if (criterion == "AIC") {
          2 * df - 2 * ll
        } else {
          df * log(length(data)) - 2 * ll
        }
      })
      
      # Store identification and criterion
      fit$distribution_mixture <- mix_name
      fit$fit_criterion_value <- val
      fit$fit_criterion_name <- criterion
      fits[[mix_name]] <- fit
    }
  }
  
  if (length(fits) == 0) {
    stop("No 3-component mixtures were successfully fitted.")
  }
  
  # Order fits by the criterion value (lowest is best)
  ordered_indices <- order(sapply(fits, function(f) f$fit_criterion_value))
  ordered_fits <- fits[ordered_indices]
  
  return(ordered_fits)
}
