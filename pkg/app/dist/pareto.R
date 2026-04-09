# Pareto Distribution Module
# All functions related to Pareto distribution fitting

#' Log-likelihood function for Pareto distribution
#' @param data Numeric vector of observations
#' @param shape Shape parameter (alpha) (default: NULL, will be estimated)
#' @param scale Scale parameter (xm, minimum value) (default: NULL, will be estimated)
#' @return Log-likelihood value
#' @export
pareto_log_likelihood <- function(data, shape = NULL, scale = NULL) {
  # Remove NA values
  data <- data[!is.na(data)]
  data <- data[data > 0]  # Pareto is only defined for positive values

  if (length(data) == 0) {
    return(-Inf)
  }

  # If parameters not provided, estimate them from data
  if (is.null(scale)) {
    scale <- min(data)
  }
  if (is.null(shape)) {
    shape <- length(data) / sum(log(data / scale))
  }

  # Compute log-likelihood
  if (shape <= 0 || scale <= 0) {
    return(-Inf)
  }

  if (any(data < scale)) {
    return(-Inf)
  }

  log_lik <- sum(log(shape) + log(scale) - (shape + 1) * log(data))
  return(log_lik)
}

#' Fit Pareto distribution using maximum likelihood estimation
#' @param data Numeric vector of observations
#' @return List containing fitted parameters (shape, scale) and log-likelihood
#' @export
pareto_fit <- function(data) {
  # Remove NA values
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) == 0) {
    stop("No valid positive data points for fitting")
  }

  if (length(data) < 2) {
    stop("Need at least 2 data points for fitting")
  }

  # MLE estimates for Pareto distribution
  scale_hat <- min(data)
  shape_hat <- length(data) / sum(log(data / scale_hat))

  if (shape_hat <= 0) {
    stop("Invalid shape parameter estimated")
  }

  log_lik <- pareto_log_likelihood(data, shape_hat, scale_hat)

  return(list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "pareto"
  ))
}

#' Fit truncated Pareto distribution using maximum likelihood estimation
#' @param data Numeric vector of observations
#' @param lower Lower truncation bound (default: NULL, uses min(data))
#' @param upper Upper truncation bound (default: Inf)
#' @return List containing fitted parameters and log-likelihood
#' @export
pareto_fit_truncated <- function(data, lower = NULL, upper = Inf) {
  # Remove NA values
  data <- data[!is.na(data)]

  # Set lower bound if not provided
  if (is.null(lower)) {
    lower <- min(data)
  }

  # Filter data within truncation bounds
  data <- data[data >= lower & data <= upper]

  if (length(data) == 0) {
    stop("No valid data points within truncation bounds")
  }

  if (length(data) < 2) {
    stop("Need at least 2 data points for fitting")
  }

  # Negative log-likelihood for truncated Pareto
  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])  # Ensure shape > 0

    if (shape <= 0) {
      return(Inf)
    }

    # For truncated Pareto, scale is the lower bound of the distribution
    # and we condition on [lower, upper]
    scale <- lower  # Scale is fixed at lower bound for truncation
    n <- length(data)

    # PDF: alpha * xm^alpha / x^(alpha+1)
    # CDF: 1 - (xm/x)^alpha
    ll <- sum(log(shape) + shape * log(scale) - (shape + 1) * log(data))

    # Normalization for truncation
    # P(a < X < b) = F(b) - F(a) = (xm/a)^alpha - (xm/b)^alpha
    if (is.finite(upper)) {
      norm <- log((scale / lower)^shape - (scale / upper)^shape)
    } else {
      norm <- log((scale / lower)^shape)
    }

    ll <- ll - n * norm

    return(-ll)
  }

  # Initial estimates
  scale_hat <- lower
  shape_hat <- length(data) / sum(log(data / scale_hat))

  # Optimize
  fit_result <- optim(
    c(log(shape_hat)),
    neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  shape_hat <- exp(fit_result$par[1])
  scale_hat <- lower  # Scale is fixed at lower bound

  log_lik <- -neg_log_likelihood(c(log(shape_hat)))

  return(list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    distribution = "truncated_pareto"
  ))
}

#' Probability density function for fitted Pareto distribution
#' @param x Numeric vector of points at which to evaluate density
#' @param fit List returned by pareto_fit
#' @return Density values
#' @export
pareto_pdf <- function(x, fit) {
  if (any(x < fit$scale)) {
    warning("Some x values are below the scale parameter")
  }
  fit$shape * (fit$scale^fit$shape) / x^(fit$shape + 1)
}

#' Cumulative distribution function for fitted Pareto distribution
#' @param q Numeric vector of quantiles
#' @param fit List returned by pareto_fit
#' @return Cumulative probabilities
#' @export
pareto_cdf <- function(q, fit) {
  pmax(0, pmin(1, 1 - (fit$scale / q)^fit$shape))
}

#' Quantile function for fitted Pareto distribution
#' @param p Numeric vector of probabilities
#' @param fit List returned by pareto_fit
#' @return Quantile values
#' @export
pareto_quantile <- function(p, fit) {
  fit$scale / (1 - p)^(1 / fit$shape)
}

#' Generate random samples from fitted Pareto distribution
#' @param n Number of samples to generate
#' @param fit List returned by pareto_fit
#' @return Numeric vector of random samples
#' @export
pareto_rand <- function(n, fit) {
  # Using inverse transform sampling
  u <- runif(n)
  fit$scale / (1 - u)^(1 / fit$shape)
}

#' Mean function for fitted Pareto distribution
#' @param fit List returned by pareto_fit
#' @return Mean value (returns Inf if shape <= 1)
#' @export
pareto_mean <- function(fit) {
  if (fit$shape <= 1) {
    return(Inf)
  }
  fit$scale * fit$shape / (fit$shape - 1)
}

#' Standard deviation function for fitted Pareto distribution
#' @param fit List returned by pareto_fit
#' @return Standard deviation value (returns Inf if shape <= 2)
#' @export
pareto_std <- function(fit) {
  if (fit$shape <= 2) {
    return(Inf)
  }
  fit$scale * sqrt(fit$shape / ((fit$shape - 1)^2 * (fit$shape - 2)))
}
