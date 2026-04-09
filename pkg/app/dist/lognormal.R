# Lognormal Distribution Module
# All functions related to lognormal distribution fitting

#' Log-likelihood function for lognormal distribution
#' @param data Numeric vector of observations
#' @param mu Mean of log-transformed data (default: NULL, will be estimated)
#' @param sigma Standard deviation of log-transformed data (default: NULL, will be estimated)
#' @return Log-likelihood value
#' @export
lognormal_log_likelihood <- function(data, mu = NULL, sigma = NULL) {
  # Remove NA values
  data <- data[!is.na(data)]
  data <- data[data > 0]  # Lognormal is only defined for positive values

  if (length(data) == 0) {
    return(-Inf)
  }

  # If parameters not provided, estimate them from data
  if (is.null(mu)) {
    mu <- mean(log(data))
  }
  if (is.null(sigma)) {
    sigma <- sd(log(data))
  }

  # Compute log-likelihood
  log_lik <- sum(dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE))
  return(log_lik)
}

#' Fit lognormal distribution using maximum likelihood estimation
#' @param data Numeric vector of observations
#' @return List containing fitted parameters (mu, sigma) and log-likelihood
#' @export
lognormal_fit <- function(data) {
  # Remove NA values
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) == 0) {
    stop("No valid positive data points for fitting")
  }

  if (length(data) < 2) {
    stop("Need at least 2 data points for fitting")
  }

  # MLE estimates for lognormal distribution
  log_data <- log(data)
  mu_hat <- mean(log_data)
  sigma_hat <- sqrt(sum((log_data - mu_hat)^2) / length(data))

  log_lik <- lognormal_log_likelihood(data, mu_hat, sigma_hat)

  return(list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "lognormal"
  ))
}

#' Fit truncated lognormal distribution using maximum likelihood estimation
#' @param data Numeric vector of observations
#' @param lower Lower truncation bound (default: 0)
#' @param upper Upper truncation bound (default: Inf)
#' @return List containing fitted parameters and log-likelihood
#' @export
lognormal_fit_truncated <- function(data, lower = 0, upper = Inf) {
  # Remove NA values
  data <- data[!is.na(data)]

  # Filter data within truncation bounds
  data <- data[data >= lower & data <= upper]

  if (length(data) == 0) {
    stop("No valid data points within truncation bounds")
  }

  if (length(data) < 2) {
    stop("Need at least 2 data points for fitting")
  }

  # Negative log-likelihood for truncated lognormal
  neg_log_likelihood <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])  # Ensure sigma > 0

    if (sigma <= 0 || is.na(mu) || is.nan(mu)) {
      return(Inf)
    }

    # Log-likelihood for truncated lognormal
    # f(x | a, b) = f(x) / (F(b) - F(a))
    ll <- sum(dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE))
    ll <- ll - length(data) * log(plnorm(upper, meanlog = mu, sdlog = sigma) -
                                plnorm(lower, meanlog = mu, sdlog = sigma))

    return(-ll)
  }

  # Initial estimates from standard MLE
  log_data <- log(data)
  initial_mu <- mean(log_data)
  initial_sigma <- sd(log_data)

  # Optimize
  fit_result <- optim(
    c(initial_mu, log(initial_sigma)),
    neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  mu_hat <- fit_result$par[1]
  sigma_hat <- exp(fit_result$par[2])

  log_lik <- -neg_log_likelihood(c(mu_hat, log(sigma_hat)))

  return(list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    distribution = "truncated_lognormal"
  ))
}

#' Probability density function for fitted lognormal distribution
#' @param x Numeric vector of points at which to evaluate density
#' @param fit List returned by lognormal_fit
#' @return Density values
#' @export
lognormal_pdf <- function(x, fit) {
  dlnorm(x, meanlog = fit$mu, sdlog = fit$sigma)
}

#' Cumulative distribution function for fitted lognormal distribution
#' @param q Numeric vector of quantiles
#' @param fit List returned by lognormal_fit
#' @return Cumulative probabilities
#' @export
lognormal_cdf <- function(q, fit) {
  plnorm(q, meanlog = fit$mu, sdlog = fit$sigma)
}

#' Quantile function for fitted lognormal distribution
#' @param p Numeric vector of probabilities
#' @param fit List returned by lognormal_fit
#' @return Quantile values
#' @export
lognormal_quantile <- function(p, fit) {
  qlnorm(p, meanlog = fit$mu, sdlog = fit$sigma)
}

#' Generate random samples from fitted lognormal distribution
#' @param n Number of samples to generate
#' @param fit List returned by lognormal_fit
#' @return Numeric vector of random samples
#' @export
lognormal_rand <- function(n, fit) {
  rlnorm(n, meanlog = fit$mu, sdlog = fit$sigma)
}

#' Mean function for fitted lognormal distribution
#' @param fit List returned by lognormal_fit
#' @return Mean value
#' @export
lognormal_mean <- function(fit) {
  exp(fit$mu + fit$sigma^2 / 2)
}

#' Standard deviation function for fitted lognormal distribution
#' @param fit List returned by lognormal_fit
#' @return Standard deviation value
#' @export
lognormal_std <- function(fit) {
  sqrt((exp(fit$sigma^2) - 1) * exp(2 * fit$mu + fit$sigma^2))
}
