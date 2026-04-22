# Rayleigh Distribution Module (1-parameter: scale, location fixed at 0)

# -------------------------------
# Log-likelihood
# -------------------------------
rayleigh_1p_log_likelihood <- function(data, sigma) {
  data <- data[!is.na(data) & data >= 0]
  if (sigma <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  # PDF: x/sigma^2 * exp(-x^2 / (2*sigma^2))
  # LogPDF: log(x) - 2*log(sigma) - x^2 / (2*sigma^2)
  sum(log(data) - 2 * log(sigma) - (data^2) / (2 * sigma^2))
}

# -------------------------------
# Fit (MLE - Closed Form)
# -------------------------------
rayleigh_1p_fit <- function(data) {
  data <- data[!is.na(data) & data >= 0]
  if (length(data) < 1) stop("Need at least 1 valid data point")

  # Closed form MLE for sigma: sqrt(mean(x^2)/2)
  sigma_hat <- sqrt(mean(data^2) / 2)

  log_lik <- rayleigh_1p_log_likelihood(data, sigma_hat)

  result <- list(
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "rayleigh_1p"
  )

  class(result) <- "rayleigh_1p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
rayleigh_1p_log_likelihood_truncated <- function(data, sigma, lower, upper) {
  if (sigma <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  ll_base <- sum(log(data) - 2 * log(sigma) - (data^2) / (2 * sigma^2))
  
  fit_tmp <- list(sigma = sigma)
  p_upper <- rayleigh_1p_cdf(upper, fit_tmp)
  p_lower <- rayleigh_1p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

rayleigh_1p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 1) stop("Need at least 1 data point within truncation bounds")

  # Use non-truncated MLE as starting point
  init <- rayleigh_1p_fit(data)

  neg_log_likelihood <- function(params) {
    sigma <- exp(params[1])
    -rayleigh_1p_log_likelihood_truncated(data, sigma, lower, upper)
  }

  fit <- optim(
    par = c(log(init$sigma)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  sigma_hat <- exp(fit$par[1])

  log_lik <- rayleigh_1p_log_likelihood_truncated(data, sigma_hat, lower, upper)

  result <- list(
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_rayleigh_1p"
  )

  class(result) <- "truncated_rayleigh_1p"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
rayleigh_1p_pdf <- function(x, fit) {
  sigma <- fit$sigma
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x >= 0
  
  if (any(valid)) {
    z <- x[valid] / sigma
    res[valid] <- (z / sigma) * exp(-0.5 * z^2)
  }
  res
}

rayleigh_1p_logpdf <- function(x, fit) {
  sigma <- fit$sigma
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x >= 0
  
  if (any(valid)) {
    z <- x[valid] / sigma
    res[valid] <- log(z) - log(sigma) - 0.5 * z^2
  }
  res
}

rayleigh_1p_cdf <- function(x, fit) {
  sigma <- fit$sigma
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x >= 0
  
  if (any(valid)) {
    z <- x[valid] / sigma
    res[valid] <- 1 - exp(-0.5 * z^2)
  }
  res[x < 0] <- 0
  res
}

rayleigh_1p_logcdf <- function(x, fit) {
  log(rayleigh_1p_cdf(x, fit))
}

rayleigh_1p_sf <- function(x, fit) {
  sigma <- fit$sigma
  
  res <- rep(1, length(x))
  valid <- !is.na(x) & x >= 0
  
  if (any(valid)) {
    z <- x[valid] / sigma
    res[valid] <- exp(-0.5 * z^2)
  }
  res[x < 0] <- 1
  res
}

rayleigh_1p_logsf <- function(x, fit) {
  sigma <- fit$sigma
  
  res <- rep(0, length(x))
  valid <- !is.na(x) & x >= 0
  
  if (any(valid)) {
    z <- x[valid] / sigma
    res[valid] <- -0.5 * z^2
  }
  res
}

rayleigh_1p_quantile <- function(p, fit) {
  sigma <- fit$sigma
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  sigma * sqrt(-2 * log(1 - p))
}

rayleigh_1p_isf <- function(p, fit) {
  rayleigh_1p_quantile(1 - p, fit)
}

rayleigh_1p_rand <- function(n, fit) {
  u <- runif(n)
  rayleigh_1p_quantile(u, fit)
}

# -------------------------------
# Moments
# -------------------------------
rayleigh_1p_moment <- function(n, fit) {
  # E[X^n] = sigma^n * 2^(n/2) * Gamma(1 + n/2)
  fit$sigma^n * 2^(n/2) * gamma(1 + n/2)
}

rayleigh_1p_mean <- function(fit) {
  fit$sigma * sqrt(pi / 2)
}

rayleigh_1p_var <- function(fit) {
  fit$sigma^2 * (4 - pi) / 2
}

rayleigh_1p_std <- function(fit) {
  sqrt(rayleigh_1p_var(fit))
}

rayleigh_1p_skew <- function(fit) {
  2 * sqrt(pi) * (pi - 3) / (4 - pi)^(1.5)
}

rayleigh_1p_kurtosis <- function(fit) {
  -(6 * pi^2 - 24 * pi + 16) / (4 - pi)^2 + 3
}

rayleigh_1p_median <- function(fit) {
  fit$sigma * sqrt(log(4))
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
rayleigh_1p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  rayleigh_1p_quantile(c(alpha, 1 - alpha), fit)
}

rayleigh_1p_entropy <- function(fit) {
  1 + log(fit$sigma / sqrt(2)) + 0.57721566 / 2
}

rayleigh_1p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * rayleigh_1p_pdf(x, fit)
  }
  # Support is [0, Inf)
  upper <- rayleigh_1p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = 0, upper = upper)
  res$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.rayleigh_1p <- function(object, ...) {
  structure(object$log_likelihood, df = 1, nobs = object$n, class = "logLik")
}

logLik.truncated_rayleigh_1p <- function(object, ...) {
  structure(object$log_likelihood, df = 1, nobs = object$n, class = "logLik")
}
