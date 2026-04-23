# Gamma Distribution Module

# -------------------------------
# Log-likelihood
# -------------------------------
gamma_2p_log_likelihood <- function(data, shape, scale) {
  data <- data[!is.na(data) & data > 0]
  if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) {
    return(-Inf)
  }

  ll <- sum(dgamma(data, shape = shape, scale = scale, log = TRUE))

  if (!is.finite(ll)) return(-Inf)
  ll
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
gamma_2p_fit <- function(data) {
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 2) stop("Need at least 2 valid positive data points")

  # Method of moments for initial guesses
  m <- mean(data)
  v <- var(data)
  initial_shape <- m^2 / v
  initial_scale <- v / m

  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])
    scale <- exp(params[2])
    -gamma_2p_log_likelihood(data, shape, scale)
  }

  fit <- optim(
    par = c(log(initial_shape), log(initial_scale)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  shape_hat <- exp(fit$par[1])
  scale_hat <- exp(fit$par[2])

  log_lik <- gamma_2p_log_likelihood(data, shape_hat, scale_hat)

  result <- list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "gamma_2p"
  )

  class(result) <- "gamma_2p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
gamma_2p_log_likelihood_truncated <- function(data, shape, scale, lower, upper) {
  if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) {
    return(-Inf)
  }
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  ll <- gamma_2p_log_likelihood(data, shape, scale)
  if (!is.finite(ll)) return(-Inf)
  
  p_upper <- pgamma(upper, shape = shape, scale = scale)
  p_lower <- pgamma(lower, shape = shape, scale = scale)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll - length(data) * log(diff)
}

gamma_2p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 2) stop("Need at least 2 data points within truncation bounds")

  init <- gamma_2p_fit(data)

  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])
    scale <- exp(params[2])
    -gamma_2p_log_likelihood_truncated(data, shape, scale, lower, upper)
  }

  fit <- optim(
    par = c(log(init$shape), log(init$scale)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  shape_hat <- exp(fit$par[1])
  scale_hat <- exp(fit$par[2])

  log_lik <- gamma_2p_log_likelihood_truncated(data, shape_hat, scale_hat, lower, upper)

  result <- list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_gamma_2p"
  )

  class(result) <- "truncated_gamma_2p"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
gamma_2p_pdf <- function(x, fit) {
  dgamma(x, shape = fit$shape, scale = fit$scale)
}

gamma_2p_logpdf <- function(x, fit) {
  dgamma(x, shape = fit$shape, scale = fit$scale, log = TRUE)
}

gamma_2p_cdf <- function(x, fit) {
  pgamma(x, shape = fit$shape, scale = fit$scale)
}

gamma_2p_logcdf <- function(x, fit) {
  pgamma(x, shape = fit$shape, scale = fit$scale, log.p = TRUE)
}

gamma_2p_sf <- function(x, fit) {
  pgamma(x, shape = fit$shape, scale = fit$scale, lower.tail = FALSE)
}

gamma_2p_logsf <- function(x, fit) {
  pgamma(x, shape = fit$shape, scale = fit$scale, lower.tail = FALSE, log.p = TRUE)
}

gamma_2p_quantile <- function(p, fit) {
  qgamma(p, shape = fit$shape, scale = fit$scale)
}

gamma_2p_isf <- function(p, fit) {
  qgamma(p, shape = fit$shape, scale = fit$scale, lower.tail = FALSE)
}

gamma_2p_rand <- function(n, fit) {
  rgamma(n, shape = fit$shape, scale = fit$scale)
}

# -------------------------------
# Moments
# -------------------------------
gamma_2p_moment <- function(n, fit) {
  fit$scale^n * gamma(fit$shape + n) / gamma(fit$shape)
}

gamma_2p_mean <- function(fit) {
  fit$shape * fit$scale
}

gamma_2p_var <- function(fit) {
  fit$shape * fit$scale^2
}

gamma_2p_std <- function(fit) {
  sqrt(gamma_2p_var(fit))
}

gamma_2p_skew <- function(fit) {
  2 / sqrt(fit$shape)
}

gamma_2p_kurtosis <- function(fit) {
  6 / fit$shape + 3
}

gamma_2p_median <- function(fit) {
  gamma_2p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
gamma_2p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  gamma_2p_quantile(c(alpha, 1 - alpha), fit)
}

gamma_2p_entropy <- function(fit) {
  fit$shape + log(fit$scale) + lgamma(fit$shape) + (1 - fit$shape) * digamma(fit$shape)
}

gamma_2p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * gamma_2p_pdf(x, fit)
  }
  # Support is [0, Inf)
  # Use high quantile for upper bound of integration
  upper <- gamma_2p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = 0, upper = upper)
  res$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.gamma_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}

logLik.truncated_gamma_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}
