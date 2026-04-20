# Weibull Distribution Module

# -------------------------------
# Log-likelihood
# -------------------------------
weibull_log_likelihood <- function(data, shape, scale) {
  data <- data[!is.na(data) & data > 0]
  if (shape <= 0 || scale <= 0) return(-Inf)
  sum(dweibull(data, shape = shape, scale = scale, log = TRUE))
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
weibull_fit <- function(data) {
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 2) stop("Need at least 2 valid positive data points")

  # Heuristic initial guesses
  initial_shape <- 1
  initial_scale <- mean(data)

  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])
    scale <- exp(params[2])
    -weibull_log_likelihood(data, shape, scale)
  }

  fit <- optim(
    par = c(log(initial_shape), log(initial_scale)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  shape_hat <- exp(fit$par[1])
  scale_hat <- exp(fit$par[2])

  log_lik <- weibull_log_likelihood(data, shape_hat, scale_hat)

  result <- list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "weibull"
  )

  class(result) <- "weibull"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
weibull_log_likelihood_truncated <- function(data, shape, scale, lower, upper) {
  if (shape <= 0 || scale <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  ll <- sum(dweibull(data, shape = shape, scale = scale, log = TRUE))
  
  p_upper <- pweibull(upper, shape = shape, scale = scale)
  p_lower <- pweibull(lower, shape = shape, scale = scale)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll - length(data) * log(diff)
}

weibull_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 2) stop("Need at least 2 data points within truncation bounds")

  init <- weibull_fit(data)

  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])
    scale <- exp(params[2])
    -weibull_log_likelihood_truncated(data, shape, scale, lower, upper)
  }

  fit <- optim(
    par = c(log(init$shape), log(init$scale)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  shape_hat <- exp(fit$par[1])
  scale_hat <- exp(fit$par[2])

  log_lik <- weibull_log_likelihood_truncated(data, shape_hat, scale_hat, lower, upper)

  result <- list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_weibull"
  )

  class(result) <- "truncated_weibull"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
weibull_pdf <- function(x, fit) {
  dweibull(x, shape = fit$shape, scale = fit$scale)
}

weibull_logpdf <- function(x, fit) {
  dweibull(x, shape = fit$shape, scale = fit$scale, log = TRUE)
}

weibull_cdf <- function(x, fit) {
  pweibull(x, shape = fit$shape, scale = fit$scale)
}

weibull_logcdf <- function(x, fit) {
  pweibull(x, shape = fit$shape, scale = fit$scale, log.p = TRUE)
}

weibull_sf <- function(x, fit) {
  pweibull(x, shape = fit$shape, scale = fit$scale, lower.tail = FALSE)
}

weibull_logsf <- function(x, fit) {
  pweibull(x, shape = fit$shape, scale = fit$scale, lower.tail = FALSE, log.p = TRUE)
}

weibull_quantile <- function(p, fit) {
  qweibull(p, shape = fit$shape, scale = fit$scale)
}

weibull_isf <- function(p, fit) {
  qweibull(p, shape = fit$shape, scale = fit$scale, lower.tail = FALSE)
}

weibull_rand <- function(n, fit) {
  rweibull(n, shape = fit$shape, scale = fit$scale)
}

# -------------------------------
# Moments
# -------------------------------
weibull_moment <- function(n, fit) {
  fit$scale^n * gamma(1 + n/fit$shape)
}

weibull_mean <- function(fit) {
  weibull_moment(1, fit)
}

weibull_var <- function(fit) {
  weibull_moment(2, fit) - weibull_mean(fit)^2
}

weibull_std <- function(fit) {
  sqrt(weibull_var(fit))
}

weibull_skew <- function(fit) {
  k <- fit$shape
  g1 <- gamma(1 + 1/k)
  g2 <- gamma(1 + 2/k)
  g3 <- gamma(1 + 3/k)
  (g3 - 3*g1*g2 + 2*g1^3) / (g2 - g1^2)^(1.5)
}

weibull_kurtosis <- function(fit) {
  k <- fit$shape
  g1 <- gamma(1 + 1/k)
  g2 <- gamma(1 + 2/k)
  g3 <- gamma(1 + 3/k)
  g4 <- gamma(1 + 4/k)
  (g4 - 4*g1*g3 + 6*g1^2*g2 - 3*g1^4) / (g2 - g1^2)^2
}

weibull_median <- function(fit) {
  fit$scale * (log(2))^(1/fit$shape)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
weibull_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  weibull_quantile(c(alpha, 1 - alpha), fit)
}

weibull_entropy <- function(fit) {
  euler_gamma <- 0.5772156649015328606
  euler_gamma * (1 - 1/fit$shape) + log(fit$scale / fit$shape) + 1
}

weibull_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * weibull_pdf(x, fit)
  }
  # Support is [0, Inf)
  # Use high quantile for upper bound of integration
  upper <- weibull_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = 0, upper = upper)
  res$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.weibull <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}

logLik.truncated_weibull <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}
