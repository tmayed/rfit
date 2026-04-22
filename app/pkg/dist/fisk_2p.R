# Fisk (Log-logistic) Distribution Module

# -------------------------------
# Log-likelihood
# -------------------------------
fisk_2p_log_likelihood <- function(data, scale, shape) {
  data <- data[!is.na(data) & data > 0]
  if (scale <= 0 || shape <= 0) return(-Inf)
  
  # Log-logistic is logistic on log-scale
  # Y = log(X) ~ Logistic(location = log(scale), scale = 1/shape)
  # f_X(x) = f_Y(log(x)) * |d/dx log(x)| = f_Y(log(x)) / x
  # log(f_X(x)) = log(f_Y(log(x))) - log(x)
  
  log_probs <- dlogis(log(data), location = log(scale), scale = 1/shape, log = TRUE) - log(data)
  sum(log_probs)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
fisk_2p_fit <- function(data) {
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 2) stop("Need at least 2 valid positive data points")

  # Heuristic initial guesses using relationship with Logistic distribution
  log_data <- log(data)
  mu_init <- mean(log_data)
  s_init <- sd(log_data) * sqrt(3) / pi
  
  initial_scale <- exp(mu_init)
  initial_shape <- 1 / max(s_init, 1e-3)

  neg_log_likelihood <- function(params) {
    scale <- exp(params[1])
    shape <- exp(params[2])
    -fisk_2p_log_likelihood(data, scale, shape)
  }

  fit <- optim(
    par = c(log(initial_scale), log(initial_shape)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  scale_hat <- exp(fit$par[1])
  shape_hat <- exp(fit$par[2])

  log_lik <- fisk_2p_log_likelihood(data, scale_hat, shape_hat)

  result <- list(
    scale = scale_hat,
    shape = shape_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "fisk_2p"
  )

  class(result) <- "fisk_2p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
fisk_2p_log_likelihood_truncated <- function(data, scale, shape, lower, upper) {
  if (scale <= 0 || shape <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  # Untruncated log-likelihood
  ll <- sum(dlogis(log(data), location = log(scale), scale = 1/shape, log = TRUE) - log(data))
  
  # Truncation adjustment
  p_upper <- plogis(log(upper), location = log(scale), scale = 1/shape)
  p_lower <- plogis(log(lower), location = log(scale), scale = 1/shape)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll - length(data) * log(diff)
}

fisk_2p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 2) stop("Need at least 2 data points within truncation bounds")

  # Use non-truncated fit as starting point
  init <- fisk_2p_fit(data)

  neg_log_likelihood <- function(params) {
    scale <- exp(params[1])
    shape <- exp(params[2])
    -fisk_2p_log_likelihood_truncated(data, scale, shape, lower, upper)
  }

  fit <- optim(
    par = c(log(init$scale), log(init$shape)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  scale_hat <- exp(fit$par[1])
  shape_hat <- exp(fit$par[2])

  log_lik <- fisk_2p_log_likelihood_truncated(data, scale_hat, shape_hat, lower, upper)

  result <- list(
    scale = scale_hat,
    shape = shape_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_fisk_2p"
  )

  class(result) <- "truncated_fisk_2p"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
fisk_2p_pdf <- function(x, fit) {
  # Avoid log(0)
  res <- numeric(length(x))
  pos <- x > 0
  if (any(pos)) {
    res[pos] <- dlogis(log(x[pos]), location = log(fit$scale), scale = 1/fit$shape) / x[pos]
  }
  res
}

fisk_2p_logpdf <- function(x, fit) {
  res <- rep(-Inf, length(x))
  pos <- x > 0
  if (any(pos)) {
    res[pos] <- dlogis(log(x[pos]), location = log(fit$scale), scale = 1/fit$shape, log = TRUE) - log(x[pos])
  }
  res
}

fisk_2p_cdf <- function(x, fit) {
  res <- numeric(length(x))
  pos <- x > 0
  if (any(pos)) {
    res[pos] <- plogis(log(x[pos]), location = log(fit$scale), scale = 1/fit$shape)
  }
  res
}

fisk_2p_logcdf <- function(x, fit) {
  res <- rep(-Inf, length(x))
  pos <- x > 0
  if (any(pos)) {
    res[pos] <- plogis(log(x[pos]), location = log(fit$scale), scale = 1/fit$shape, log.p = TRUE)
  }
  res
}

fisk_2p_sf <- function(x, fit) {
  res <- ones <- rep(1, length(x))
  pos <- x > 0
  if (any(pos)) {
    res[pos] <- plogis(log(x[pos]), location = log(fit$scale), scale = 1/fit$shape, lower.tail = FALSE)
  }
  res
}

fisk_2p_logsf <- function(x, fit) {
  res <- rep(0, length(x))
  pos <- x > 0
  if (any(pos)) {
    res[pos] <- plogis(log(x[pos]), location = log(fit$scale), scale = 1/fit$shape, lower.tail = FALSE, log.p = TRUE)
  }
  res
}

fisk_2p_quantile <- function(p, fit) {
  exp(qlogis(p, location = log(fit$scale), scale = 1/fit$shape))
}

fisk_2p_isf <- function(p, fit) {
  exp(qlogis(p, location = log(fit$scale), scale = 1/fit$shape, lower.tail = FALSE))
}

fisk_2p_rand <- function(n, fit) {
  exp(rlogis(n, location = log(fit$scale), scale = 1/fit$shape))
}

# -------------------------------
# Moments
# -------------------------------
fisk_2p_moment <- function(n, fit) {
  # E[X^n] = scale^n * (n*pi/shape) / sin(n*pi/shape)
  # Valid only for n < shape
  if (n >= fit$shape) return(NA)
  
  b <- pi / fit$shape
  fit$scale^n * (n * b) / sin(n * b)
}

fisk_2p_mean <- function(fit) {
  fisk_2p_moment(1, fit)
}

fisk_2p_var <- function(fit) {
  if (fit$shape <= 2) return(NA)
  fisk_2p_moment(2, fit) - fisk_2p_mean(fit)^2
}

fisk_2p_std <- function(fit) {
  sqrt(fisk_2p_var(fit))
}

fisk_2p_skew <- function(fit) {
  if (fit$shape <= 3) return(NA)
  
  # Standardized 3rd moment formula for Fisk
  # Using relationship with Logistic distribution via moments of log-logistic
  # Alternatively, use general formula:
  # skew = [E(X^3) - 3*mu*sigma^2 - mu^3] / sigma^3
  
  m1 <- fisk_2p_moment(1, fit)
  m2 <- fisk_2p_moment(2, fit)
  m3 <- fisk_2p_moment(3, fit)
  
  var_val <- m2 - m1^2
  (m3 - 3*m1*m2 + 2*m1^3) / (var_val^1.5)
}

fisk_2p_kurtosis <- function(fit) {
  if (fit$shape <= 4) return(NA)
  
  m1 <- fisk_2p_moment(1, fit)
  m2 <- fisk_2p_moment(2, fit)
  m3 <- fisk_2p_moment(3, fit)
  m4 <- fisk_2p_moment(4, fit)
  
  var_val <- m2 - m1^2
  (m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4) / (var_val^2)
}

fisk_2p_median <- function(fit) {
  fit$scale
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
fisk_2p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  fisk_2p_quantile(c(alpha, 1 - alpha), fit)
}

fisk_2p_entropy <- function(fit) {
  # Entropy = ln(scale) - ln(shape) + 2
  log(fit$scale) - log(fit$shape) + 2
}

fisk_2p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * fisk_2p_pdf(x, fit)
  }
  # Support is (0, Inf)
  # Use high quantile for upper bound of integration if necessary
  upper <- fisk_2p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = 0, upper = upper)
  res$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.fisk_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}

logLik.truncated_fisk_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}
