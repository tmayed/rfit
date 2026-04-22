# Rayleigh Distribution Module (2-parameter: location and scale)

# -------------------------------
# Log-likelihood
# -------------------------------
rayleigh_2p_log_likelihood <- function(data, mu, sigma) {
  data <- data[!is.na(data) & data >= mu]
  if (sigma <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  # PDF: (x-mu)/sigma^2 * exp(-(x-mu)^2 / (2*sigma^2))
  # LogPDF: log(x-mu) - 2*log(sigma) - (x-mu)^2 / (2*sigma^2)
  sum(log(data - mu) - 2 * log(sigma) - (data - mu)^2 / (2 * sigma^2))
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
rayleigh_2p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 2) stop("Need at least 2 valid data points")

  # Heuristic initial guesses
  mu_init <- min(data) - 0.01 * sd(data)
  sigma_init <- sqrt(mean((data - mu_init)^2) / 2)

  neg_log_likelihood <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])
    
    if (any(data < mu)) return(1e15)
    
    -rayleigh_2p_log_likelihood(data, mu, sigma)
  }

  fit <- optim(
    par = c(mu_init, log(sigma_init)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 2000)
  )

  mu_hat <- fit$par[1]
  sigma_hat <- exp(fit$par[2])

  log_lik <- rayleigh_2p_log_likelihood(data, mu_hat, sigma_hat)

  result <- list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "rayleigh_2p"
  )

  class(result) <- "rayleigh_2p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
rayleigh_2p_log_likelihood_truncated <- function(data, mu, sigma, lower, upper) {
  if (sigma <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  # Base log-likelihood for Rayleigh on (mu, Inf)
  # For points in data, we need them >= mu
  if (any(data < mu)) return(-Inf)
  
  ll_base <- sum(log(data - mu) - 2 * log(sigma) - (data - mu)^2 / (2 * sigma^2))
  
  p_upper <- rayleigh_2p_cdf(upper, list(mu = mu, sigma = sigma))
  p_lower <- rayleigh_2p_cdf(lower, list(mu = mu, sigma = sigma))
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

rayleigh_2p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 2) stop("Need at least 2 data points within truncation bounds")

  init <- rayleigh_2p_fit(data)

  neg_log_likelihood <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])
    
    if (mu > lower) return(1e15)
    
    -rayleigh_2p_log_likelihood_truncated(data, mu, sigma, lower, upper)
  }

  fit <- optim(
    par = c(init$mu, log(init$sigma)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 2000)
  )

  mu_hat <- fit$par[1]
  sigma_hat <- exp(fit$par[2])

  log_lik <- rayleigh_2p_log_likelihood_truncated(data, mu_hat, sigma_hat, lower, upper)

  result <- list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_rayleigh_2p"
  )

  class(result) <- "truncated_rayleigh_2p"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
rayleigh_2p_pdf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    z <- (x[valid] - mu) / sigma
    res[valid] <- (z / sigma) * exp(-0.5 * z^2)
  }
  res
}

rayleigh_2p_logpdf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    z <- (x[valid] - mu) / sigma
    res[valid] <- log(z) - log(sigma) - 0.5 * z^2
  }
  res
}

rayleigh_2p_cdf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    z <- (x[valid] - mu) / sigma
    res[valid] <- 1 - exp(-0.5 * z^2)
  }
  res[x < mu] <- 0
  res
}

rayleigh_2p_logcdf <- function(x, fit) {
  log(rayleigh_2p_cdf(x, fit))
}

rayleigh_2p_sf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- rep(1, length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    z <- (x[valid] - mu) / sigma
    res[valid] <- exp(-0.5 * z^2)
  }
  res[x < mu] <- 1
  res
}

rayleigh_2p_logsf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- rep(0, length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    z <- (x[valid] - mu) / sigma
    res[valid] <- -0.5 * z^2
  }
  res
}

rayleigh_2p_quantile <- function(p, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  mu + sigma * sqrt(-2 * log(1 - p))
}

rayleigh_2p_isf <- function(p, fit) {
  rayleigh_2p_quantile(1 - p, fit)
}

rayleigh_2p_rand <- function(n, fit) {
  u <- runif(n)
  rayleigh_2p_quantile(u, fit)
}

# -------------------------------
# Moments
# -------------------------------
rayleigh_2p_moment <- function(n, fit) {
  # This is complex for non-central moments with location shift.
  # Let's use numerical expectation.
  rayleigh_2p_expect(function(x) x^n, fit)
}

rayleigh_2p_mean <- function(fit) {
  fit$mu + fit$sigma * sqrt(pi / 2)
}

rayleigh_2p_var <- function(fit) {
  fit$sigma^2 * (4 - pi) / 2
}

rayleigh_2p_std <- function(fit) {
  sqrt(rayleigh_2p_var(fit))
}

rayleigh_2p_skew <- function(fit) {
  2 * sqrt(pi) * (pi - 3) / (4 - pi)^(1.5)
}

rayleigh_2p_kurtosis <- function(fit) {
  -(6 * pi^2 - 24 * pi + 16) / (4 - pi)^2 + 3
}

rayleigh_2p_median <- function(fit) {
  fit$mu + fit$sigma * sqrt(log(4))
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
rayleigh_2p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  rayleigh_2p_quantile(c(alpha, 1 - alpha), fit)
}

rayleigh_2p_entropy <- function(fit) {
  1 + log(fit$sigma / sqrt(2)) + 0.57721566 / 2
}

rayleigh_2p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * rayleigh_2p_pdf(x, fit)
  }
  # Support is [mu, Inf)
  upper <- rayleigh_2p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = fit$mu, upper = upper)
  res$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.rayleigh_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}

logLik.truncated_rayleigh_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}
