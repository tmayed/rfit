# Lognormal Distribution Module (Improved + S3 compatible)

# -------------------------------
# Log-likelihood (pure)
# -------------------------------
lognormal_2p_log_likelihood <- function(data, mu, sigma) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) == 0 || sigma <= 0) {
    return(-Inf)
  }

  sum(dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE))
}


# -------------------------------
# Fit standard lognormal (MLE)
# -------------------------------
lognormal_2p_fit <- function(data) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) < 2) {
    stop("Need at least 2 valid positive data points")
  }

  log_data <- log(data)

  mu_hat <- mean(log_data)
  sigma_hat <- sqrt(mean((log_data - mu_hat)^2))  # MLE (divide by n)

  log_lik <- lognormal_2p_log_likelihood(data, mu_hat, sigma_hat)

  result <- list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "lognormal_2p"
  )

  class(result) <- "lognormal_2p"
  result
}


# -------------------------------
# Truncated lognormal (MLE via optim)
# -------------------------------
lognormal_2p_log_likelihood_truncated <- function(data, mu, sigma, lower, upper) {
  if (sigma <= 0 || lower < 0) return(-Inf)

  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)

  ll <- sum(dlnorm(data, meanlog = mu, sdlog = sigma, log = TRUE))

  p_upper <- plnorm(upper, meanlog = mu, sdlog = sigma)
  p_lower <- plnorm(lower, meanlog = mu, sdlog = sigma)

  diff <- p_upper - p_lower
  if (diff <= 0) return(-Inf)

  ll - length(data) * log(diff)
}


lognormal_2p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data)]

  if (lower < 0) {
    stop("Lower bound must be >= 0 for lognormal")
  }

  data <- data[data >= lower & data <= upper]

  if (length(data) < 2) {
    stop("Need at least 2 data points within truncation bounds")
  }

  log_data <- log(data)

  initial_mu <- mean(log_data)
  initial_sigma <- sd(log_data)

  neg_log_likelihood <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])

    if (sigma < 1e-8) return(Inf)

    -lognormal_2p_log_likelihood_truncated(data, mu, sigma, lower, upper)
  }

  fit <- optim(
    par = c(initial_mu, log(initial_sigma)),
    fn = neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  mu_hat <- fit$par[1]
  sigma_hat <- exp(fit$par[2])

  log_lik <- lognormal_2p_log_likelihood_truncated(
    data, mu_hat, sigma_hat, lower, upper
  )

  result <- list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_lognormal_2p"
  )

  class(result) <- "truncated_lognormal_2p"
  result
}


# -------------------------------
# S3: logLik (enables AIC/BIC)
# -------------------------------
logLik.lognormal_2p <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 2,
    nobs = object$n,
    class = "logLik"
  )
}

logLik.truncated_lognormal_2p <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 2,
    nobs = object$n,
    class = "logLik"
  )
}


# -------------------------------
# PDF
# -------------------------------
lognormal_2p_pdf <- function(x, fit) {
  out <- dlnorm(x, meanlog = fit$mu, sdlog = fit$sigma)
  out[x <= 0] <- 0
  out
}


# -------------------------------
# CDF
# -------------------------------
lognormal_2p_cdf <- function(x, fit) {
  ifelse(x <= 0, 0, plnorm(x, meanlog = fit$mu, sdlog = fit$sigma))
}


# -------------------------------
# Quantile
# -------------------------------
lognormal_2p_quantile <- function(p, fit) {
  if (any(p < 0 | p > 1)) {
    stop("Probabilities must be in [0,1]")
  }

  qlnorm(p, meanlog = fit$mu, sdlog = fit$sigma)
}


# -------------------------------
# Random generation
# -------------------------------
lognormal_2p_rand <- function(n, fit) {
  rlnorm(n, meanlog = fit$mu, sdlog = fit$sigma)
}


# -------------------------------
# Mean
# -------------------------------
lognormal_2p_mean <- function(fit) {
  exp(fit$mu + fit$sigma^2 / 2)
}


# -------------------------------
# Standard deviation
# -------------------------------
lognormal_2p_std <- function(fit) {
  sqrt((exp(fit$sigma^2) - 1) * exp(2 * fit$mu + fit$sigma^2))
}


# -------------------------------
# Survival Function
# -------------------------------
lognormal_2p_sf <- function(x, fit) {
  plnorm(x, meanlog = fit$mu, sdlog = fit$sigma, lower.tail = FALSE)
}


# -------------------------------
# Inverse Survival Function
# -------------------------------
lognormal_2p_isf <- function(p, fit) {
  qlnorm(p, meanlog = fit$mu, sdlog = fit$sigma, lower.tail = FALSE)
}


# -------------------------------
# Log PDF
# -------------------------------
lognormal_2p_logpdf <- function(x, fit) {
  dlnorm(x, meanlog = fit$mu, sdlog = fit$sigma, log = TRUE)
}


# -------------------------------
# Log CDF
# -------------------------------
lognormal_2p_logcdf <- function(x, fit) {
  plnorm(x, meanlog = fit$mu, sdlog = fit$sigma, log.p = TRUE)
}


# -------------------------------
# Log SF
# -------------------------------
lognormal_2p_logsf <- function(x, fit) {
  plnorm(x, meanlog = fit$mu, sdlog = fit$sigma, lower.tail = FALSE, log.p = TRUE)
}


# -------------------------------
# Variance
# -------------------------------
lognormal_2p_var <- function(fit) {
  (exp(fit$sigma^2) - 1) * exp(2 * fit$mu + fit$sigma^2)
}


# -------------------------------
# Moment
# -------------------------------
lognormal_2p_moment <- function(n, fit) {
  exp(n * fit$mu + (n^2 * fit$sigma^2) / 2)
}


# -------------------------------
# Skewness
# -------------------------------
lognormal_2p_skew <- function(fit) {
  (exp(fit$sigma^2) + 2) * sqrt(exp(fit$sigma^2) - 1)
}


# -------------------------------
# Kurtosis
# -------------------------------
lognormal_2p_kurtosis <- function(fit) {
  exp(4 * fit$sigma^2) + 2 * exp(3 * fit$sigma^2) + 3 * exp(2 * fit$sigma^2) - 3
}


# -------------------------------
# Median
# -------------------------------
lognormal_2p_median <- function(fit) {
  exp(fit$mu)
}


# -------------------------------
# Interval
# -------------------------------
lognormal_2p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  qlnorm(c(alpha, 1 - alpha), meanlog = fit$mu, sdlog = fit$sigma)
}


# -------------------------------
# Entropy
# -------------------------------
lognormal_2p_entropy <- function(fit) {
  log(fit$sigma * exp(fit$mu + 0.5) * sqrt(2 * pi))
}


# -------------------------------
# Expect
# -------------------------------
lognormal_2p_expect <- function(func, fit, ...) {
  # Numerical integration for arbitrary function
  integrand <- function(x) {
    func(x, ...) * lognormal_2p_pdf(x, fit)
  }
  
  # For lognormal, support is (0, Inf)
  # We use a large upper bound if Inf causes issues, or try integrate with Inf
  res <- integrate(integrand, lower = 0, upper = Inf)
  res$value
}
