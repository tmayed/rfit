# Johnson SL Distribution Module
# The Johnson SL distribution is a three-parameter family of probability distributions
# that is equivalent to the Lognormal distribution with a location shift.
# Z = gamma + delta * log(x - xi) ~ N(0, 1)

# -------------------------------
# Log-likelihood
# -------------------------------
johnsonsl_3p_log_likelihood <- function(data, gamma, delta, xi) {
  data <- data[!is.na(data)]
  if (delta <= 0) return(-Inf)
  if (any(data <= xi)) return(-Inf)
  
  fit <- list(gamma = gamma, delta = delta, xi = xi)
  log_dens <- johnsonsl_3p_logpdf(data, fit)
  
  sum(log_dens)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
johnsonsl_3p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 3) stop("Need at least 3 valid data points")
  
  # Initial guesses
  xi_init <- min(data) - 0.1 * sd(data)
  log_data_shifted <- log(data - xi_init)
  delta_init <- 1 / sd(log_data_shifted)
  gamma_init <- -mean(log_data_shifted) * delta_init
  
  neg_log_likelihood <- function(params) {
    gamma  <- params[1]
    delta  <- exp(params[2])
    xi     <- params[3]
    
    if (any(data <= xi)) return(1e10)
    
    -johnsonsl_3p_log_likelihood(data, gamma, delta, xi)
  }
  
  fit_opt <- optim(
    par = c(gamma_init, log(delta_init), xi_init),
    fn = neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 2000)
  )
  
  gamma_hat <- fit_opt$par[1]
  delta_hat <- exp(fit_opt$par[2])
  xi_hat    <- fit_opt$par[3]
  
  log_lik <- johnsonsl_3p_log_likelihood(data, gamma_hat, delta_hat, xi_hat)
  
  result <- list(
    gamma = gamma_hat,
    delta = delta_hat,
    xi = xi_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "johnsonsl_3p"
  )
  
  class(result) <- "johnsonsl_3p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
johnsonsl_3p_log_likelihood_truncated <- function(data, gamma, delta, xi, lower, upper) {
  if (delta <= 0) return(-Inf)
  if (any(data <= xi)) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  fit_tmp <- list(gamma = gamma, delta = delta, xi = xi)
  ll <- sum(johnsonsl_3p_logpdf(data, fit_tmp))
  
  p_upper <- johnsonsl_3p_cdf(upper, fit_tmp)
  p_lower <- johnsonsl_3p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll - length(data) * log(diff)
}

johnsonsl_3p_fit_truncated <- function(data, lower = -Inf, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 3) stop("Need at least 3 data points within truncation bounds")

  init <- johnsonsl_3p_fit(data)

  neg_log_likelihood <- function(params) {
    gamma <- params[1]
    delta <- exp(params[2])
    xi    <- params[3]
    
    -johnsonsl_3p_log_likelihood_truncated(data, gamma, delta, xi, lower, upper)
  }

  fit_opt <- optim(
    par = c(init$gamma, log(init$delta), init$xi),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  gamma_hat <- fit_opt$par[1]
  delta_hat <- exp(fit_opt$par[2])
  xi_hat    <- fit_opt$par[3]

  log_lik <- johnsonsl_3p_log_likelihood_truncated(data, gamma_hat, delta_hat, xi_hat, lower, upper)

  result <- list(
    gamma = gamma_hat,
    delta = delta_hat,
    xi = xi_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit_opt$convergence,
    distribution = "truncated_johnsonsl_3p"
  )

  class(result) <- "truncated_johnsonsl_3p"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
johnsonsl_3p_pdf <- function(x, fit) {
  exp(johnsonsl_3p_logpdf(x, fit))
}

johnsonsl_3p_logpdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  
  out <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > xi
  
  if (any(valid)) {
    z <- log(x[valid] - xi)
    trans_z <- gamma + delta * z
    out[valid] <- log(delta) - 0.5 * log(2 * pi) - z - 0.5 * trans_z^2
  }
  out
}

johnsonsl_3p_cdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x > xi
  
  if (any(valid)) {
    z <- log(x[valid] - xi)
    trans_z <- gamma + delta * z
    res[valid] <- pnorm(trans_z)
  }
  
  res[x <= xi] <- 0
  res
}

johnsonsl_3p_logcdf <- function(x, fit) {
  log(johnsonsl_3p_cdf(x, fit))
}

johnsonsl_3p_sf <- function(x, fit) {
  1 - johnsonsl_3p_cdf(x, fit)
}

johnsonsl_3p_logsf <- function(x, fit) {
  log(johnsonsl_3p_sf(x, fit))
}

johnsonsl_3p_quantile <- function(p, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  z_norm <- qnorm(p)
  # log(x - xi) = (z_norm - gamma) / delta
  xi + exp((z_norm - gamma) / delta)
}

johnsonsl_3p_isf <- function(p, fit) {
  johnsonsl_3p_quantile(1 - p, fit)
}

johnsonsl_3p_rand <- function(n, fit) {
  z_norm <- rnorm(n)
  fit$xi + exp((z_norm - fit$gamma) / fit$delta)
}

# -------------------------------
# Moments
# -------------------------------
johnsonsl_3p_moment <- function(n, fit) {
  # For Lognormal (SL), moments exist analytically:
  # E[(X-xi)^n] = exp(n * mu + 0.5 * n^2 * sigma^2)
  # where mu = -gamma/delta and sigma = 1/delta
  mu <- -fit$gamma / fit$delta
  sigma <- 1 / fit$delta
  
  # Non-central moment E[X^n] = E[( (X-xi) + xi )^n]
  # Using simulation for simplicity and correctness of E[X^n]
  set.seed(42)
  mean(johnsonsl_3p_rand(100000, fit)^n)
}

johnsonsl_3p_mean <- function(fit) {
  mu <- -fit$gamma / fit$delta
  sigma <- 1 / fit$delta
  fit$xi + exp(mu + 0.5 * sigma^2)
}

johnsonsl_3p_var <- function(fit) {
  mu <- -fit$gamma / fit$delta
  sigma <- 1 / fit$delta
  (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)
}

johnsonsl_3p_std <- function(fit) {
  sqrt(johnsonsl_3p_var(fit))
}

johnsonsl_3p_skew <- function(fit) {
  sigma <- 1 / fit$delta
  (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
}

johnsonsl_3p_kurtosis <- function(fit) {
  sigma <- 1 / fit$delta
  exp(4 * sigma^2) + 2 * exp(3 * sigma^2) + 3 * exp(2 * sigma^2) - 3
}

johnsonsl_3p_median <- function(fit) {
  johnsonsl_3p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
johnsonsl_3p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  johnsonsl_3p_quantile(c(alpha, 1 - alpha), fit)
}

johnsonsl_3p_entropy <- function(fit) {
  mu <- -fit$gamma / fit$delta
  sigma <- 1 / fit$delta
  # Entropy of Lognormal is mu + 0.5 + log(sigma * sqrt(2*pi))
  mu + 0.5 + log(sigma * sqrt(2 * pi))
}

johnsonsl_3p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * johnsonsl_3p_pdf(x, fit)
  }
  lower <- fit$xi
  upper <- johnsonsl_3p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}

# -------------------------------
# S3: logLik
# -------------------------------
logLik.johnsonsl_3p <- function(object, ...) {
  structure(object$log_likelihood, df = 3, nobs = object$n, class = "logLik")
}

logLik.truncated_johnsonsl_3p <- function(object, ...) {
  structure(object$log_likelihood, df = 3, nobs = object$n, class = "logLik")
}
