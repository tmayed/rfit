# Johnson SU Distribution Module
# The Johnson SU distribution is a four-parameter family of probability distributions
# transformed from the normal distribution.
# Z = gamma + delta * asinh((x - xi) / lambda) ~ N(0, 1)

# -------------------------------
# Log-likelihood
# -------------------------------
johnsonsu_4p_log_likelihood <- function(data, gamma, delta, xi, lambda) {
  data <- data[!is.na(data)]
  if (delta <= 0 || lambda <= 0) return(-Inf)
  
  fit <- list(gamma = gamma, delta = delta, xi = xi, lambda = lambda)
  log_dens <- johnsonsu_4p_logpdf(data, fit)
  
  sum(log_dens)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
johnsonsu_4p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 4) stop("Need at least 4 valid data points")
  
  # Initial guesses
  xi_init <- median(data)
  lambda_init <- sd(data)
  gamma_init <- 0
  delta_init <- 1
  
  neg_log_likelihood <- function(params) {
    gamma  <- params[1]
    delta  <- exp(params[2])
    xi     <- params[3]
    lambda <- exp(params[4])
    
    -johnsonsu_4p_log_likelihood(data, gamma, delta, xi, lambda)
  }
  
  fit_opt <- optim(
    par = c(gamma_init, log(delta_init), xi_init, log(lambda_init)),
    fn = neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 2000)
  )
  
  gamma_hat  <- fit_opt$par[1]
  delta_hat  <- exp(fit_opt$par[2])
  xi_hat     <- fit_opt$par[3]
  lambda_hat <- exp(fit_opt$par[4])
  
  log_lik <- johnsonsu_4p_log_likelihood(data, gamma_hat, delta_hat, xi_hat, lambda_hat)
  
  result <- list(
    gamma = gamma_hat,
    delta = delta_hat,
    xi = xi_hat,
    lambda = lambda_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "johnsonsu_4p"
  )
  
  class(result) <- "johnsonsu_4p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
johnsonsu_4p_log_likelihood_truncated <- function(data, gamma, delta, xi, lambda, lower, upper) {
  if (delta <= 0 || lambda <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  fit_tmp <- list(gamma = gamma, delta = delta, xi = xi, lambda = lambda)
  ll <- sum(johnsonsu_4p_logpdf(data, fit_tmp))
  
  p_upper <- johnsonsu_4p_cdf(upper, fit_tmp)
  p_lower <- johnsonsu_4p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll - length(data) * log(diff)
}

johnsonsu_4p_fit_truncated <- function(data, lower = -Inf, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 4) stop("Need at least 4 data points within truncation bounds")

  init <- johnsonsu_4p_fit(data)

  neg_log_likelihood <- function(params) {
    gamma  <- params[1]
    delta  <- exp(params[2])
    xi     <- params[3]
    lambda <- exp(params[4])
    
    -johnsonsu_4p_log_likelihood_truncated(data, gamma, delta, xi, lambda, lower, upper)
  }

  fit_opt <- optim(
    par = c(init$gamma, log(init$delta), init$xi, log(init$lambda)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  gamma_hat  <- fit_opt$par[1]
  delta_hat  <- exp(fit_opt$par[2])
  xi_hat     <- fit_opt$par[3]
  lambda_hat <- exp(fit_opt$par[4])

  log_lik <- johnsonsu_4p_log_likelihood_truncated(data, gamma_hat, delta_hat, xi_hat, lambda_hat, lower, upper)

  result <- list(
    gamma = gamma_hat,
    delta = delta_hat,
    xi = xi_hat,
    lambda = lambda_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit_opt$convergence,
    distribution = "truncated_johnsonsu_4p"
  )

  class(result) <- "truncated_johnsonsu_4p"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
johnsonsu_4p_pdf <- function(x, fit) {
  exp(johnsonsu_4p_logpdf(x, fit))
}

johnsonsu_4p_logpdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  z <- (x - xi) / lambda
  sq_z <- sqrt(1 + z^2)
  trans_z <- gamma + delta * asinh(z)
  
  log_delta <- log(delta)
  log_lambda <- log(lambda)
  log_sqrt_2pi <- 0.5 * log(2 * pi)
  
  log_f <- log_delta - log_lambda - log_sqrt_2pi - log(sq_z) - 0.5 * trans_z^2
  log_f[!is.finite(log_f)] <- -Inf
  log_f
}

johnsonsu_4p_cdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  z <- (x - xi) / lambda
  trans_z <- gamma + delta * asinh(z)
  pnorm(trans_z)
}

johnsonsu_4p_logcdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  z <- (x - xi) / lambda
  trans_z <- gamma + delta * asinh(z)
  pnorm(trans_z, log.p = TRUE)
}

johnsonsu_4p_sf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  z <- (x - xi) / lambda
  trans_z <- gamma + delta * asinh(z)
  pnorm(trans_z, lower.tail = FALSE)
}

johnsonsu_4p_logsf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  z <- (x - xi) / lambda
  trans_z <- gamma + delta * asinh(z)
  pnorm(trans_z, lower.tail = FALSE, log.p = TRUE)
}

johnsonsu_4p_quantile <- function(p, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  z_norm <- qnorm(p)
  fit$xi + fit$lambda * sinh((z_norm - fit$gamma) / fit$delta)
}

johnsonsu_4p_isf <- function(p, fit) {
  johnsonsu_4p_quantile(1 - p, fit)
}

johnsonsu_4p_rand <- function(n, fit) {
  z_norm <- rnorm(n)
  fit$xi + fit$lambda * sinh((z_norm - fit$gamma) / fit$delta)
}

# -------------------------------
# Moments
# -------------------------------
johnsonsu_4p_moment <- function(n, fit) {
  # Numerical integration for general non-central moments
  johnsonsu_4p_expect(function(x) x^n, fit)
}

johnsonsu_4p_mean <- function(fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  w <- exp(delta^(-2))
  xi - lambda * sqrt(w) * sinh(gamma / delta)
}

johnsonsu_4p_var <- function(fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  lambda <- fit$lambda
  
  w <- exp(delta^(-2))
  0.5 * lambda^2 * (w - 1) * (w * cosh(2 * gamma / delta) + 1)
}

johnsonsu_4p_std <- function(fit) {
  sqrt(johnsonsu_4p_var(fit))
}

johnsonsu_4p_skew <- function(fit) {
  # Analytical skewness exists but it is very complex.
  # Fallback to simulation for robustness.
  set.seed(42)
  x <- johnsonsu_4p_rand(100000, fit)
  m3 <- mean((x - mean(x))^3)
  m3 / (sd(x)^3)
}

johnsonsu_4p_kurtosis <- function(fit) {
  set.seed(42)
  x <- johnsonsu_4p_rand(100000, fit)
  m4 <- mean((x - mean(x))^4)
  m4 / (var(x)^2)
}

johnsonsu_4p_median <- function(fit) {
  johnsonsu_4p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
johnsonsu_4p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  johnsonsu_4p_quantile(c(alpha, 1 - alpha), fit)
}

johnsonsu_4p_entropy <- function(fit) {
  integrand <- function(x) {
    p <- johnsonsu_4p_pdf(x, fit)
    ifelse(p > 0, -p * log(p), 0)
  }
  # Determine integration bounds from quantiles
  lower <- johnsonsu_4p_quantile(1e-10, fit)
  upper <- johnsonsu_4p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}

johnsonsu_4p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * johnsonsu_4p_pdf(x, fit)
  }
  lower <- johnsonsu_4p_quantile(1e-10, fit)
  upper <- johnsonsu_4p_quantile(1 - 1e-10, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}

# -------------------------------
# S3: logLik
# -------------------------------
logLik.johnsonsu_4p <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}

logLik.truncated_johnsonsu_4p <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}
