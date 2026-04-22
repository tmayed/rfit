# Nakagami Distribution Module (2-parameter: m and omega)
# m: shape parameter (m >= 0.5)
# omega: spread parameter (omega > 0)
# PDF: f(x; m, omega) = [2 * m^m / (gamma(m) * omega^m)] * x^(2m-1) * exp(-m * x^2 / omega)

# -------------------------------
# Log-likelihood
# -------------------------------
nakagami_2p_log_likelihood <- function(data, m, omega) {
  data <- data[!is.na(data) & data > 0]
  if (m < 0.001 || omega <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  # LogPDF: log(2) + m*log(m) - lgamma(m) - m*log(omega) + (2m-1)*log(x) - m*x^2/omega
  sum(log(2) + m * log(m) - lgamma(m) - m * log(omega) + (2 * m - 1) * log(data) - (m * data^2) / omega)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
nakagami_2p_fit <- function(data) {
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 2) stop("Need at least 2 valid positive data points")

  # Heuristic initial guesses
  # omega is E[X^2]
  omega_init <- mean(data^2)
  # m can be estimated via (E[X^2])^2 / Var(X^2)
  m_init <- max(0.5, omega_init^2 / var(data^2))

  neg_log_likelihood <- function(params) {
    m <- exp(params[1])
    omega <- exp(params[2])
    
    # Nakagami usually requires m >= 0.5, but we allow slightly lower for flexibility
    if (m < 0.1) return(1e15)
    
    -nakagami_2p_log_likelihood(data, m, omega)
  }

  fit <- optim(
    par = c(log(m_init), log(omega_init)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  m_hat <- exp(fit$par[1])
  omega_hat <- exp(fit$par[2])

  log_lik <- nakagami_2p_log_likelihood(data, m_hat, omega_hat)

  result <- list(
    m = m_hat,
    omega = omega_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "nakagami_2p"
  )

  class(result) <- "nakagami_2p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
nakagami_2p_log_likelihood_truncated <- function(data, m, omega, lower, upper) {
  if (m <= 0 || omega <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  fit_tmp <- list(m = m, omega = omega)
  ll_base <- sum(nakagami_2p_logpdf(data, fit_tmp))
  
  p_upper <- nakagami_2p_cdf(upper, fit_tmp)
  p_lower <- nakagami_2p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

nakagami_2p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 2) stop("Need at least 2 data points within truncation bounds")

  init <- nakagami_2p_fit(data)

  neg_log_likelihood <- function(params) {
    m <- exp(params[1])
    omega <- exp(params[2])
    -nakagami_2p_log_likelihood_truncated(data, m, omega, lower, upper)
  }

  fit <- optim(
    par = c(log(init$m), log(init$omega)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  res <- list(
    m = exp(fit$par[1]),
    omega = exp(fit$par[2]),
    log_likelihood = -fit$value,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_nakagami_2p"
  )
  class(res) <- "truncated_nakagami_2p"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
nakagami_2p_pdf <- function(x, fit) {
  exp(nakagami_2p_logpdf(x, fit))
}

nakagami_2p_logpdf <- function(x, fit) {
  m <- fit$m
  omega <- fit$omega
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > 0
  
  if (any(valid)) {
    xv <- x[valid]
    # Use Gamma relationship: X^2 ~ Gamma(m, omega/m)
    # log f_X(x) = log(f_Y(x^2) * 2x)
    res[valid] <- dgamma(xv^2, shape = m, scale = omega / m, log = TRUE) + log(2 * xv)
  }
  res
}

nakagami_2p_cdf <- function(x, fit) {
  m <- fit$m
  omega <- fit$omega
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x > 0
  
  if (any(valid)) {
    # F_X(x) = P(X <= x) = P(X^2 <= x^2) = F_Y(x^2)
    res[valid] <- pgamma(x[valid]^2, shape = m, scale = omega / m)
  }
  res
}

nakagami_2p_logcdf <- function(x, fit) {
  log(nakagami_2p_cdf(x, fit))
}

nakagami_2p_sf <- function(x, fit) {
  1 - nakagami_2p_cdf(x, fit)
}

nakagami_2p_logsf <- function(x, fit) {
  log(nakagami_2p_sf(x, fit))
}

nakagami_2p_quantile <- function(p, fit) {
  m <- fit$m
  omega <- fit$omega
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  # x = sqrt(qgamma(p, shape=m, scale=omega/m))
  sqrt(qgamma(p, shape = m, scale = omega / m))
}

nakagami_2p_isf <- function(p, fit) {
  nakagami_2p_quantile(1 - p, fit)
}

nakagami_2p_rand <- function(n, fit) {
  sqrt(rgamma(n, shape = fit$m, scale = fit$omega / fit$m))
}

# -------------------------------
# Moments
# -------------------------------
nakagami_2p_moment <- function(n, fit) {
  m <- fit$m
  omega <- fit$omega
  # E[X^n] = [gamma(m + n/2) / gamma(m)] * (omega/m)^(n/2)
  (gamma(m + n / 2) / gamma(m)) * (omega / m)^(n / 2)
}

nakagami_2p_mean <- function(fit) {
  nakagami_2p_moment(1, fit)
}

nakagami_2p_var <- function(fit) {
  fit$omega - nakagami_2p_mean(fit)^2
}

nakagami_2p_std <- function(fit) {
  sqrt(nakagami_2p_var(fit))
}

nakagami_2p_skew <- function(fit) {
  m <- fit$m
  o <- fit$omega
  mean_x <- nakagami_2p_mean(fit)
  var_x  <- nakagami_2p_var(fit)
  mu3 <- nakagami_2p_moment(3, fit) - 3 * mean_x * var_x - mean_x^3
  mu3 / var_x^(1.5)
}

nakagami_2p_kurtosis <- function(fit) {
  # Kurtosis calculation is involved, let's use expectation
  m <- nakagami_2p_mean(fit)
  v <- nakagami_2p_var(fit)
  nakagami_2p_expect(function(x) ((x - m) / sqrt(v))^4, fit)
}

nakagami_2p_median <- function(fit) {
  nakagami_2p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
nakagami_2p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  nakagami_2p_quantile(c(alpha, 1 - alpha), fit)
}

nakagami_2p_entropy <- function(fit) {
  m <- fit$m
  omega <- fit$omega
  # Entropy(Nakagami) = log(gamma(m)) - (m-1/2)*psi(m) + m + log(omega/(2*m))
  lgamma(m) - (m - 0.5) * digamma(m) + m + log(omega / (2 * m))
}

nakagami_2p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * nakagami_2p_pdf(x, fit)
  }
  upper <- nakagami_2p_quantile(1 - 1e-10, fit)
  integrate(integrand, lower = 0, upper = upper)$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.nakagami_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}

logLik.truncated_nakagami_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}
