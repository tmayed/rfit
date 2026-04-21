# Generalized Beta of the Second Kind (GB2) Distribution Module

# -------------------------------
# Log-likelihood (pure)
# -------------------------------
gb2_log_likelihood <- function(data, a, b, p, q) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) == 0 || a <= 0 || b <= 0 || p <= 0 || q <= 0) {
    return(-Inf)
  }

  # log(f(x)) = log(a) + (a*p - 1)*log(x) - a*p*log(b) - lbeta(p, q) - (p + q)*log(1 + (x/b)^a)
  ll <- sum(log(a) + (a * p - 1) * log(data) - (a * p) * log(b) -
              lbeta(p, q) - (p + q) * log(1 + (data / b)^a))

  return(ll)
}


# -------------------------------
# Fit GB2 (MLE via optim)
# -------------------------------
gb2_fit <- function(data) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) < 5) {
    stop("Need at least 5 valid positive data points for GB2 estimation")
  }

  # Initial guesses
  # b is often close to the median
  initial_b <- median(data)
  # a, p, q can start at 1
  initial_a <- 1
  initial_p <- 1
  initial_q <- 1

  neg_log_likelihood <- function(params) {
    a <- exp(params[1])
    b <- exp(params[2])
    p <- exp(params[3])
    q <- exp(params[4])

    # Avoid numerical issues
    if (any(c(a, b, p, q) < 1e-8)) return(Inf)

    -gb2_log_likelihood(data, a, b, p, q)
  }

  fit <- optim(
    par = log(c(initial_a, initial_b, initial_p, initial_q)),
    fn = neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  a_hat <- exp(fit$par[1])
  b_hat <- exp(fit$par[2])
  p_hat <- exp(fit$par[3])
  q_hat <- exp(fit$par[4])

  log_lik <- gb2_log_likelihood(data, a_hat, b_hat, p_hat, q_hat)

  result <- list(
    a = a_hat,
    b = b_hat,
    p = p_hat,
    q = q_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "gb2",
    convergence = fit$convergence
  )

  class(result) <- "gb2"
  result
}


# -------------------------------
# Truncated GB2 (MLE via optim)
# -------------------------------
gb2_log_likelihood_truncated <- function(data, a, b, p, q, lower, upper) {
  if (a <= 0 || b <= 0 || p <= 0 || q <= 0 || lower < 0) return(-Inf)

  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)

  # Full log-likelihood
  ll <- sum(log(a) + (a * p - 1) * log(data) - (a * p) * log(b) -
              lbeta(p, q) - (p + q) * log(1 + (data / b)^a))

  # Normalization
  # F(x) = pbeta( (x/b)^a / (1 + (x/b)^a), p, q)
  
  f_cdf <- function(x) {
    if (x <= 0) return(0)
    if (x == Inf) return(1)
    
    # Calculate y = (x/b)^a / (1 + (x/b)^a)
    # y = exp(a * (log(x) - log(b))) / (1 + exp(a * (log(x) - log(b))))
    # For very large (x/b)^a, y approaches 1
    log_xb_a <- a * (log(x) - log(b))
    
    if (log_xb_a > 100) { # log_xb_a > 100 => xb_a is massive
      y <- 1
    } else {
      xb_a <- exp(log_xb_a)
      y <- xb_a / (1 + xb_a)
    }
    
    pbeta(y, p, q)
  }

  p_upper <- f_cdf(upper)
  p_lower <- f_cdf(lower)

  diff <- p_upper - p_lower
  if (diff <= 1e-10) return(-Inf)

  ll - length(data) * log(diff)
}


gb2_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data)]
  data <- data[data >= lower & data <= upper]

  if (length(data) < 5) {
    stop("Need at least 5 data points within truncation bounds")
  }

  # Initial guesses
  initial_b <- median(data)
  initial_a <- 1
  initial_p <- 1
  initial_q <- 1

  neg_log_likelihood <- function(params) {
    a <- exp(params[1])
    b <- exp(params[2])
    p <- exp(params[3])
    q <- exp(params[4])

    if (any(c(a, b, p, q) < 1e-8)) return(Inf)

    -gb2_log_likelihood_truncated(data, a, b, p, q, lower, upper)
  }

  fit <- optim(
    par = log(c(initial_a, initial_b, initial_p, initial_q)),
    fn = neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 1000)
  )

  a_hat <- exp(fit$par[1])
  b_hat <- exp(fit$par[2])
  p_hat <- exp(fit$par[3])
  q_hat <- exp(fit$par[4])

  log_lik <- gb2_log_likelihood_truncated(
    data, a_hat, b_hat, p_hat, q_hat, lower, upper
  )

  result <- list(
    a = a_hat,
    b = b_hat,
    p = p_hat,
    q = q_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_gb2"
  )

  class(result) <- "truncated_gb2"
  result
}


# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
gb2_pdf <- function(x, fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  out <- numeric(length(x))
  valid <- x > 0 & !is.na(x)

  # Use log-domain for numerical stability
  # log(f(x)) = log(a) + (a*p - 1)*log(x) - (a*p)*log(b) - lbeta(p, q) - (p + q)*log(1 + (x/b)^a)
  log_pdf <- log(a) + (a * p - 1) * log(x[valid]) - (a * p) * log(b) - 
             lbeta(p, q) - (p + q) * log(1 + (x[valid] / b)^a)
  
  out[valid] <- exp(log_pdf)
  out
}


gb2_cdf <- function(x, fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  res <- numeric(length(x))
  valid <- !is.na(x) & x > 0
  
  # y = (x / b)^a / (1 + (x / b)^a)
  # For stability: log_xb_a = a * (log(x) - log(b))
  log_xb_a <- a * (log(x[valid]) - log(b))
  
  xb_a <- exp(log_xb_a)
  y <- xb_a / (1 + xb_a)
  # Handle Inf/large values
  y[log_xb_a > 100] <- 1
  
  res[valid] <- pbeta(y, p, q)
  res[x <= 0 & !is.na(x)] <- 0
  res
}


gb2_quantile <- function(prob, fit) {
  if (any(prob < 0 | prob > 1)) {
    stop("Probabilities must be in [0,1]")
  }

  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  y <- qbeta(prob, p, q)
  # x = b * (y / (1 - y))^(1 / a)
  
  res <- numeric(length(prob))
  res[prob == 1] <- Inf
  res[prob == 0] <- 0
  
  valid <- prob > 0 & prob < 1
  res[valid] <- b * (y[valid] / (1 - y[valid]))^(1 / a)
  res
}


gb2_rand <- function(n, fit) {
  u <- runif(n)
  gb2_quantile(u, fit)
}


# -------------------------------
# Survival / Inverse Survival / Logs
# -------------------------------
gb2_sf <- function(x, fit) {
  1 - gb2_cdf(x, fit)
}

gb2_isf <- function(p, fit) {
  gb2_quantile(1 - p, fit)
}

gb2_logpdf <- function(x, fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  out <- rep(-Inf, length(x))
  valid <- x > 0 & !is.na(x)

  out[valid] <- log(a) + (a * p - 1) * log(x[valid]) - (a * p) * log(b) - 
                lbeta(p, q) - (p + q) * log(1 + (x[valid] / b)^a)
  out
}

gb2_logcdf <- function(x, fit) {
  log(gb2_cdf(x, fit))
}

gb2_logsf <- function(x, fit) {
  log(gb2_sf(x, fit))
}


# -------------------------------
# Moments
# -------------------------------
gb2_mean <- function(fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  if (a * q <= 1) return(Inf)

  exp(log(b) + lbeta(p + 1 / a, q - 1 / a) - lbeta(p, q))
}

gb2_var <- function(fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  if (a * q <= 2) return(Inf)

  log_mean_sq <- 2 * log(b) + lbeta(p + 2 / a, q - 2 / a) - lbeta(p, q)
  mean_sq <- exp(log_mean_sq)
  m <- gb2_mean(fit)
  max(0, mean_sq - m^2)
}

gb2_std <- function(fit) {
  sqrt(gb2_var(fit))
}

gb2_moment <- function(n, fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q

  if (a * q <= n) return(Inf)
  
  # E[X^n] = b^n * B(p + n/a, q - n/a) / B(p, q)
  exp(n * log(b) + lbeta(p + n / a, q - n / a) - lbeta(p, q))
}

gb2_skew <- function(fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q
  
  if (a * q <= 3) return(NA)
  
  m1 <- gb2_moment(1, fit)
  m2 <- gb2_moment(2, fit)
  m3 <- gb2_moment(3, fit)
  sigma <- sqrt(m2 - m1^2)
  
  (m3 - 3*m1*sigma^2 - m1^3) / sigma^3
}

gb2_kurtosis <- function(fit) {
  a <- fit$a
  b <- fit$b
  p <- fit$p
  q <- fit$q
  
  if (a * q <= 4) return(NA)
  
  m1 <- gb2_moment(1, fit)
  m2 <- gb2_moment(2, fit)
  m3 <- gb2_moment(3, fit)
  m4 <- gb2_moment(4, fit)
  sigma2 <- m2 - m1^2
  
  # Central 4th moment: E[(X-mu)^4] = E[X^4] - 4*mu*E[X^3] + 6*mu^2*E[X^2] - 3*mu^4
  mu <- m1
  mu4 <- m4 - 4*mu*m3 + 6*mu^2*m2 - 3*mu^4
  
  mu4 / sigma2^2
}

gb2_median <- function(fit) {
  gb2_quantile(0.5, fit)
}

gb2_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  gb2_quantile(c(alpha, 1 - alpha), fit)
}


# -------------------------------
# Entropy & Expect
# -------------------------------
gb2_entropy <- function(fit) {
  # Entropy for GB2 has a complex analytical form
  # Using numerical integration for reliability
  -gb2_expect(function(x) gb2_logpdf(x, fit), fit)
}

gb2_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * gb2_pdf(x, fit)
  }
  
  # GB2 support is (0, Inf)
  res <- integrate(integrand, lower = 0, upper = Inf)
  res$value
}


# -------------------------------
# S3: logLik (enables AIC/BIC)
# -------------------------------
logLik.gb2 <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 4,
    nobs = object$n,
    class = "logLik"
  )
}

logLik.truncated_gb2 <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 4,
    nobs = object$n,
    class = "logLik"
  )
}
