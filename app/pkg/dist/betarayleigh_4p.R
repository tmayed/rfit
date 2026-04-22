# Beta-Rayleigh Distribution Module (4-parameter: a, b, mu, sigma)
# The Beta-Rayleigh distribution is derived from the Beta-G class.
# F(x) = I_G(x)(a, b), where G(x) is the Rayleigh CDF.

# -------------------------------
# Log-likelihood
# -------------------------------
betarayleigh_4p_log_likelihood <- function(data, a, b, mu, sigma) {
  data <- data[!is.na(data) & data >= mu]
  if (a <= 0 || b <= 0 || sigma <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  sum(betarayleigh_4p_logpdf(data, list(a = a, b = b, mu = mu, sigma = sigma)))
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
betarayleigh_4p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 4) stop("Need at least 4 valid data points")

  # Heuristic initial guesses
  mu_init <- min(data) - 0.01 * sd(data)
  # Start with a=1, b=1 (standard Rayleigh)
  a_init <- 1.0
  b_init <- 1.0
  sigma_init <- sqrt(mean((data - mu_init)^2) / 2)

  neg_log_likelihood <- function(params) {
    a <- exp(params[1])
    b <- exp(params[2])
    mu <- params[3]
    sigma <- exp(params[4])
    
    if (any(data < mu)) return(1e15)
    
    val <- -betarayleigh_4p_log_likelihood(data, a, b, mu, sigma)
    if (!is.finite(val)) return(1e15)
    val
  }

  fit <- optim(
    par = c(log(a_init), log(b_init), mu_init, log(sigma_init)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 3000)
  )

  a_hat <- exp(fit$par[1])
  b_hat <- exp(fit$par[2])
  mu_hat <- fit$par[3]
  sigma_hat <- exp(fit$par[4])

  log_lik <- betarayleigh_4p_log_likelihood(data, a_hat, b_hat, mu_hat, sigma_hat)

  result <- list(
    a = a_hat,
    b = b_hat,
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "betarayleigh_4p"
  )

  class(result) <- "betarayleigh_4p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
betarayleigh_4p_log_likelihood_truncated <- function(data, a, b, mu, sigma, lower, upper) {
  if (a <= 0 || b <= 0 || sigma <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  if (any(data < mu)) return(-Inf)
  
  fit_tmp <- list(a = a, b = b, mu = mu, sigma = sigma)
  ll_base <- sum(betarayleigh_4p_logpdf(data, fit_tmp))
  
  p_upper <- betarayleigh_4p_cdf(upper, fit_tmp)
  p_lower <- betarayleigh_4p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

betarayleigh_4p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 4) stop("Need at least 4 data points within truncation bounds")

  init <- betarayleigh_4p_fit(data)

  neg_log_likelihood <- function(params) {
    a <- exp(params[1])
    b <- exp(params[2])
    mu <- params[3]
    sigma <- exp(params[4])
    
    if (mu > lower) return(1e15)
    
    -betarayleigh_4p_log_likelihood_truncated(data, a, b, mu, sigma, lower, upper)
  }

  fit <- optim(
    par = c(log(init$a), log(init$b), init$mu, log(init$sigma)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 3000)
  )

  res <- list(
    a = exp(fit$par[1]),
    b = exp(fit$par[2]),
    mu = fit$par[3],
    sigma = exp(fit$par[4]),
    log_likelihood = -fit$value,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_betarayleigh_4p"
  )
  class(res) <- "truncated_betarayleigh_4p"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
betarayleigh_4p_pdf <- function(x, fit) {
  exp(betarayleigh_4p_logpdf(x, fit))
}

betarayleigh_4p_logpdf <- function(x, fit) {
  a <- fit$a
  b <- fit$b
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    xv <- x[valid]
    # G(x) = Rayleigh CDF
    z <- (xv - mu) / sigma
    G <- 1 - exp(-0.5 * z^2)
    G <- pmax(1e-15, pmin(1 - 1e-15, G)) # Numerical stability
    
    # g(x) = Rayleigh PDF
    log_g <- log(z) - log(sigma) - 0.5 * z^2
    
    # Beta-G PDF: 1/B(a,b) * G^(a-1) * (1-G)^(b-1) * g
    res[valid] <- (a - 1) * log(G) + (b - 1) * log(1 - G) + log_g - lbeta(a, b)
  }
  res
}

betarayleigh_4p_cdf <- function(x, fit) {
  a <- fit$a
  b <- fit$b
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x >= mu
  
  if (any(valid)) {
    z <- (x[valid] - mu) / sigma
    G <- 1 - exp(-0.5 * z^2)
    res[valid] <- pbeta(G, a, b)
  }
  res[x < mu] <- 0
  res
}

betarayleigh_4p_logcdf <- function(x, fit) {
  log(betarayleigh_4p_cdf(x, fit))
}

betarayleigh_4p_sf <- function(x, fit) {
  1 - betarayleigh_4p_cdf(x, fit)
}

betarayleigh_4p_logsf <- function(x, fit) {
  log(betarayleigh_4p_sf(x, fit))
}

betarayleigh_4p_quantile <- function(p, fit) {
  a <- fit$a
  b <- fit$b
  mu <- fit$mu
  sigma <- fit$sigma
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  # G^-1(qbeta(p, a, b))
  G_inv_p <- qbeta(p, a, b)
  mu + sigma * sqrt(-2 * log(1 - G_inv_p))
}

betarayleigh_4p_isf <- function(p, fit) {
  betarayleigh_4p_quantile(1 - p, fit)
}

betarayleigh_4p_rand <- function(n, fit) {
  # G^-1(rbeta(n, a, b))
  u <- rbeta(n, fit$a, fit$b)
  fit$mu + fit$sigma * sqrt(-2 * log(1 - u))
}

# -------------------------------
# Moments
# -------------------------------
betarayleigh_4p_moment <- function(n, fit) {
  betarayleigh_4p_expect(function(x) x^n, fit)
}

betarayleigh_4p_mean <- function(fit) {
  betarayleigh_4p_expect(function(x) x, fit)
}

betarayleigh_4p_var <- function(fit) {
  betarayleigh_4p_expect(function(x) x^2, fit) - betarayleigh_4p_mean(fit)^2
}

betarayleigh_4p_std <- function(fit) {
  sqrt(betarayleigh_4p_var(fit))
}

betarayleigh_4p_skew <- function(fit) {
  m <- betarayleigh_4p_mean(fit)
  s <- betarayleigh_4p_std(fit)
  betarayleigh_4p_expect(function(x) ((x - m) / s)^3, fit)
}

betarayleigh_4p_kurtosis <- function(fit) {
  m <- betarayleigh_4p_mean(fit)
  s <- betarayleigh_4p_std(fit)
  betarayleigh_4p_expect(function(x) ((x - m) / s)^4, fit)
}

betarayleigh_4p_median <- function(fit) {
  betarayleigh_4p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
betarayleigh_4p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  betarayleigh_4p_quantile(c(alpha, 1 - alpha), fit)
}

betarayleigh_4p_entropy <- function(fit) {
  integrand <- function(x) {
    p <- betarayleigh_4p_pdf(x, fit)
    ifelse(p > 0, -p * log(p), 0)
  }
  upper <- betarayleigh_4p_quantile(1 - 1e-10, fit)
  integrate(integrand, lower = fit$mu, upper = upper)$value
}

betarayleigh_4p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * betarayleigh_4p_pdf(x, fit)
  }
  upper <- betarayleigh_4p_quantile(1 - 1e-10, fit)
  integrate(integrand, lower = fit$mu, upper = upper)$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.betarayleigh_4p <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}

logLik.truncated_betarayleigh_4p <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}
