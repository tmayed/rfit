# Lévy Distribution Module (2-parameter: mu and sigma)
# PDF: f(x; mu, sigma) = sqrt(sigma/(2*pi)) * exp(-sigma/(2*(x-mu))) / (x-mu)^(3/2)
# CDF: F(x; mu, sigma) = erfc(sqrt(sigma/(2*(x-mu))))

# -------------------------------
# Log-likelihood
# -------------------------------
levy_2p_log_likelihood <- function(data, mu, sigma) {
  data <- data[!is.na(data) & data > mu]
  if (sigma <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  # LogPDF: 0.5*log(sigma/(2*pi)) - sigma/(2*(x-mu)) - 1.5*log(x-mu)
  sum(0.5 * log(sigma / (2 * pi)) - sigma / (2 * (data - mu)) - 1.5 * log(data - mu))
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
levy_2p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 2) stop("Need at least 2 valid data points")

  # Heuristic initial guesses
  # mu must be < min(data)
  mu_init <- min(data) - 0.01 * sd(data)
  
  # Closed form for sigma given mu: sigma = n / sum(1/(x-mu))
  get_sigma_mle <- function(mu, d) {
    diffs <- d - mu
    n <- length(d[d > mu])
    if (n == 0) return(1.0)
    n / sum(1 / diffs[diffs > 0])
  }
  
  sigma_init <- get_sigma_mle(mu_init, data)

  neg_log_likelihood <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])
    
    if (any(data <= mu)) return(1e15)
    
    -levy_2p_log_likelihood(data, mu, sigma)
  }

  # Nelder-Mead is usually safer for location bounds
  fit <- optim(
    par = c(mu_init, log(sigma_init)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 2000)
  )

  mu_hat <- fit$par[1]
  sigma_hat <- exp(fit$par[2])

  log_lik <- levy_2p_log_likelihood(data, mu_hat, sigma_hat)

  result <- list(
    mu = mu_hat,
    sigma = sigma_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "levy_2p"
  )

  class(result) <- "levy_2p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
levy_2p_log_likelihood_truncated <- function(data, mu, sigma, lower, upper) {
  if (sigma <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  if (any(data <= mu)) return(-Inf)
  
  fit_tmp <- list(mu = mu, sigma = sigma)
  ll_base <- sum(levy_2p_logpdf(data, fit_tmp))
  
  p_upper <- levy_2p_cdf(upper, fit_tmp)
  p_lower <- levy_2p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

levy_2p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 2) stop("Need at least 2 data points within truncation bounds")

  init <- levy_2p_fit(data)

  neg_log_likelihood <- function(params) {
    mu <- params[1]
    sigma <- exp(params[2])
    
    if (mu >= lower) return(1e15)
    
    -levy_2p_log_likelihood_truncated(data, mu, sigma, lower, upper)
  }

  fit <- optim(
    par = c(init$mu, log(init$sigma)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 2000)
  )

  res <- list(
    mu = fit$par[1],
    sigma = exp(fit$par[2]),
    log_likelihood = -fit$value,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_levy_2p"
  )
  class(res) <- "truncated_levy_2p"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
levy_2p_pdf <- function(x, fit) {
  exp(levy_2p_logpdf(x, fit))
}

levy_2p_logpdf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > mu
  
  if (any(valid)) {
    xv <- x[valid]
    res[valid] <- 0.5 * log(sigma / (2 * pi)) - sigma / (2 * (xv - mu)) - 1.5 * log(xv - mu)
  }
  res
}

levy_2p_cdf <- function(x, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x > mu
  
  if (any(valid)) {
    # F(x) = erfc(sqrt(sigma/(2(x-mu)))) = 2 * (1 - Phi(sqrt(sigma/(x-mu))))
    z <- sqrt(sigma / (x[valid] - mu))
    res[valid] <- 2 * pnorm(z, lower.tail = FALSE)
  }
  res[x <= mu] <- 0
  res
}

levy_2p_logcdf <- function(x, fit) {
  log(levy_2p_cdf(x, fit))
}

levy_2p_sf <- function(x, fit) {
  1 - levy_2p_cdf(x, fit)
}

levy_2p_logsf <- function(x, fit) {
  log(levy_2p_sf(x, fit))
}

levy_2p_quantile <- function(p, fit) {
  mu <- fit$mu
  sigma <- fit$sigma
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  # x = mu + sigma / (Phi^-1(1 - p/2))^2
  res <- mu + sigma / (qnorm(1 - p / 2)^2)
  res[p == 0] <- mu
  res[p == 1] <- Inf
  res
}

levy_2p_isf <- function(p, fit) {
  levy_2p_quantile(1 - p, fit)
}

levy_2p_rand <- function(n, fit) {
  u <- runif(n)
  levy_2p_quantile(u, fit)
}

# -------------------------------
# Moments (Lévy has no finite mean or variance)
# -------------------------------
levy_2p_moment <- function(n, fit) {
  if (n >= 0.5) return(Inf)
  levy_2p_expect(function(x) x^n, fit)
}

levy_2p_mean <- function(fit) {
  Inf
}

levy_2p_var <- function(fit) {
  Inf
}

levy_2p_std <- function(fit) {
  Inf
}

levy_2p_skew <- function(fit) {
  NA
}

levy_2p_kurtosis <- function(fit) {
  NA
}

levy_2p_median <- function(fit) {
  levy_2p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
levy_2p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  levy_2p_quantile(c(alpha, 1 - alpha), fit)
}

levy_2p_entropy <- function(fit) {
  # (1 + 3*euler_gamma)/2 + log(16*pi*sigma^2)/2
  # But we use mu-sigma convention. sigma is the scale c.
  0.5 * (1 + 3 * 0.57721566) + 0.5 * log(16 * pi * fit$sigma^2)
}

levy_2p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * levy_2p_pdf(x, fit)
  }
  # Heavy tail, use very high quantile
  upper <- levy_2p_quantile(1 - 1e-8, fit)
  integrate(integrand, lower = fit$mu, upper = upper)$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.levy_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}

logLik.truncated_levy_2p <- function(object, ...) {
  structure(object$log_likelihood, df = 2, nobs = object$n, class = "logLik")
}
