# Johnson SB Distribution Module
# The Johnson SB distribution is a four-parameter family of probability distributions
# defined on a bounded interval [xi, xi + lambda].
# Z = gamma + delta * log((x - xi) / (xi + lambda - x)) ~ N(0, 1)

# -------------------------------
# Log-likelihood
# -------------------------------
johnsonsb_log_likelihood <- function(data, gamma, delta, xi, lambda) {
  data <- data[!is.na(data)]
  if (delta <= 0 || lambda <= 0) return(-Inf)
  
  # Support check
  if (any(data <= xi | data >= xi + lambda)) return(-Inf)
  
  fit <- list(gamma = gamma, delta = delta, xi = xi, lambda = lambda)
  log_dens <- johnsonsb_logpdf(data, fit)
  
  sum(log_dens)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
johnsonsb_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 4) stop("Need at least 4 valid data points")
  
  # Heuristic for xi and lambda (support)
  xi_hat <- min(data) - 0.01 * (max(data) - min(data))
  lambda_hat <- (max(data) - min(data)) + 0.02 * (max(data) - min(data))
  
  # Initial guesses for shape
  gamma_init <- 0
  delta_init <- 1
  
  neg_log_likelihood <- function(params) {
    gamma  <- params[1]
    delta  <- exp(params[2])
    xi     <- params[3]
    lambda <- exp(params[4])
    
    # Constraints
    if (any(data <= xi | data >= xi + lambda)) return(1e10)
    
    -johnsonsb_log_likelihood(data, gamma, delta, xi, lambda)
  }
  
  fit_opt <- optim(
    par = c(gamma_init, log(delta_init), xi_hat, log(lambda_hat)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead", # Nelder-Mead is often more robust for bounded distributions
    control = list(maxit = 5000)
  )
  
  gamma_hat  <- fit_opt$par[1]
  delta_hat  <- exp(fit_opt$par[2])
  xi_hat     <- fit_opt$par[3]
  lambda_hat <- exp(fit_opt$par[4])
  
  log_lik <- johnsonsb_log_likelihood(data, gamma_hat, delta_hat, xi_hat, lambda_hat)
  
  result <- list(
    gamma = gamma_hat,
    delta = delta_hat,
    xi = xi_hat,
    lambda = lambda_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "johnsonsb"
  )
  
  class(result) <- "johnsonsb"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
johnsonsb_fit_truncated <- function(data, lower = -Inf, upper = Inf) {
  # Johnson SB is already on a bounded interval [xi, xi + lambda].
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 4) stop("Need at least 4 data points within bounds")
  
  # Standard fit on the restricted data
  res <- johnsonsb_fit(data)
  res$distribution <- "truncated_johnsonsb"
  class(res) <- "truncated_johnsonsb"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
johnsonsb_pdf <- function(x, fit) {
  exp(johnsonsb_logpdf(x, fit))
}

johnsonsb_logpdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  out <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > xi & x < xi + lambda
  
  if (any(valid)) {
    u <- (x[valid] - xi) / lambda
    # log f(x) = log(delta) - log(lambda) - log(u) - log(1-u) - log(sqrt(2pi)) - 0.5 * (gamma + delta * log(u/(1-u)))^2
    z <- gamma + delta * log(u / (1 - u))
    out[valid] <- log(delta) - log(lambda) - log(u) - log(1 - u) - 0.5 * log(2 * pi) - 0.5 * z^2
  }
  out
}

johnsonsb_cdf <- function(x, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x > xi & x < xi + lambda
  
  if (any(valid)) {
    u <- (x[valid] - xi) / lambda
    z <- gamma + delta * log(u / (1 - u))
    res[valid] <- pnorm(z)
  }
  
  res[x >= xi + lambda] <- 1
  res[x <= xi] <- 0
  res
}

johnsonsb_logcdf <- function(x, fit) {
  log(johnsonsb_cdf(x, fit))
}

johnsonsb_sf <- function(x, fit) {
  1 - johnsonsb_cdf(x, fit)
}

johnsonsb_logsf <- function(x, fit) {
  log(johnsonsb_sf(x, fit))
}

johnsonsb_quantile <- function(p, fit) {
  gamma <- fit$gamma
  delta <- fit$delta
  xi <- fit$xi
  lambda <- fit$lambda
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  z_norm <- qnorm(p)
  # log(u/(1-u)) = (z_norm - gamma) / delta
  # u/(1-u) = exp((z_norm - gamma) / delta) = K
  # u = K / (1 + K) = 1 / (1 + 1/K)
  K_inv <- exp(-(z_norm - gamma) / delta)
  u <- 1 / (1 + K_inv)
  xi + lambda * u
}

johnsonsb_isf <- function(p, fit) {
  johnsonsb_quantile(1 - p, fit)
}

johnsonsb_rand <- function(n, fit) {
  u_norm <- runif(n)
  johnsonsb_quantile(u_norm, fit)
}

# -------------------------------
# Moments
# -------------------------------
johnsonsb_moment <- function(n, fit) {
  johnsonsb_expect(function(x) x^n, fit)
}

johnsonsb_mean <- function(fit) {
  johnsonsb_expect(function(x) x, fit)
}

johnsonsb_var <- function(fit) {
  m2 <- johnsonsb_moment(2, fit)
  m1 <- johnsonsb_mean(fit)
  m2 - m1^2
}

johnsonsb_std <- function(fit) {
  sqrt(johnsonsb_var(fit))
}

johnsonsb_skew <- function(fit) {
  set.seed(42)
  x <- johnsonsb_rand(100000, fit)
  m3 <- mean((x - mean(x))^3)
  m3 / (sd(x)^3)
}

johnsonsb_kurtosis <- function(fit) {
  set.seed(42)
  x <- johnsonsb_rand(100000, fit)
  m4 <- mean((x - mean(x))^4)
  m4 / (var(x)^2)
}

johnsonsb_median <- function(fit) {
  johnsonsb_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
johnsonsb_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  johnsonsb_quantile(c(alpha, 1 - alpha), fit)
}

johnsonsb_entropy <- function(fit) {
  integrand <- function(x) {
    p <- johnsonsb_pdf(x, fit)
    ifelse(p > 0, -p * log(p), 0)
  }
  # Use exact support for integration
  integrate(integrand, lower = fit$xi, upper = fit$xi + fit$lambda)$value
}

johnsonsb_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * johnsonsb_pdf(x, fit)
  }
  integrate(integrand, lower = fit$xi, upper = fit$xi + fit$lambda)$value
}

# -------------------------------
# S3: logLik
# -------------------------------
logLik.johnsonsb <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}

logLik.truncated_johnsonsb <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}
