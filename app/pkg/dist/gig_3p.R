# Generalized Inverse Gaussian (GIG) Distribution Module (3-parameter: lambda, chi, psi)
# PDF: f(x) = [(psi/chi)^(lambda/2) / (2 * K_lambda(sqrt(chi*psi)))] * x^(lambda-1) * exp(-(chi/x + psi*x)/2)
# where K_lambda is the modified Bessel function of the second kind.

# -------------------------------
# Log-likelihood
# -------------------------------
gig_3p_log_likelihood <- function(data, lambda, chi, psi) {
  data <- data[!is.na(data) & data > 0]
  if (chi <= 0 || psi <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  sum(gig_3p_logpdf(data, list(lambda = lambda, chi = chi, psi = psi)))
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
gig_3p_fit <- function(data) {
  data <- data[!is.na(data) & data > 0]
  if (length(data) < 3) stop("Need at least 3 valid positive data points")

  # Heuristic initial guesses
  # lambda = 0 is a common starting point (Inverse Gaussian variant)
  lambda_init <- 0.5
  # chi and psi can be estimated from moments
  # For simplicity, we use values that make the peak near the data median
  m <- median(data)
  chi_init <- m
  psi_init <- 1/m

  neg_log_likelihood <- function(params) {
    lambda <- params[1]
    chi <- exp(params[2])
    psi <- exp(params[3])
    
    val <- -gig_3p_log_likelihood(data, lambda, chi, psi)
    if (!is.finite(val)) return(1e15)
    val
  }

  fit <- optim(
    par = c(lambda_init, log(chi_init), log(psi_init)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  res <- list(
    lambda = fit$par[1],
    chi = exp(fit$par[2]),
    psi = exp(fit$par[3]),
    log_likelihood = -fit$value,
    n = length(data),
    distribution = "gig_3p"
  )

  class(res) <- "gig_3p"
  res
}

# -------------------------------
# Truncated Fit
# -------------------------------
gig_3p_log_likelihood_truncated <- function(data, lambda, chi, psi, lower, upper) {
  if (chi <= 0 || psi <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  
  fit_tmp <- list(lambda = lambda, chi = chi, psi = psi)
  ll_base <- sum(gig_3p_logpdf(data, fit_tmp))
  
  p_upper <- gig_3p_cdf(upper, fit_tmp)
  p_lower <- gig_3p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

gig_3p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 3) stop("Need at least 3 data points within truncation bounds")

  init <- gig_3p_fit(data)

  neg_log_likelihood <- function(params) {
    lambda <- params[1]
    chi <- exp(params[2])
    psi <- exp(params[3])
    -gig_3p_log_likelihood_truncated(data, lambda, chi, psi, lower, upper)
  }

  fit <- optim(
    par = c(init$lambda, log(init$chi), log(init$psi)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )

  res <- list(
    lambda = fit$par[1],
    chi = exp(fit$par[2]),
    psi = exp(fit$par[3]),
    log_likelihood = -fit$value,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_gig_3p"
  )
  class(res) <- "truncated_gig_3p"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
gig_3p_pdf <- function(x, fit) {
  exp(gig_3p_logpdf(x, fit))
}

gig_3p_logpdf <- function(x, fit) {
  lambda <- fit$lambda
  chi <- fit$chi
  psi <- fit$psi
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > 0
  
  if (any(valid)) {
    xv <- x[valid]
    omega <- sqrt(chi * psi)
    
    # Standard GIG density formula
    # Constant = (psi/chi)^(lambda/2) / (2 * besselK(omega, lambda))
    log_const <- (lambda / 2) * (log(psi) - log(chi)) - log(2) - log(besselK(omega, lambda))
    
    res[valid] <- log_const + (lambda - 1) * log(xv) - (chi / xv + psi * xv) / 2
  }
  res
}

gig_3p_cdf <- function(x, fit) {
  res <- numeric(length(x))
  valid <- !is.na(x) & x > 0
  
  if (any(valid)) {
    # No closed form CDF, use numerical integration
    for (i in which(valid)) {
      res[i] <- tryCatch({
        integrate(function(t) gig_3p_pdf(t, fit), lower = 0, upper = x[i])$value
      }, error = function(e) 1.0) # Fallback to 1.0 if it blows up
    }
  }
  pmax(0, pmin(1, res))
}

gig_3p_logcdf <- function(x, fit) {
  log(gig_3p_cdf(x, fit))
}

gig_3p_sf <- function(x, fit) {
  1 - gig_3p_cdf(x, fit)
}

gig_3p_logsf <- function(x, fit) {
  log(gig_3p_sf(x, fit))
}

gig_3p_quantile <- function(p, fit) {
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  res <- numeric(length(p))
  for (i in seq_along(p)) {
    if (p[i] == 0) { res[i] <- 0; next }
    if (p[i] == 1) { res[i] <- Inf; next }
    
    target <- p[i]
    f <- function(x) gig_3p_cdf(x, fit) - target
    
    # Improved initial bounds using mean/mode properties
    # GIG can be very heavy or very light tailed. 
    # Let's use a wide search with extension allowed.
    low <- 1e-10
    high <- 10
    
    # Expand high bound if necessary
    # Check if target is high
    if (target > 0.99) {
       while(gig_3p_cdf(high, fit) < target && high < 1e12) high <- high * 10
    }
    
    res[i] <- uniroot(f, lower = 0, upper = high, extendInt = "yes", tol = 1e-8)$root
  }
  res
}

gig_3p_isf <- function(p, fit) {
  gig_3p_quantile(1 - p, fit)
}

gig_3p_rand <- function(n, fit) {
  # Probability integral transform
  u <- runif(n)
  gig_3p_quantile(u, fit)
}

# -------------------------------
# Moments
# -------------------------------
gig_3p_moment <- function(n, fit) {
  lambda <- fit$lambda
  chi <- fit$chi
  psi <- fit$psi
  omega <- sqrt(chi * psi)
  
  # E[X^n] = (chi/psi)^(n/2) * K_{lambda+n}(omega) / K_lambda(omega)
  (chi / psi)^(n / 2) * besselK(omega, lambda + n) / besselK(omega, lambda)
}

gig_3p_mean <- function(fit) {
  gig_3p_moment(1, fit)
}

gig_3p_var <- function(fit) {
  gig_3p_moment(2, fit) - gig_3p_mean(fit)^2
}

gig_3p_std <- function(fit) {
  sqrt(gig_3p_var(fit))
}

gig_3p_skew <- function(fit) {
  m <- gig_3p_mean(fit)
  v <- gig_3p_var(fit)
  (gig_3p_moment(3, fit) - 3 * m * v - m^3) / v^(1.5)
}

gig_3p_kurtosis <- function(fit) {
  m <- gig_3p_mean(fit)
  v <- gig_3p_var(fit)
  (gig_3p_moment(4, fit) - 4 * m * gig_3p_moment(3, fit) + 6 * m^2 * gig_3p_moment(2, fit) - 3 * m^4) / v^2
}

gig_3p_median <- function(fit) {
  gig_3p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
gig_3p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  gig_3p_quantile(c(alpha, 1 - alpha), fit)
}

gig_3p_entropy <- function(fit) {
  gig_3p_expect(function(x) -log(gig_3p_pdf(x, fit)), fit)
}

gig_3p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * gig_3p_pdf(x, fit)
  }
  upper <- gig_3p_quantile(1 - 1e-8, fit)
  integrate(integrand, lower = 0, upper = upper)$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.gig_3p <- function(object, ...) {
  structure(object$log_likelihood, df = 3, nobs = object$n, class = "logLik")
}

logLik.truncated_gig_3p <- function(object, ...) {
  structure(object$log_likelihood, df = 3, nobs = object$n, class = "logLik")
}
