# Birnbaum-Saunders (Fatigue-Life) Distribution Module (3-parameter: alpha, beta, mu)
# alpha: shape parameter (alpha > 0)
# beta: scale parameter (beta > 0)
# mu: location parameter (x > mu)
# CDF: F(x) = Phi( [sqrt((x-mu)/beta) - sqrt(beta/(x-mu))] / alpha )

# -------------------------------
# Log-likelihood
# -------------------------------
birnbaumsaunders_3p_log_likelihood <- function(data, alpha, beta, mu) {
  data <- data[!is.na(data) & data > mu]
  if (alpha <= 0 || beta <= 0) return(-Inf)
  
  n <- length(data)
  if (n == 0) return(-Inf)
  
  sum(birnbaumsaunders_3p_logpdf(data, list(alpha = alpha, beta = beta, mu = mu)))
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
birnbaumsaunders_3p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 3) stop("Need at least 3 valid data points")

  # Heuristic initial guesses
  mu_init <- min(data) - 0.01 * sd(data)
  d_shifted <- data - mu_init
  
  # For mu=0, beta_hat is sqrt(mean(x) * 1/mean(1/x))
  beta_init <- sqrt(mean(d_shifted) / mean(1/d_shifted))
  # alpha_hat is sqrt(mean(x/beta) + mean(beta/x) - 2)
  alpha_init <- sqrt(mean(d_shifted/beta_init) + mean(beta_init/d_shifted) - 2)

  neg_log_likelihood <- function(params) {
    alpha <- exp(params[1])
    beta <- exp(params[2])
    mu <- params[3]
    
    if (any(data <= mu)) return(1e15)
    
    val <- -birnbaumsaunders_3p_log_likelihood(data, alpha, beta, mu)
    if (!is.finite(val)) return(1e15)
    val
  }

  fit <- optim(
    par = c(log(alpha_init), log(beta_init), mu_init),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 3000)
  )

  res <- list(
    alpha = exp(fit$par[1]),
    beta = exp(fit$par[2]),
    mu = fit$par[3],
    log_likelihood = -fit$value,
    n = length(data),
    distribution = "birnbaumsaunders_3p"
  )

  class(res) <- "birnbaumsaunders_3p"
  res
}

# -------------------------------
# Truncated Fit
# -------------------------------
birnbaumsaunders_3p_log_likelihood_truncated <- function(data, alpha, beta, mu, lower, upper) {
  if (alpha <= 0 || beta <= 0) return(-Inf)
  
  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)
  if (any(data <= mu)) return(-Inf)
  
  fit_tmp <- list(alpha = alpha, beta = beta, mu = mu)
  ll_base <- sum(birnbaumsaunders_3p_logpdf(data, fit_tmp))
  
  p_upper <- birnbaumsaunders_3p_cdf(upper, fit_tmp)
  p_lower <- birnbaumsaunders_3p_cdf(lower, fit_tmp)
  
  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)
  
  ll_base - length(data) * log(diff)
}

birnbaumsaunders_3p_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 3) stop("Need at least 3 data points within truncation bounds")

  init <- birnbaumsaunders_3p_fit(data)

  neg_log_likelihood <- function(params) {
    alpha <- exp(params[1])
    beta <- exp(params[2])
    mu <- params[3]
    
    if (mu >= lower) return(1e15)
    
    -birnbaumsaunders_3p_log_likelihood_truncated(data, alpha, beta, mu, lower, upper)
  }

  fit <- optim(
    par = c(log(init$alpha), log(init$beta), init$mu),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 3000)
  )

  res <- list(
    alpha = exp(fit$par[1]),
    beta = exp(fit$par[2]),
    mu = fit$par[3],
    log_likelihood = -fit$value,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_birnbaumsaunders_3p"
  )
  class(res) <- "truncated_birnbaumsaunders_3p"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
birnbaumsaunders_3p_pdf <- function(x, fit) {
  exp(birnbaumsaunders_3p_logpdf(x, fit))
}

birnbaumsaunders_3p_logpdf <- function(x, fit) {
  alpha <- fit$alpha
  beta <- fit$beta
  mu <- fit$mu
  
  res <- rep(-Inf, length(x))
  valid <- !is.na(x) & x > mu
  
  if (any(valid)) {
    xv <- x[valid] - mu
    z <- (sqrt(xv/beta) - sqrt(beta/xv)) / alpha
    
    # f(x) = [ (sqrt(x/beta) + sqrt(beta/x)) / (2*alpha*x) ] * phi(z)
    term1 <- log(sqrt(xv/beta) + sqrt(beta/xv))
    term2 <- log(2 * alpha * xv)
    log_phi <- dnorm(z, log = TRUE)
    
    res[valid] <- term1 - term2 + log_phi
  }
  res
}

birnbaumsaunders_3p_cdf <- function(x, fit) {
  alpha <- fit$alpha
  beta <- fit$beta
  mu <- fit$mu
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x > mu
  
  if (any(valid)) {
    xv <- x[valid] - mu
    z <- (sqrt(xv/beta) - sqrt(beta/xv)) / alpha
    res[valid] <- pnorm(z)
  }
  res[x <= mu] <- 0
  res
}

birnbaumsaunders_3p_logcdf <- function(x, fit) {
  log(birnbaumsaunders_3p_cdf(x, fit))
}

birnbaumsaunders_3p_sf <- function(x, fit) {
  1 - birnbaumsaunders_3p_cdf(x, fit)
}

birnbaumsaunders_3p_logsf <- function(x, fit) {
  log(birnbaumsaunders_3p_sf(x, fit))
}

birnbaumsaunders_3p_quantile <- function(p, fit) {
  alpha <- fit$alpha
  beta <- fit$beta
  mu <- fit$mu
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  # z = Phi^-1(p)
  # x = mu + beta * [1 + 0.5*(alpha*z)^2 + alpha*z * sqrt(1 + 0.25*(alpha*z)^2)]
  z <- qnorm(p)
  az <- alpha * z
  res <- mu + beta * (1 + 0.5 * az^2 + az * sqrt(1 + 0.25 * az^2))
  res[p == 0] <- mu
  res[p == 1] <- Inf
  res
}

birnbaumsaunders_3p_isf <- function(p, fit) {
  birnbaumsaunders_3p_quantile(1 - p, fit)
}

birnbaumsaunders_3p_rand <- function(n, fit) {
  u <- runif(n)
  birnbaumsaunders_3p_quantile(u, fit)
}

# -------------------------------
# Moments
# -------------------------------
birnbaumsaunders_3p_mean <- function(fit) {
  # E[X] = mu + beta * (1 + alpha^2 / 2)
  fit$mu + fit$beta * (1 + fit$alpha^2 / 2)
}

birnbaumsaunders_3p_var <- function(fit) {
  # Var(X) = (alpha*beta)^2 * (1 + 5/4 * alpha^2)
  (fit$alpha * fit$beta)^2 * (1 + 1.25 * fit$alpha^2)
}

birnbaumsaunders_3p_std <- function(fit) {
  sqrt(birnbaumsaunders_3p_var(fit))
}

birnbaumsaunders_3p_moment <- function(n, fit) {
  birnbaumsaunders_3p_expect(function(x) x^n, fit)
}

birnbaumsaunders_3p_skew <- function(fit) {
  a <- fit$alpha
  # Skewness = [4*a*(1 + 11/8 * a^2)] / (1 + 5/4 * a^2)^1.5
  (4 * a * (1 + 1.375 * a^2)) / (1 + 1.25 * a^2)^1.5
}

birnbaumsaunders_3p_kurtosis <- function(fit) {
  a <- fit$alpha
  # Kurtosis = 3 + [3*a^2 * (1 + 5/2 * a^2 + 7/16 * a^4)] / (1 + 5/4 * a^2)^2
  # This formula is complex, let's use expectation for Kurtosis to be safe or re-verify.
  # Let's use expectation.
  m <- birnbaumsaunders_3p_mean(fit)
  v <- birnbaumsaunders_3p_var(fit)
  birnbaumsaunders_3p_expect(function(x) ((x - m) / sqrt(v))^4, fit)
}

birnbaumsaunders_3p_median <- function(fit) {
  # Median is mu + beta
  fit$mu + fit$beta
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
birnbaumsaunders_3p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  birnbaumsaunders_3p_quantile(c(alpha, 1 - alpha), fit)
}

birnbaumsaunders_3p_entropy <- function(fit) {
  birnbaumsaunders_3p_expect(function(x) -birnbaumsaunders_3p_logpdf(x, fit), fit)
}

birnbaumsaunders_3p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * birnbaumsaunders_3p_pdf(x, fit)
  }
  upper <- birnbaumsaunders_3p_quantile(1 - 1e-10, fit)
  integrate(integrand, lower = fit$mu, upper = upper)$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.birnbaumsaunders_3p <- function(object, ...) {
  structure(object$log_likelihood, df = 3, nobs = object$n, class = "logLik")
}

logLik.truncated_birnbaumsaunders_3p <- function(object, ...) {
  structure(object$log_likelihood, df = 3, nobs = object$n, class = "logLik")
}
