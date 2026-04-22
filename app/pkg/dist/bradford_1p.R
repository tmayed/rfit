# Bradford Distribution Module
# The Bradford distribution is a probability distribution defined on the interval [lower, upper].
# It is often used to model the distribution of scientific journals among different fields.

# -------------------------------
# Log-likelihood
# -------------------------------
bradford_1p_log_likelihood <- function(data, shape, lower, upper) {
  data <- data[!is.na(data)]
  data <- data[data >= lower & data <= upper]
  
  if (length(data) == 0 || shape <= 0 || lower >= upper) {
    return(-Inf)
  }
  
  fit <- list(shape = shape, lower = lower, upper = upper)
  log_dens <- bradford_1p_logpdf(data, fit)
  
  sum(log_dens)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
bradford_1p_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 1) stop("Need at least 1 valid data point")
  
  # For Bradford, the support [lower, upper] are often known or estimated by min/max
  lower_hat <- min(data)
  upper_hat <- max(data)
  
  # If data is constant, we can't fit
  if (lower_hat == upper_hat) {
    upper_hat <- lower_hat + 1
  }
  
  # Initial guess for shape
  shape_init <- 1.0
  
  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])
    # We fix lower and upper as min/max for standard MLE of support
    -bradford_1p_log_likelihood(data, shape, lower_hat, upper_hat)
  }
  
  fit_opt <- optim(
    par = c(log(shape_init)),
    fn = neg_log_likelihood,
    method = "BFGS"
  )
  
  shape_hat <- exp(fit_opt$par[1])
  
  log_lik <- bradford_1p_log_likelihood(data, shape_hat, lower_hat, upper_hat)
  
  result <- list(
    shape = shape_hat,
    lower = lower_hat,
    upper = upper_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "bradford_1p"
  )
  
  class(result) <- "bradford_1p"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
bradford_1p_fit_truncated <- function(data, lower_bound = -Inf, upper_bound = Inf) {
  # Bradford is already on a bounded interval [lower, upper].
  # Truncating it further just restricts that interval.
  data <- data[!is.na(data) & data >= lower_bound & data <= upper_bound]
  if (length(data) < 1) stop("Need data within truncation bounds")
  
  # Standard fit on the truncated data
  res <- bradford_1p_fit(data)
  res$distribution <- "truncated_bradford_1p"
  class(res) <- "truncated_bradford_1p"
  res
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
bradford_1p_pdf <- function(x, fit) {
  exp(bradford_1p_logpdf(x, fit))
}

bradford_1p_logpdf <- function(x, fit) {
  c <- fit$shape
  a <- fit$lower
  b <- fit$upper
  L <- b - a
  
  out <- rep(-Inf, length(x))
  valid <- !is.na(x) & x >= a & x <= b
  
  if (any(valid)) {
    # f(x) = c / (log(1+c) * (L + c*(x-a)))
    out[valid] <- log(c) - log(log(1 + c)) - log(L + c * (x[valid] - a))
  }
  out
}

bradford_1p_cdf <- function(x, fit) {
  c <- fit$shape
  a <- fit$lower
  b <- fit$upper
  L <- b - a
  
  res <- numeric(length(x))
  valid <- !is.na(x) & x >= a & x <= b
  
  if (any(valid)) {
    res[valid] <- log(1 + c * (x[valid] - a) / L) / log(1 + c)
  }
  
  res[x > b] <- 1
  res[x < a] <- 0
  res
}

bradford_1p_logcdf <- function(x, fit) {
  log(bradford_1p_cdf(x, fit))
}

bradford_1p_sf <- function(x, fit) {
  1 - bradford_1p_cdf(x, fit)
}

bradford_1p_logsf <- function(x, fit) {
  log(bradford_1p_sf(x, fit))
}

bradford_1p_quantile <- function(p, fit) {
  c <- fit$shape
  a <- fit$lower
  b <- fit$upper
  L <- b - a
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  a + (L / c) * ((1 + c)^p - 1)
}

bradford_1p_isf <- function(p, fit) {
  bradford_1p_quantile(1 - p, fit)
}

bradford_1p_rand <- function(n, fit) {
  u <- runif(n)
  bradford_1p_quantile(u, fit)
}

# -------------------------------
# Moments
# -------------------------------
bradford_1p_moment <- function(n, fit) {
  # Use numerical integration for general moments
  bradford_1p_expect(function(x) x^n, fit)
}

bradford_1p_mean <- function(fit) {
  c <- fit$shape
  a <- fit$lower
  b <- fit$upper
  L <- b - a
  
  # Standard mean on [0,1] is (c - log(1+c)) / (c * log(1+c))
  # General mean is a + L * E[X_standard]
  m_std <- (c - log(1 + c)) / (c * log(1 + c))
  a + L * m_std
}

bradford_1p_var <- function(fit) {
  c <- fit$shape
  a <- fit$lower
  b <- fit$upper
  L <- b - a
  
  # Standard variance on [0,1]
  # E[X^2] = (c * (c-2) + 2*log(1+c)) / (2 * c^2 * log(1+c)) --- No, let's re-derive or integrate
  m2_std <- bradford_1p_moment(2, list(shape=c, lower=0, upper=1))
  m1_std <- (c - log(1 + c)) / (c * log(1 + c))
  v_std <- m2_std - m1_std^2
  
  L^2 * v_std
}

bradford_1p_std <- function(fit) {
  sqrt(bradford_1p_var(fit))
}

bradford_1p_skew <- function(fit) {
  x <- bradford_1p_rand(100000, fit)
  m3 <- mean((x - mean(x))^3)
  m3 / (sd(x)^3)
}

bradford_1p_kurtosis <- function(fit) {
  x <- bradford_1p_rand(100000, fit)
  m4 <- mean((x - mean(x))^4)
  m4 / (var(x)^2)
}

bradford_1p_median <- function(fit) {
  bradford_1p_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
bradford_1p_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  bradford_1p_quantile(c(alpha, 1 - alpha), fit)
}

bradford_1p_entropy <- function(fit) {
  c <- fit$shape
  a <- fit$lower
  b <- fit$upper
  L <- b - a
  # Entropy of standard Bradford(c) is log(log(1+c)/c) + log(1+c)/2 + 1/2? No.
  # Let's use numerical integration for safety.
  integrand <- function(x) {
    p <- bradford_1p_pdf(x, fit)
    ifelse(p > 0, -p * log(p), 0)
  }
  integrate(integrand, lower = a, upper = b)$value
}

bradford_1p_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * bradford_1p_pdf(x, fit)
  }
  integrate(integrand, lower = fit$lower, upper = fit$upper)$value
}

# -------------------------------
# S3: logLik
# -------------------------------
logLik.bradford_1p <- function(object, ...) {
  structure(object$log_likelihood, df = 1, nobs = object$n, class = "logLik")
}

logLik.truncated_bradford_1p <- function(object, ...) {
  structure(object$log_likelihood, df = 1, nobs = object$n, class = "logLik")
}
