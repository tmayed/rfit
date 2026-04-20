# Pareto Distribution Module (Improved)

# -------------------------------
# Log-likelihood (standard Pareto)
# -------------------------------
pareto_log_likelihood <- function(data, shape, scale) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) == 0 || shape <= 0 || scale <= 0) {
    return(-Inf)
  }

  if (any(data < scale)) {
    return(-Inf)
  }

  sum(log(shape) + shape * log(scale) - (shape + 1) * log(data))
}


# -------------------------------
# Fit standard Pareto (closed-form MLE)
# -------------------------------
pareto_fit <- function(data) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) < 2) {
    stop("Need at least 2 valid positive data points")
  }

  scale_hat <- min(data)

  denom <- sum(log(data / scale_hat))
  if (denom <= 0) {
    stop("Numerical instability in shape estimation")
  }

  shape_hat <- length(data) / denom

  log_lik <- pareto_log_likelihood(data, shape_hat, scale_hat)

  result <- list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "pareto"
  )

  class(result) <- "pareto"
  result
}


# -------------------------------
# Truncated Pareto (proper MLE)
# -------------------------------
pareto_log_likelihood_truncated <- function(data, shape, scale, lower, upper) {
  if (shape <= 0 || scale <= 0) return(-Inf)
  if (scale > lower) return(-Inf)  # IMPORTANT constraint

  data <- data[data >= lower & data <= upper]

  if (length(data) == 0) return(-Inf)
  if (any(data < scale)) return(-Inf)

  n <- length(data)

  # Standard log-likelihood
  ll <- sum(log(shape) + shape * log(scale) - (shape + 1) * log(data))

  # Truncation normalization
  if (is.finite(upper)) {
    norm <- log((scale / lower)^shape - (scale / upper)^shape)
  } else {
    norm <- log((scale / lower)^shape)
  }

  ll - n * norm
}


pareto_fit_truncated <- function(data, lower = NULL, upper = Inf) {
  data <- data[!is.na(data)]

  if (is.null(lower)) {
    lower <- min(data)
  }

  data <- data[data >= lower & data <= upper]

  if (length(data) < 2) {
    stop("Need at least 2 data points within truncation bounds")
  }

  # Initial guesses
  scale_init <- min(data) * 0.9  # allow scale < lower
  shape_init <- length(data) / sum(log(data / min(data)))

  neg_log_likelihood <- function(params) {
    shape <- exp(params[1])
    scale <- exp(params[2])

    # enforce scale <= lower
    if (scale > lower) return(Inf)

    -pareto_log_likelihood_truncated(data, shape, scale, lower, upper)
  }

  fit <- optim(
    par = c(log(shape_init), log(scale_init)),
    fn = neg_log_likelihood,
    method = "L-BFGS-B",
    lower = c(log(1e-8), log(1e-8)),
    upper = c(log(1e6), log(lower)),
    control = list(maxit = 1000)
  )

  shape_hat <- exp(fit$par[1])
  scale_hat <- exp(fit$par[2])

  log_lik <- pareto_log_likelihood_truncated(
    data, shape_hat, scale_hat, lower, upper
  )

  result <- list(
    shape = shape_hat,
    scale = scale_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_pareto"
  )

  class(result) <- "truncated_pareto"
  result
}


# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
pareto_pdf <- function(x, fit) {
  out <- numeric(length(x))
  valid <- !is.na(x) & x >= fit$scale

  out[valid] <- fit$shape * (fit$scale^fit$shape) / x[valid]^(fit$shape + 1)
  out
}


pareto_cdf <- function(x, fit) {
  ifelse(x < fit$scale, 0, 1 - (fit$scale / x)^fit$shape)
}


pareto_quantile <- function(p, fit) {
  if (any(p < 0 | p > 1)) {
    stop("Probabilities must be in [0,1]")
  }

  fit$scale / (1 - p)^(1 / fit$shape)
}


pareto_rand <- function(n, fit) {
  u <- runif(n)
  fit$scale / (1 - u)^(1 / fit$shape)
}


# -------------------------------
# Survival / Inverse Survival / Logs
# -------------------------------
pareto_sf <- function(x, fit) {
  ifelse(x < fit$scale, 1, (fit$scale / x)^fit$shape)
}

pareto_isf <- function(p, fit) {
  fit$scale / p^(1 / fit$shape)
}

pareto_logpdf <- function(x, fit) {
  out <- rep(-Inf, length(x))
  valid <- !is.na(x) & x >= fit$scale
  out[valid] <- log(fit$shape) + fit$shape * log(fit$scale) - (fit$shape + 1) * log(x[valid])
  out
}

pareto_logcdf <- function(x, fit) {
  log(pareto_cdf(x, fit))
}

pareto_logsf <- function(x, fit) {
  ifelse(x < fit$scale, 0, fit$shape * (log(fit$scale) - log(x)))
}


# -------------------------------
# Moments
# -------------------------------
pareto_mean <- function(fit) {
  if (fit$shape <= 1) return(Inf)
  fit$scale * fit$shape / (fit$shape - 1)
}

pareto_var <- function(fit) {
  if (fit$shape <= 2) return(Inf)
  fit$scale^2 * fit$shape / ((fit$shape - 1)^2 * (fit$shape - 2))
}

pareto_std <- function(fit) {
  sqrt(pareto_var(fit))
}

pareto_moment <- function(n, fit) {
  if (fit$shape <= n) return(Inf)
  fit$shape * fit$scale^n / (fit$shape - n)
}

pareto_skew <- function(fit) {
  if (fit$shape <= 3) return(NA)
  (2 * (1 + fit$shape) / (fit$shape - 3)) * sqrt((fit$shape - 2) / fit$shape)
}

pareto_kurtosis <- function(fit) {
  if (fit$shape <= 4) return(NA)
  3 * (fit$shape - 2) * (3 * fit$shape^2 + fit$shape + 2) / (fit$shape * (fit$shape - 3) * (fit$shape - 4))
}

pareto_median <- function(fit) {
  fit$scale * 2^(1 / fit$shape)
}

pareto_interval <- function(level, fit) {
  alpha <- (1 - level) / 2
  pareto_quantile(c(alpha, 1 - alpha), fit)
}


# -------------------------------
# Entropy & Expect
# -------------------------------
pareto_entropy <- function(fit) {
  log(fit$scale / fit$shape) + 1 / fit$shape + 1
}

pareto_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * pareto_pdf(x, fit)
  }
  # Pareto support is [scale, Inf)
  res <- integrate(integrand, lower = fit$scale, upper = Inf)
  res$value
}


# -------------------------------
# S3: logLik (enables AIC/BIC)
# -------------------------------
logLik.pareto <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 2,
    nobs = object$n,
    class = "logLik"
  )
}

logLik.truncated_pareto <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 2,
    nobs = object$n,
    class = "logLik"
  )
}
