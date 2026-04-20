# Double Pareto-Lognormal (dPLN) Distribution Module
# Based on Reed and Jorgensen (2004)

# -------------------------------
# Helper: A(theta, nu, tau)
# -------------------------------
.dpln_A <- function(theta, nu, tau) {
  exp(theta * nu + (theta^2 * tau^2) / 2)
}

# -------------------------------
# Log-likelihood (pure)
# -------------------------------
dpln_log_likelihood <- function(data, alpha, beta, nu, tau) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) == 0 || alpha <= 0 || beta <= 0 || tau <= 0) {
    return(-Inf)
  }

  log_dens <- dpln_logpdf(data, list(alpha = alpha, beta = beta, nu = nu, tau = tau))
  
  if (any(is.infinite(log_dens))) {
    # Fallback for extreme cases
    log_dens[is.infinite(log_dens) & log_dens > 0] <- 1e10
    log_dens[is.infinite(log_dens) & log_dens < 0] <- -1e10
  }
  
  sum(log_dens)
}


# -------------------------------
# Fit dPLN (MLE via optim)
# -------------------------------
dpln_fit <- function(data) {
  data <- data[!is.na(data)]
  data <- data[data > 0]

  if (length(data) < 4) {
    stop("Need at least 4 valid positive data points for dPLN (4 parameters)")
  }

  # Initial guesses based on log-data
  log_data <- log(data)
  initial_nu <- mean(log_data)
  initial_tau <- sd(log_data)
  initial_alpha <- 2.0
  initial_beta <- 2.0

  neg_log_likelihood <- function(params) {
    alpha <- exp(params[1])
    beta  <- exp(params[2])
    nu    <- params[3]
    tau   <- exp(params[4])

    if (is.na(tau) || tau < 1e-8 || alpha < 1e-8 || beta < 1e-8) return(1e20)

    ll <- dpln_log_likelihood(data, alpha, beta, nu, tau)
    if (!is.finite(ll)) return(1e20)
    -ll
  }

  fit <- optim(
    par = c(log(initial_alpha), log(initial_beta), initial_nu, log(initial_tau)),
    fn = neg_log_likelihood,
    method = "L-BFGS-B",
    lower = c(log(1e-4), log(1e-4), -20, log(1e-4)),
    upper = c(log(20), log(20), 20, log(10)),
    control = list(maxit = 2000)
  )

  alpha_hat <- exp(fit$par[1])
  beta_hat  <- exp(fit$par[2])
  nu_hat    <- fit$par[3]
  tau_hat   <- exp(fit$par[4])

  log_lik <- dpln_log_likelihood(data, alpha_hat, beta_hat, nu_hat, tau_hat)

  result <- list(
    alpha = alpha_hat,
    beta = beta_hat,
    nu = nu_hat,
    tau = tau_hat,
    log_likelihood = log_lik,
    n = length(data),
    distribution = "dpln"
  )

  class(result) <- "dpln"
  result
}


# -------------------------------
# Truncated dPLN (MLE via optim)
# -------------------------------
dpln_log_likelihood_truncated <- function(data, alpha, beta, nu, tau, lower, upper) {
  if (alpha <= 0 || beta <= 0 || tau <= 0 || lower < 0) return(-Inf)

  data <- data[data >= lower & data <= upper]
  if (length(data) == 0) return(-Inf)

  fit_tmp <- list(alpha = alpha, beta = beta, nu = nu, tau = tau)
  ll <- sum(dpln_logpdf(data, fit_tmp))

  p_upper <- dpln_cdf(upper, fit_tmp)
  p_lower <- dpln_cdf(lower, fit_tmp)

  diff <- p_upper - p_lower
  if (diff <= 1e-12) return(-Inf)

  ll - length(data) * log(diff)
}


dpln_fit_truncated <- function(data, lower = 0, upper = Inf) {
  data <- data[!is.na(data)]
  data <- data[data >= lower & data <= upper]

  if (length(data) < 4) {
    stop("Need at least 4 data points within truncation bounds")
  }

  log_data <- log(data)
  initial_nu <- mean(log_data)
  initial_tau <- sd(log_data)
  initial_alpha <- 2.0
  initial_beta <- 2.0

  neg_log_likelihood <- function(params) {
    alpha <- exp(params[1])
    beta  <- exp(params[2])
    nu    <- params[3]
    tau   <- exp(params[4])

    if (tau < 1e-8 || alpha < 1e-8 || beta < 1e-8) return(Inf)

    -dpln_log_likelihood_truncated(data, alpha, beta, nu, tau, lower, upper)
  }

  fit <- optim(
    par = c(log(initial_alpha), log(initial_beta), initial_nu, log(initial_tau)),
    fn = neg_log_likelihood,
    method = "BFGS",
    control = list(maxit = 2000)
  )

  alpha_hat <- exp(fit$par[1])
  beta_hat  <- exp(fit$par[2])
  nu_hat    <- fit$par[3]
  tau_hat   <- exp(fit$par[4])

  log_lik <- dpln_log_likelihood_truncated(
    data, alpha_hat, beta_hat, nu_hat, tau_hat, lower, upper
  )

  result <- list(
    alpha = alpha_hat,
    beta = beta_hat,
    nu = nu_hat,
    tau = tau_hat,
    log_likelihood = log_lik,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_dpln"
  )

  class(result) <- "truncated_dpln"
  result
}


# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
dpln_pdf <- function(x, fit) {
  exp(dpln_logpdf(x, fit))
}

dpln_logpdf <- function(x, fit) {
  alpha <- fit$alpha
  beta <- fit$beta
  nu <- fit$nu
  tau <- fit$tau

  out <- rep(-Inf, length(x))
  valid <- x > 0 & !is.na(x)
  xv <- x[valid]
  lx <- log(xv)

  ln_term1 <- (-alpha - 1) * lx + alpha * nu + (alpha^2 * tau^2) / 2 + 
              pnorm((lx - nu - alpha * tau^2) / tau, log.p = TRUE)
  ln_term2 <- (beta - 1) * lx - beta * nu + (beta^2 * tau^2) / 2 + 
              pnorm(-(lx - nu + beta * tau^2) / tau, log.p = TRUE)
  
  # log-sum-exp trick
  ln_max <- pmax(ln_term1, ln_term2)
  ln_sum <- ln_max + log(exp(ln_term1 - ln_max) + exp(ln_term2 - ln_max))
  ln_sum[is.nan(ln_sum)] <- -Inf
  
  out[valid] <- log(alpha * beta / (alpha + beta)) + ln_sum
  out
}

dpln_cdf <- function(x, fit) {
  alpha <- fit$alpha
  beta <- fit$beta
  nu <- fit$nu
  tau <- fit$tau

  out <- numeric(length(x))
  valid <- x > 0 & !is.na(x)
  if (!any(valid)) return(out)
  
  xv <- x[valid]
  lx <- log(xv)

  # term1 = P(N(nu, tau^2) <= log(x))
  term1 <- pnorm((lx - nu) / tau)
  
  # Stable term2: beta/(alpha+beta) * exp(alpha*nu + alpha^2*tau^2/2 - alpha*log(xv)) * pnorm(...)
  log_term2 <- log(beta) - log(alpha + beta) + alpha * nu + (alpha^2 * tau^2) / 2 - 
               alpha * lx + pnorm((lx - nu - alpha * tau^2) / tau, log.p = TRUE)
  term2 <- exp(log_term2)
  
  # Stable term3: alpha/(alpha+beta) * exp(-beta*nu + beta^2*tau^2/2 + beta*log(xv)) * pnorm(...)
  log_term3 <- log(alpha) - log(alpha + beta) - beta * nu + (beta^2 * tau^2) / 2 + 
               beta * lx + pnorm(-(lx - nu + beta * tau^2) / tau, log.p = TRUE)
  term3 <- exp(log_term3)
  
  # Heuristic: if exp() overflowed to Inf, it's likely a numerical artifact 
  # for a term that should be cancelling out or reaching a boundary.
  # Given the property that CDF must be in [0,1], we clamp.
  res_valid <- term1 - term2 + term3
  res_valid[is.na(res_valid) | is.nan(res_valid)] <- term1[is.na(res_valid) | is.nan(res_valid)] # Fallback to lognormal base
  
  out[valid] <- res_valid
  
  # Ensure bounds
  out[x <= 0] <- 0
  out[out < 0] <- 0
  out[out > 1] <- 1
  out
}

dpln_logcdf <- function(x, fit) {
  log(dpln_cdf(x, fit))
}

dpln_sf <- function(x, fit) {
  1 - dpln_cdf(x, fit)
}

dpln_logsf <- function(x, fit) {
  log(dpln_sf(x, fit))
}

dpln_quantile <- function(p, fit) {
  if (any(p < 0 | p > 1)) {
    stop("Probabilities must be in [0,1]")
  }

  res <- sapply(p, function(prob) {
    if (is.na(prob)) return(NA)
    if (prob == 0) return(0)
    if (prob == 1) return(Inf)
    
    # Use uniroot to find where CDF(x) - p = 0
    # Search range: use lognormal quantiles as bounds
    lower_bound <- qlnorm(max(1e-10, prob/10), meanlog = fit$nu, sdlog = fit$tau)
    upper_bound <- qlnorm(min(1-1e-10, 1 - (1-prob)/10), meanlog = fit$nu, sdlog = fit$tau)
    
    # Expand bounds if necessary
    while(dpln_cdf(lower_bound, fit) > prob && lower_bound > 1e-20) {
      lower_bound <- lower_bound / 10
    }
    while(dpln_cdf(upper_bound, fit) < prob && upper_bound < 1e20) {
      upper_bound <- upper_bound * 10
    }

    tryCatch({
      uniroot(function(x) dpln_cdf(x, fit) - prob, 
              lower = lower_bound, 
              upper = upper_bound, 
              tol = 1e-8)$root
    }, error = function(e) {
      # Fallback to optim if uniroot fails
      optim(par = fit$nu, 
            fn = function(lx) (dpln_cdf(exp(lx), fit) - prob)^2,
            method = "BFGS")$par |> exp()
    })
  })
  res
}

dpln_isf <- function(p, fit) {
  dpln_quantile(1 - p, fit)
}

dpln_rand <- function(n, fit) {
  # X = exp(nu + tau * Z + E1 - E2)
  # where Z ~ N(0,1), E1 ~ Exp(alpha), E2 ~ Exp(beta)
  z  <- rnorm(n, mean = 0, sd = 1)
  e1 <- rexp(n, rate = fit$alpha)
  e2 <- rexp(n, rate = fit$beta)
  
  exp(fit$nu + fit$tau * z + e1 - e2)
}


# -------------------------------
# Moments
# -------------------------------
dpln_moment <- function(n, fit) {
  if (fit$alpha <= n) return(Inf)
  exp(n * fit$nu + (n^2 * fit$tau^2) / 2) * (fit$alpha / (fit$alpha - n)) * (fit$beta / (fit$beta + n))
}

dpln_mean <- function(fit) {
  dpln_moment(1, fit)
}

dpln_var <- function(fit) {
  if (fit$alpha <= 2) return(Inf)
  dpln_moment(2, fit) - dpln_mean(fit)^2
}

dpln_std <- function(fit) {
  sqrt(dpln_var(fit))
}

dpln_skew <- function(fit) {
  if (fit$alpha <= 3) return(NA)
  m1 <- dpln_moment(1, fit)
  m2 <- dpln_moment(2, fit)
  m3 <- dpln_moment(3, fit)
  sigma <- sqrt(m2 - m1^2)
  (m3 - 3*m1*m2 + 2*m1^3) / sigma^3
}

dpln_kurtosis <- function(fit) {
  if (fit$alpha <= 4) return(NA)
  m1 <- dpln_moment(1, fit)
  m2 <- dpln_moment(2, fit)
  m3 <- dpln_moment(3, fit)
  m4 <- dpln_moment(4, fit)
  sigma2 <- m2 - m1^2
  (m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4) / sigma2^2
}

dpln_median <- function(fit) {
  dpln_quantile(0.5, fit)
}


# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
dpln_interval <- function(level, fit) {
  alpha_val <- (1 - level) / 2
  dpln_quantile(c(alpha_val, 1 - alpha_val), fit)
}

dpln_entropy <- function(fit) {
  integrand <- function(x) {
    pdf_val <- dpln_pdf(x, fit)
    res <- numeric(length(x))
    pos <- !is.na(pdf_val) & pdf_val > 0
    res[pos] <- -pdf_val[pos] * log(pdf_val[pos])
    res
  }
  lower <- dpln_quantile(1e-7, fit)
  upper <- dpln_quantile(1 - 1e-7, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}

dpln_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * dpln_pdf(x, fit)
  }
  lower <- dpln_quantile(1e-7, fit)
  upper <- dpln_quantile(1 - 1e-7, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}


# -------------------------------
# S3: logLik (enables AIC/BIC)
# -------------------------------
logLik.dpln <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 4,
    nobs = object$n,
    class = "logLik"
  )
}

logLik.truncated_dpln <- function(object, ...) {
  structure(
    object$log_likelihood,
    df = 4,
    nobs = object$n,
    class = "logLik"
  )
}
