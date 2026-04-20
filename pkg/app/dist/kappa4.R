# Kappa 4 Distribution Module (Hosking 1994)

# -------------------------------
# Log-likelihood
# -------------------------------
kappa4_log_likelihood <- function(data, xi, alpha, k, h) {
  data <- data[!is.na(data)]
  
  if (alpha <= 0 || h <= 0) return(-Inf)
  
  log_dens <- kappa4_logpdf(data, list(xi = xi, alpha = alpha, k = k, h = h))
  
  if (any(is.infinite(log_dens))) {
    # Fallback for extreme cases
    log_dens[is.infinite(log_dens) & log_dens > 0] <- 1e10
    log_dens[is.infinite(log_dens) & log_dens < 0] <- -1e10
  }
  
  sum(log_dens)
}

# -------------------------------
# Fit (MLE via optim)
# -------------------------------
kappa4_fit <- function(data) {
  data <- data[!is.na(data)]
  if (length(data) < 5) stop("Need more data points for 4-parameter fit")
  
  xi_init <- min(data) - 0.1
  alpha_init <- sd(data)
  k_init <- 0.1
  h_init <- 1.0
  
  neg_log_likelihood <- function(params) {
    xi <- params[1]
    alpha <- exp(params[2])
    k <- params[3]
    h <- exp(params[4])
    
    if (abs(k) < 1e-6) k <- 1e-6
    
    val <- -kappa4_log_likelihood(data, xi, alpha, k, h)
    if (!is.finite(val)) return(1e10)
    val
  }
  
  fit <- optim(
    par = c(xi_init, log(alpha_init), k_init, log(h_init)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 2000)
  )
  
  xi_hat <- fit$par[1]
  alpha_hat <- exp(fit$par[2])
  k_hat <- fit$par[3]
  h_hat <- exp(fit$par[4])
  
  log_lik <- kappa4_log_likelihood(data, xi_hat, alpha_hat, k_hat, h_hat)
  
  result <- list(
    xi = xi_hat,
    alpha = alpha_hat,
    k = k_hat,
    h = h_hat,
    log_likelihood = log_lik,
    n = length(data),
    convergence = fit$convergence,
    distribution = "kappa4"
  )
  
  class(result) <- "kappa4"
  result
}

# -------------------------------
# Truncated Fit
# -------------------------------
kappa4_fit_truncated <- function(data, lower = -Inf, upper = Inf) {
  data <- data[!is.na(data) & data >= lower & data <= upper]
  if (length(data) < 5) stop("Need more data points within bounds")
  
  init <- kappa4_fit(data)
  
  neg_log_likelihood <- function(params) {
    xi <- params[1]
    alpha <- exp(params[2])
    k <- params[3]
    h <- exp(params[4])
    
    if (abs(k) < 1e-6) k <- 1e-6
    
    ll_base <- kappa4_log_likelihood(data, xi, alpha, k, h)
    
    fit_tmp <- list(xi=xi, alpha=alpha, k=k, h=h)
    F_upper <- if (is.finite(upper)) kappa4_cdf(upper, fit_tmp) else 1
    F_lower <- if (is.finite(lower)) kappa4_cdf(lower, fit_tmp) else 0
    
    diff <- F_upper - F_lower
    if (diff <= 1e-10) return(1e10)
    
    -(ll_base - length(data) * log(diff))
  }
  
  fit <- optim(
    par = c(init$xi, log(init$alpha), init$k, log(init$h)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 2000)
  )
  
  result <- list(
    xi = fit$par[1],
    alpha = exp(fit$par[2]),
    k = fit$par[3],
    h = exp(fit$par[4]),
    log_likelihood = -fit$value,
    n = length(data),
    lower = lower,
    upper = upper,
    convergence = fit$convergence,
    distribution = "truncated_kappa4"
  )
  class(result) <- "truncated_kappa4"
  result
}

# -------------------------------
# PDF / CDF / Quantile / Rand
# -------------------------------
kappa4_pdf <- function(x, fit) {
  exp(kappa4_logpdf(x, fit))
}

kappa4_logpdf <- function(x, fit) {
  xi <- fit$xi
  alpha <- fit$alpha
  k <- fit$k
  h <- fit$h
  
  out <- rep(-Inf, length(x))
  if (alpha <= 0 || h <= 0) return(out)
  
  y <- 1 - k * (x - xi) / alpha
  valid <- !is.na(y) & y > 0
  
  if (any(valid)) {
    # Check if k is zero (Gumbel limit-ish, though Kappa 4 is different)
    if (abs(k) < 1e-8) {
      # Use limit k -> 0
      z <- (x[valid] - xi) / alpha
      exp_z <- exp(-z)
      term <- 1 - h * exp_z
      pos_term <- term > 0
      
      if (any(pos_term)) {
        out_valid <- numeric(sum(valid))
        out_valid[pos_term] <- -log(alpha) - z[pos_term] + (1/h - 1) * log(term[pos_term])
        out[valid] <- out_valid
      }
    } else {
      term <- 1 - h * y[valid]^(1/k)
      pos_term <- !is.na(term) & term > 0
      
      if (any(pos_term)) {
        out_valid <- rep(-Inf, sum(valid))
        # f(x) = alpha^-1 * y^(1/k - 1) * [F(x)]^(1 - h)
        # F(x) = (1 - h * y^(1/k))^(1/h)
        # log f(x) = -log alpha + (1/k - 1)log y + (1 - h)/h * log(1 - h * y^(1/k))
        out_valid[pos_term] <- -log(alpha) + (1/k - 1) * log(y[valid][pos_term]) + 
                               ((1 - h) / h) * log(term[pos_term])
        out[valid] <- out_valid
      }
    }
  }
  out
}

kappa4_cdf <- function(x, fit) {
  xi <- fit$xi
  alpha <- fit$alpha
  k <- fit$k
  h <- fit$h
  
  res <- numeric(length(x))
  y <- 1 - k * (x - xi) / alpha
  valid_y <- !is.na(y) & y > 0
  
  if (any(valid_y)) {
    if (abs(k) < 1e-8) {
      z <- (x[valid_y] - xi) / alpha
      term <- 1 - h * exp(-z)
    } else {
      term <- 1 - h * y[valid_y]^(1/k)
    }
    
    valid_term <- !is.na(term) & term > 0
    temp_res <- numeric(sum(valid_y))
    if (abs(h) > 1e-10) {
      temp_res[valid_term] <- (term[valid_term])^(1/h)
    }
    temp_res[!valid_term] <- 0
    res[valid_y] <- temp_res
  }
  
  # Support depends on k
  if (k > 0) {
    res[y <= 0 & !is.na(y)] <- 1
  } else if (k < 0) {
    res[y <= 0 & !is.na(y)] <- 0
  }
  
  res <- pmax(0, pmin(1, res))
  res[is.na(res)] <- 0
  res
}

kappa4_logcdf <- function(x, fit) {
  log(kappa4_cdf(x, fit))
}

kappa4_sf <- function(x, fit) {
  1 - kappa4_cdf(x, fit)
}

kappa4_logsf <- function(x, fit) {
  log(kappa4_sf(x, fit))
}

kappa4_quantile <- function(p, fit) {
  xi <- fit$xi
  alpha <- fit$alpha
  k <- fit$k
  h <- fit$h
  
  p <- pmax(1e-10, pmin(1-1e-10, p))
  
  if (abs(k) < 1e-8) {
    # Limit k -> 0
    return(xi - alpha * log((1 - p^h) / h))
  }
  
  xi + (alpha / k) * (1 - ((1 - p^h) / h)^k)
}

kappa4_isf <- function(p, fit) {
  kappa4_quantile(1 - p, fit)
}

kappa4_rand <- function(n, fit) {
  u <- runif(n)
  kappa4_quantile(u, fit)
}

# -------------------------------
# Moments
# -------------------------------
kappa4_moment <- function(n, fit) {
  set.seed(42)
  mean(kappa4_rand(100000, fit)^n)
}

kappa4_mean <- function(fit) {
  set.seed(42)
  mean(kappa4_rand(100000, fit))
}

kappa4_var <- function(fit) {
  set.seed(42)
  x <- kappa4_rand(100000, fit)
  var(x)
}

kappa4_std <- function(fit) {
  sqrt(kappa4_var(fit))
}

kappa4_skew <- function(fit) {
  set.seed(42)
  x <- kappa4_rand(100000, fit)
  m3 <- mean((x - mean(x))^3)
  m3 / (sd(x)^3)
}

kappa4_kurtosis <- function(fit) {
  set.seed(42)
  x <- kappa4_rand(100000, fit)
  m4 <- mean((x - mean(x))^4)
  m4 / (var(x)^2)
}

kappa4_median <- function(fit) {
  kappa4_quantile(0.5, fit)
}

# -------------------------------
# Interval / Entropy / Expect
# -------------------------------
kappa4_interval <- function(level, fit) {
  alpha_val <- (1 - level) / 2
  kappa4_quantile(c(alpha_val, 1 - alpha_val), fit)
}

kappa4_entropy <- function(fit) {
  integrand <- function(x) {
    pdf_val <- kappa4_pdf(x, fit)
    res <- numeric(length(x))
    pos <- !is.na(pdf_val) & pdf_val > 0
    res[pos] <- -pdf_val[pos] * log(pdf_val[pos])
    res
  }
  # Determine integration bounds from quantiles
  lower <- kappa4_quantile(1e-7, fit)
  upper <- kappa4_quantile(1 - 1e-7, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}

kappa4_expect <- function(func, fit, ...) {
  integrand <- function(x) {
    func(x, ...) * kappa4_pdf(x, fit)
  }
  lower <- kappa4_quantile(1e-7, fit)
  upper <- kappa4_quantile(1 - 1e-7, fit)
  res <- integrate(integrand, lower = lower, upper = upper, rel.tol = 1e-6)
  res$value
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.kappa4 <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}

logLik.truncated_kappa4 <- function(object, ...) {
  structure(object$log_likelihood, df = 4, nobs = object$n, class = "logLik")
}
