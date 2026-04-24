#' @title fit_ln_gamma_dpnl
#' @description Fit a three-component mixture of lognormal, gamma, and dPLN
#'   distributions by maximum likelihood.
#' @param data A numeric vector of positive observations.
#' @param initial_ln Optional named list with lognormal starting values
#'   `mu` and `sigma`.
#' @param initial_gamma Optional named list with gamma starting values
#'   `shape` and `scale`.
#' @param initial_dpnl Optional named list with dPLN starting values
#'   `alpha`, `beta`, `nu`, and `tau`.
#' @param n_starts Number of optimization starts. Defaults to 5.
#' @param maxit Maximum optimizer iterations per start. Defaults to 3000.
#' @return A fitted mixture object of class `ln_gamma_dpnl_mixture`.
#' @details
#' This implementation is fully self-contained. It uses a multi-start strategy
#' with a baseline initialized from data-driven heuristics. Parameters are 
#' optimized on an unconstrained internal scale. The function requires at 
#' least one start to converge (code 0) or it will stop with an error.
#'
#' @section Design Constraints & Maintenance Mandate:
#' To ensure the mathematical integrity and physical interpretability of the 
#' mixture model, the following constraints are STRICTLY enforced and MUST be 
#' preserved in all future updates:
#' \itemize{
#'   \item \bold{alpha > 1}: Required to ensure the dPLN component mean is finite and well-defined.
#'   \item \bold{beta > 1}: Required to ensure the dPLN PDF $f(0) = 0$.
#'   \item \bold{gamma shape > 1}: Required to ensure the Gamma PDF $f(0) = 0$.
#' }
#' These conditions are implemented using a `1 + exp(theta)` transformation in 
#' the unpacking logic. Any modification to the parameter mapping must 
#' explicitly maintain these lower bounds.
#' @export
fit_ln_gamma_dpnl <- function(data,
                          initial_ln = NULL,
                          initial_gamma = NULL,
                          initial_dpnl = NULL,
                          n_starts = 5L,
                          maxit = 3000L) {
  # --- 1. Data Cleaning ---
  data <- data[is.finite(data) & !is.na(data) & data > 0]
  if (length(data) < 10L) {
    stop("Need at least 10 positive observations.")
  }

  n_starts <- as.integer(n_starts[1])
  maxit <- as.integer(maxit[1])
  if (is.na(n_starts) || n_starts < 1L) stop("`n_starts` must be an integer >= 1.")
  if (is.na(maxit) || maxit < 1L) stop("`maxit` must be an integer >= 1.")
  
  # --- 2. Initializer Validation ---
  if (!is.null(initial_ln)) {
    if (!is.list(initial_ln) || is.null(initial_ln$mu) || is.null(initial_ln$sigma))
      stop("`initial_ln` must be a list with `mu` and `sigma`.")
    if (!is.finite(initial_ln$mu) || !is.finite(initial_ln$sigma) || initial_ln$sigma <= 0)
      stop("`initial_ln` parameters must be finite and `sigma > 0`.")
  }
  if (!is.null(initial_gamma)) {
    if (!is.list(initial_gamma) || is.null(initial_gamma$shape) || is.null(initial_gamma$scale))
      stop("`initial_gamma` must be a list with `shape` and `scale`.")
    if (!is.finite(initial_gamma$shape) || !is.finite(initial_gamma$scale) || initial_gamma$shape <= 1 || initial_gamma$scale <= 0)
      stop("`initial_gamma` parameters must be finite, `shape > 1` and `scale > 0`.")
  }
  if (!is.null(initial_dpnl)) {
    if (!is.list(initial_dpnl) || is.null(initial_dpnl$alpha) || is.null(initial_dpnl$beta) || is.null(initial_dpnl$nu) || is.null(initial_dpnl$tau))
      stop("`initial_dpnl` must be a list with `alpha`, `beta`, `nu`, and `tau`.")
    if (!is.finite(initial_dpnl$alpha) || !is.finite(initial_dpnl$beta) || !is.finite(initial_dpnl$nu) || !is.finite(initial_dpnl$tau) ||
        initial_dpnl$alpha <= 1 || initial_dpnl$beta <= 1 || initial_dpnl$tau <= 0)
      stop("`initial_dpnl` parameters must be finite, `alpha > 1`, `beta > 1`, and `tau > 0`.")
  }

  # --- 3. Self-Contained Initialization ---
  data_sorted <- sort(data)
  n <- length(data)
  idx1 <- 1:floor(n/3)
  idx2 <- (floor(n/3)+1):floor(2*n/3)
  idx3 <- (floor(2*n/3)+1):n
  
  ln_mu <- mean(log(data_sorted[idx1]))
  ln_sigma <- stats::sd(log(data_sorted[idx1]))
  if (is.na(ln_sigma) || ln_sigma < 1e-3) ln_sigma <- 0.1
  if (!is.null(initial_ln)) {
    ln_mu <- initial_ln$mu
    ln_sigma <- initial_ln$sigma
  }

  m2 <- mean(data_sorted[idx2])
  v2 <- stats::var(data_sorted[idx2])
  if (is.na(v2) || v2 < 1e-6) v2 <- m2^2 * 0.1
  g_shape <- max(m2^2 / v2, 1.1)
  g_scale <- v2 / m2
  if (!is.null(initial_gamma)) {
    g_shape <- initial_gamma$shape
    g_scale <- initial_gamma$scale
  }

  d_nu <- mean(log(data_sorted[idx3]))
  d_tau <- stats::sd(log(data_sorted[idx3]))
  if (is.na(d_tau) || d_tau < 1e-3) d_tau <- 0.1
  d_alpha <- 2
  d_beta <- 2
  if (!is.null(initial_dpnl)) {
    d_alpha <- initial_dpnl$alpha
    d_beta <- initial_dpnl$beta
    d_nu <- initial_dpnl$nu
    d_tau <- initial_dpnl$tau
  }

  # --- 4. Optimization Setup ---
  base_start <- c(
    0, 0,
    ln_mu, log(ln_sigma),
    log(g_shape - 1), log(g_scale),
    log(d_alpha - 1), log(d_beta - 1),
    d_nu, log(d_tau)
  )

  starts <- list(base_start)
  if (n_starts > 1) {
    set.seed(42)
    for (i in 2:n_starts) {
      starts[[i]] <- base_start + stats::rnorm(length(base_start), sd = 0.25)
    }
  }

  # --- 5. Run Optimizations ---
  best_fit <- NULL
  best_value <- Inf

  for (start in starts) {
    fit_try <- tryCatch(
      optim(
        par = start,
        fn = .fit_ln_gamma_dpnl_neg_log_likelihood,
        data = data,
        method = "Nelder-Mead",
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_try) && is.finite(fit_try$value) && fit_try$convergence == 0) {
      if (fit_try$value < best_value) {
        best_fit <- fit_try
        best_value <- fit_try$value
      }
    }
  }

  if (is.null(best_fit)) {
    stop("Optimization failed to converge for all starting values.")
  }

  # --- 6. Result Assembly ---
  final <- .fit_ln_gamma_dpnl_unpack(best_fit$par)
  log_lik <- .fit_ln_gamma_dpnl_mixture_log_likelihood(best_fit$par, data, penalty = FALSE)

  result <- list(
    weights = stats::setNames(final$weights, c("lognormal", "gamma", "dpln")),
    components = list(
      lognormal = final$lognormal,
      gamma = final$gamma,
      dpln = final$dpln
    ),
    log_likelihood = log_lik,
    n = length(data),
    distribution = "ln_gamma_dpnl_mixture",
    convergence = best_fit$convergence,
    params_internal = best_fit$par
  )

  class(result) <- "ln_gamma_dpnl_mixture"
  result
}

#' @export
logLik.ln_gamma_dpnl_mixture <- function(object, ...) {
  structure(object$log_likelihood, df = 10, nobs = object$n, class = "logLik")
}

.fit_ln_gamma_dpnl_unpack <- function(params) {
  # Logit weights
  w_logit <- pmax(-20, pmin(20, params[1:2]))
  w_raw <- c(w_logit, 0)
  weights <- exp(w_raw) / sum(exp(w_raw))

  # Parameters with 1 + exp() constraints
  list(
    weights = weights,
    lognormal = list(
      mu = params[3],
      sigma = exp(pmax(-15, pmin(15, params[4])))
    ),
    gamma = list(
      shape = 1 + exp(pmax(-15, pmin(15, params[5]))),
      scale = exp(pmax(-15, pmin(15, params[6])))
    ),
    dpln = list(
      alpha = 1 + exp(pmax(-15, pmin(15, params[7]))),
      beta = 1 + exp(pmax(-15, pmin(15, params[8]))),
      nu = params[9],
      tau = exp(pmax(-15, pmin(15, params[10])))
    )
  )
}

.fit_ln_gamma_dpnl_mixture_log_likelihood <- function(params, data, penalty = TRUE) {
  fit <- .fit_ln_gamma_dpnl_unpack(params)
  
  log_ln <- dlnorm(data, meanlog = fit$lognormal$mu, sdlog = fit$lognormal$sigma, log = TRUE)
  log_g  <- dgamma(data, shape = fit$gamma$shape, scale = fit$gamma$scale, log = TRUE)
  
  lx <- log(data)
  alpha <- fit$dpln$alpha; beta <- fit$dpln$beta; nu <- fit$dpln$nu; tau <- fit$dpln$tau
  
  z1 <- (lx - nu - alpha * tau^2) / tau
  z2 <- -(lx - nu + beta * tau^2) / tau
  
  # Log-domain terms for dPLN
  ln_term1 <- (-alpha - 1) * lx + alpha * nu + (alpha^2 * tau^2) / 2 + pnorm(z1, log.p = TRUE)
  ln_term2 <- (beta - 1) * lx - beta * nu + (beta^2 * tau^2) / 2 + pnorm(z2, log.p = TRUE)
  
  ln_max_d <- pmax(ln_term1, ln_term2)
  log_dp <- log(alpha * beta / (alpha + beta)) + ln_max_d + log(exp(ln_term1 - ln_max_d) + exp(ln_term2 - ln_max_d))
  
  comp_log <- cbind(
    log(fit$weights[1]) + log_ln,
    log(fit$weights[2]) + log_g,
    log(fit$weights[3]) + log_dp
  )
  
  row_max <- pmax(comp_log[,1], comp_log[,2], comp_log[,3])
  valid_rows <- is.finite(row_max)
  
  ll_vec <- rep(-1e20, length(data))
  if (any(valid_rows)) {
    ll_vec[valid_rows] <- row_max[valid_rows] + log(
      exp(comp_log[valid_rows, 1] - row_max[valid_rows]) + 
      exp(comp_log[valid_rows, 2] - row_max[valid_rows]) + 
      exp(comp_log[valid_rows, 3] - row_max[valid_rows])
    )
  }
  
  ll <- sum(ll_vec)
  if (!is.finite(ll)) return(-1e20)
  
  if (penalty) {
    # Regularization to keep parameters in sane regions and encourage 
    # identifiability by penalizing extreme overlaps or scale collapses.
    return(ll - 5e-4 * sum(params^2))
  }
  ll
}

.fit_ln_gamma_dpnl_neg_log_likelihood <- function(params, data) {
  - .fit_ln_gamma_dpnl_mixture_log_likelihood(params, data, penalty = TRUE)
}
