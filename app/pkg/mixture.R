# Mixture Distribution Module - ROBUST VERSION

# -------------------------------
# Parameter Definitions
# -------------------------------
.dist_params <- list(
  lognormal = list(
    names = c("mu", "sigma"),
    to_internal = function(p) c(p$mu, log(p$sigma)),
    from_internal = function(v) list(mu = v[1], sigma = exp(v[2]))
  ),
  normal = list(
    names = c("mean", "sd"),
    to_internal = function(p) c(p$mean, log(p$sd)),
    from_internal = function(v) list(mean = v[1], sd = exp(v[2]))
  ),
  weibull = list(
    names = c("shape", "scale"),
    to_internal = function(p) c(log(p$shape), log(p$scale)),
    from_internal = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  pareto = list(
    names = c("shape", "scale"),
    to_internal = function(p) c(log(p$shape), log(p$scale)),
    from_internal = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  dpln = list(
    names = c("alpha", "beta", "nu", "tau"),
    to_internal = function(p) c(log(p$alpha), log(p$beta), p$nu, log(p$tau)),
    from_internal = function(v) list(alpha = exp(v[1]), beta = exp(v[2]), nu = v[3], tau = exp(v[4]))
  ),
  gb2 = list(
    names = c("a", "b", "p", "q"),
    to_internal = function(p) c(log(p$a), log(p$b), log(p$p), log(p$q)),
    from_internal = function(v) list(a = exp(v[1]), b = exp(v[2]), p = exp(v[3]), q = exp(v[4]))
  ),
  kappa4 = list(
    names = c("xi", "alpha", "k", "h"),
    to_internal = function(p) c(p$xi, log(p$alpha), p$k, log(p$h)),
    from_internal = function(v) list(xi = v[1], alpha = exp(v[2]), k = v[3], h = exp(v[4]))
  ),
  bradford = list(
    names = c("shape", "lower", "upper"),
    to_internal = function(p) c(log(p$shape), p$lower, p$upper),
    from_internal = function(v) list(shape = exp(v[1]), lower = v[2], upper = v[3])
  ),
  fisk = list(
    names = c("scale", "shape"),
    to_internal = function(p) c(log(p$scale), log(p$shape)),
    from_internal = function(v) list(scale = exp(v[1]), shape = exp(v[2]))
  ),
  johnsonsb = list(
    names = c("gamma", "delta", "xi", "lambda"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi, log(p$lambda)),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3], lambda = exp(v[4]))
  ),
  johnsonsl = list(
    names = c("gamma", "delta", "xi"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3])
  ),
  johnsonsu = list(
    names = c("gamma", "delta", "xi", "lambda"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi, log(p$lambda)),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3], lambda = exp(v[4]))
  )
)

# -------------------------------
# Mixture Fit
# -------------------------------
mixture_fit <- function(data, dist_names) {
  data <- data[!is.na(data) & data > 0]
  K <- length(dist_names)
  
  # Initialization: Simple split
  initial_dist_params_list <- list()
  param_counts <- numeric(K)
  data_sorted <- sort(data)
  
  for (i in 1:K) {
    name <- dist_names[i]
    fit_func <- get(paste0(name, "_fit"))
    # Split data into K chunks
    idx <- floor((i-1)*length(data)/K + 1) : floor(i*length(data)/K)
    fit_obj <- tryCatch(fit_func(data_sorted[idx]), error = function(e) fit_func(data))
    initial_dist_params_list[[i]] <- .dist_params[[name]]$to_internal(fit_obj)
    param_counts[i] <- length(initial_dist_params_list[[i]])
  }
  
  unpack <- function(params) {
    w_raw <- c(params[1:(K-1)], 0)
    w_raw <- pmax(-20, pmin(20, w_raw))
    weights <- exp(w_raw) / sum(exp(w_raw))
    dist_fits <- list()
    curr <- K
    for (i in 1:K) {
        p_vec <- params[curr:(curr + param_counts[i] - 1)]
        # STRICT BOUNDS ON INTERNAL PARAMS
        p_vec <- pmax(-8, pmin(8, p_vec))
        dist_fits[[i]] <- .dist_params[[dist_names[i]]]$from_internal(p_vec)
        dist_fits[[i]]$distribution <- dist_names[i]
        curr <- curr + param_counts[i]
    }
    list(weights = weights, dist_fits = dist_fits)
  }
  
  neg_log_likelihood <- function(params) {
    unpacked <- tryCatch(unpack(params), error = function(e) NULL)
    if (is.null(unpacked)) return(1e15)
    
    total_dens <- numeric(length(data))
    for (i in 1:K) {
      pdf_func <- get(paste0(dist_names[i], "_pdf"))
      dens <- tryCatch(pdf_func(data, unpacked$dist_fits[[i]]), error = function(e) rep(0, length(data)))
      dens[!is.finite(dens)] <- 0
      dens[dens > 100] <- 100 # AGGRESSIVE CAP to prevent spikes
      total_dens <- total_dens + unpacked$weights[i] * dens
    }
    
    total_dens[total_dens < 1e-10] <- 1e-10
    -sum(log(total_dens))
  }

  fit <- optim(
    par = c(rep(0, K-1), unlist(initial_dist_params_list)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 3000)
  )
  
  final <- unpack(fit$par)
  result <- list(
    weights = final$weights, components = final$dist_fits,
    log_likelihood = -fit$value, n = length(data),
    distribution = "mixture", dist_names = dist_names,
    convergence = fit$convergence
  )
  class(result) <- "mixture"
  result
}

# -------------------------------
# PDF / CDF / etc
# -------------------------------
mixture_pdf <- function(x, fit) {
  out <- numeric(length(x))
  for (i in seq_along(fit$components)) {
    pdf_func <- get(paste0(fit$dist_names[i], "_pdf"))
    dens <- pdf_func(x, fit$components[[i]])
    dens[!is.finite(dens)] <- 0
    out <- out + fit$weights[i] * dens
  }
  out
}

mixture_cdf <- function(x, fit) {
  out <- numeric(length(x))
  for (i in seq_along(fit$components)) {
    cdf_func <- get(paste0(fit$dist_names[i], "_cdf"))
    probs <- cdf_func(x, fit$components[[i]])
    if (any(!is.finite(probs))) {
        med <- get(paste0(fit$dist_names[i], "_median"))(fit$components[[i]])
        probs[!is.finite(probs) & x > med] <- 1
        probs[!is.finite(probs)] <- 0
    }
    out <- out + fit$weights[i] * probs
  }
  pmax(0, pmin(1, out))
}

mixture_rand <- function(n, fit) {
  comp_indices <- sample(1:length(fit$components), n, replace = TRUE, prob = fit$weights)
  out <- numeric(n)
  for (i in seq_along(fit$components)) {
    idx <- which(comp_indices == i)
    if (length(idx) > 0) {
      rand_func <- get(paste0(fit$dist_names[i], "_rand"))
      out[idx] <- rand_func(length(idx), fit$components[[i]])
    }
  }
  out
}

mixture_mean <- function(fit) {
  out <- 0
  for (i in seq_along(fit$components)) {
    mean_func <- get(paste0(fit$dist_names[i], "_mean"))
    m <- mean_func(fit$components[[i]])
    if (is.finite(m)) out <- out + fit$weights[i] * m
  }
  out
}

# -------------------------------
# S3 logLik
# -------------------------------
logLik.mixture <- function(object, ...) {
  df <- length(object$weights) - 1
  for (comp in object$components) df <- df + (length(comp) - 1)
  structure(object$log_likelihood, df = df, nobs = object$n, class = "logLik")
}
