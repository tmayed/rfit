# Mixture Distribution Module - ROBUST VERSION

# -------------------------------
# Parameter Definitions
# -------------------------------
# Sourced from definitions.R
if (!exists(".DIST_REGISTRY")) {
  current_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) ".")
  base_dir <- dirname(current_file)
  if (base_dir == "." || is.na(base_dir)) {
    source("definitions.R")
  } else {
    source(file.path(base_dir, "definitions.R"))
  }
}

# -------------------------------
# Mixture Fit
# -------------------------------
#' @param fixed_weights Optional numeric vector of weights (must sum to 1)
mixture_fit <- function(data, dist_names, initial_fit = NULL, quiet = FALSE, ordered = TRUE, fixed_weights = NULL) {
  data <- data[!is.na(data) & data > 0]
  K <- length(dist_names)
  
  if (!is.null(fixed_weights)) {
    if (length(fixed_weights) != K) stop("length(fixed_weights) must equal length(dist_names)")
    if (abs(sum(fixed_weights) - 1) > 1e-6) stop("fixed_weights must sum to 1")
  }

  # 1. Initialization
  initial_dist_params_list <- list()
  param_counts <- numeric(K)
  
  if (!is.null(initial_fit) && all(initial_fit$dist_names == dist_names)) {
    if (!quiet) cat("Using warm start from initial_fit...\n")
    for (i in 1:K) {
      initial_dist_params_list[[i]] <- .DIST_REGISTRY[[dist_names[i]]]$to_internal(initial_fit$components[[i]])
      param_counts[i] <- length(initial_dist_params_list[[i]])
    }
    # Convert weights to multinomial logit scale (only if not fixed)
    if (is.null(fixed_weights)) {
      initial_weights_raw <- log(initial_fit$weights[1:(K-1)] / initial_fit$weights[K])
    }
  } else {
    data_sorted <- sort(data)
    for (i in 1:K) {
      name <- dist_names[i]
      fit_func <- get(paste0(name, "_fit"))
      # Split data into K chunks
      idx <- floor((i-1)*length(data)/K + 1) : floor(i*length(data)/K)
      fit_obj <- tryCatch(fit_func(data_sorted[idx]), error = function(e) fit_func(data))
      initial_dist_params_list[[i]] <- .DIST_REGISTRY[[name]]$to_internal(fit_obj)
      param_counts[i] <- length(initial_dist_params_list[[i]])
    }
    initial_weights_raw <- rep(0, K-1)
  }
  
  unpack <- function(params) {
    if (is.null(fixed_weights)) {
      w_raw <- c(params[1:(K-1)], 0)
      w_raw <- pmax(-15, pmin(15, w_raw))
      weights <- exp(w_raw) / sum(exp(w_raw))
      curr <- K
    } else {
      weights <- fixed_weights
      curr <- 1
    }
    
    dist_fits <- list()
    for (i in 1:K) {
        p_vec <- params[curr:(curr + param_counts[i] - 1)]
        p_vec <- pmax(-10, pmin(10, p_vec))
        dist_fits[[i]] <- .DIST_REGISTRY[[dist_names[i]]]$from_internal(p_vec)
        dist_fits[[i]]$distribution <- dist_names[i]
        curr <- curr + param_counts[i]
    }
    list(weights = weights, dist_fits = dist_fits)
  }
  
  neg_log_likelihood <- function(params) {
    unpacked <- tryCatch(unpack(params), error = function(e) NULL)
    if (is.null(unpacked)) return(1e15)
    
    total_dens <- numeric(length(data))
    comp_means <- numeric(K)
    
    for (i in 1:K) {
      pdf_func <- get(paste0(dist_names[i], "_pdf"))
      dens <- tryCatch(pdf_func(data, unpacked$dist_fits[[i]]), error = function(e) rep(0, length(data)))
      dens[!is.finite(dens)] <- 0
      total_dens <- total_dens + unpacked$weights[i] * dens
      
      if (ordered) {
        mean_func <- get(paste0(dist_names[i], "_mean"))
        comp_means[i] <- tryCatch(mean_func(unpacked$dist_fits[[i]]), error = function(e) Inf)
      }
    }
    
    total_dens[total_dens < 1e-12] <- 1e-12
    ll <- sum(log(total_dens))
    
    penalty <- 0
    if (ordered && K > 1) {
      for (i in 1:(K-1)) {
        if (is.finite(comp_means[i]) && is.finite(comp_means[i+1])) {
          if (comp_means[i] > comp_means[i+1]) {
            penalty <- penalty + 1e5 * (log(comp_means[i]) - log(comp_means[i+1]))^2
          }
        }
      }
    }
    
    penalty <- penalty + 1e-4 * sum(params^2)
    -ll + penalty
  }

  initial_par <- if (is.null(fixed_weights)) {
    c(initial_weights_raw, unlist(initial_dist_params_list))
  } else {
    unlist(initial_dist_params_list)
  }

  fit <- optim(
    par = initial_par,
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
  
  final <- unpack(fit$par)
  
  # Calculate TRUE log-likelihood without penalties
  total_dens <- numeric(length(data))
  for (i in 1:K) {
    pdf_func <- get(paste0(dist_names[i], "_pdf"))
    dens <- tryCatch(pdf_func(data, final$dist_fits[[i]]), error = function(e) rep(0, length(data)))
    dens[!is.finite(dens)] <- 0
    total_dens <- total_dens + final$weights[i] * dens
  }
  total_dens[total_dens < 1e-12] <- 1e-12
  true_ll <- sum(log(total_dens))

  result <- list(
    weights = final$weights, components = final$dist_fits,
    log_likelihood = true_ll, 
    penalized_ll = -fit$value,
    n = length(data),
    distribution = "mixture", dist_names = dist_names,
    convergence = fit$convergence,
    params_internal = fit$par
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

mixture_quantile <- function(p, fit) {
  # Numerical root finding for mixture quantile
  out <- numeric(length(p))
  for (i in seq_along(p)) {
    if (p[i] <= 0) {
      out[i] <- 0
      next
    }
    if (p[i] >= 1) {
      out[i] <- Inf # Or some large value if bounded
      next
    }
    
    target <- p[i]
    f <- function(x) mixture_cdf(x, fit) - target
    
    # Heuristic for bounds: use min/max of components' quantiles
    low <- 0
    high <- 1e6 # Default large
    
    # Try to find better bounds
    try({
      comp_qs <- sapply(seq_along(fit$components), function(j) {
        q_func <- get(paste0(fit$dist_names[j], "_quantile"))
        q_func(p[i], fit$components[[j]])
      })
      low <- min(comp_qs) * 0.1
      high <- max(comp_qs) * 10
    }, silent = TRUE)
    
    # Expand high bound if needed
    while(f(high) < 0 && high < 1e12) high <- high * 10
    
    res <- tryCatch({
      uniroot(f, lower = low, upper = high, extendInt = "yes")$root
    }, error = function(e) NA)
    out[i] <- res
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
