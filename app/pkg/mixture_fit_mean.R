# Mixture Distribution Module with Mean Constraint

# -------------------------------
# Parameter Definitions
# -------------------------------
# Sourced from definitions.R
if (!exists(".dist_params")) {
  current_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) ".")
  base_dir <- dirname(current_file)
  if (base_dir == "." || is.na(base_dir)) {
    source("definitions.R")
  } else {
    source(file.path(base_dir, "definitions.R"))
  }
}

# -------------------------------
# Mixture Fit with Mean Constraint
# -------------------------------
mixture_fit_mean <- function(data, dist_names) {
  data <- data[!is.na(data) & data > 0]
  empirical_mean <- mean(data)
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
    theo_mean <- 0
    mean_finite <- TRUE
    
    for (i in 1:K) {
      # PDF calculation
      pdf_func <- get(paste0(dist_names[i], "_pdf"))
      dens <- tryCatch(pdf_func(data, unpacked$dist_fits[[i]]), error = function(e) rep(0, length(data)))
      dens[!is.finite(dens)] <- 0
      dens[dens > 100] <- 100 # AGGRESSIVE CAP to prevent spikes
      total_dens <- total_dens + unpacked$weights[i] * dens
      
      # Mean calculation for constraint
      mean_func <- get(paste0(dist_names[i], "_mean"))
      m <- mean_func(unpacked$dist_fits[[i]])
      if (is.finite(m)) {
        theo_mean <- theo_mean + unpacked$weights[i] * m
      } else {
        mean_finite <- FALSE
      }
    }
    
    total_dens[total_dens < 1e-10] <- 1e-10
    ll <- sum(log(total_dens))
    
    # Penalty term for mean constraint
    # We want Theoretical Mean = Empirical Mean
    if (!mean_finite) {
      penalty <- 1e12 # Huge penalty for infinite mean
    } else {
      # Squared error penalty. 1e6 is usually large enough to enforce constraint
      # while allowing the optimizer to move.
      penalty <- 1e6 * (theo_mean - empirical_mean)^2
    }
    
    -ll + penalty
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
    convergence = fit$convergence,
    empirical_mean = empirical_mean
  )
  class(result) <- "mixture"
  result
}

# -------------------------------
# PDF / CDF / etc (Re-used from mixture.R logic)
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
    else return(Inf)
  }
  out
}
