# Mixture Distribution Module with Mean Constraint

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
# Mixture Fit with Mean Constraint
# -------------------------------
#' @param initial_fit Optional previous fit object to use as starting parameters (warm start)
#' @param quiet Logical, if TRUE suppresses terminal output
#' @param ordered Logical, if TRUE enforces mean(comp 1) < mean(comp 2) < ...
mixture_fit_mean <- function(data, dist_names, initial_fit = NULL, quiet = FALSE, ordered = TRUE) {
  data <- data[!is.na(data) & data > 0]
  empirical_mean <- mean(data)
  K <- length(dist_names)
  
  # 1. Initialization
  initial_dist_params_list <- list()
  param_counts <- numeric(K)
  
  if (!is.null(initial_fit) && all(initial_fit$dist_names == dist_names)) {
    if (!quiet) cat("Using warm start from initial_fit...\n")
    for (i in 1:K) {
      initial_dist_params_list[[i]] <- .DIST_REGISTRY[[dist_names[i]]]$to_internal(initial_fit$components[[i]])
      param_counts[i] <- length(initial_dist_params_list[[i]])
    }
    initial_weights_raw <- log(initial_fit$weights[1:(K-1)] / initial_fit$weights[K])
  } else {
    data_sorted <- sort(data)
    for (i in 1:K) {
      name <- dist_names[i]
      fit_func <- get(paste0(name, "_fit"))
      idx <- floor((i-1)*length(data)/K + 1) : floor(i*length(data)/K)
      fit_obj <- tryCatch(fit_func(data_sorted[idx]), error = function(e) fit_func(data))
      initial_dist_params_list[[i]] <- .DIST_REGISTRY[[name]]$to_internal(fit_obj)
      param_counts[i] <- length(initial_dist_params_list[[i]])
    }
    initial_weights_raw <- rep(0, K-1)
  }
  
  unpack <- function(params) {
    w_raw <- c(params[1:(K-1)], 0)
    w_raw <- pmax(-15, pmin(15, w_raw))
    weights <- exp(w_raw) / sum(exp(w_raw))
    dist_fits <- list()
    curr <- K
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
    theo_mean <- 0
    comp_means <- numeric(K)
    mean_finite <- TRUE
    
    for (i in 1:K) {
      pdf_func <- get(paste0(dist_names[i], "_pdf"))
      dens <- tryCatch(pdf_func(data, unpacked$dist_fits[[i]]), error = function(e) rep(0, length(data)))
      dens[!is.finite(dens)] <- 0
      total_dens <- total_dens + unpacked$weights[i] * dens
      
      mean_func <- get(paste0(dist_names[i], "_mean"))
      m <- tryCatch(mean_func(unpacked$dist_fits[[i]]), error = function(e) Inf)
      comp_means[i] <- m
      if (is.finite(m)) {
        theo_mean <- theo_mean + unpacked$weights[i] * m
      } else {
        mean_finite <- FALSE
      }
    }
    
    total_dens[total_dens < 1e-12] <- 1e-12
    ll <- sum(log(total_dens))
    
    penalty <- 0
    # Mean constraint penalty
    if (!mean_finite) {
      penalty <- 1e12
    } else {
      penalty <- 1e6 * (theo_mean - empirical_mean)^2
    }
    
    # Ordering penalty
    if (ordered && K > 1) {
      for (i in 1:(K-1)) {
        if (is.finite(comp_means[i]) && is.finite(comp_means[i+1])) {
          if (comp_means[i] > comp_means[i+1]) {
            penalty <- penalty + 1e5 * (log(comp_means[i]) - log(comp_means[i+1]))^2
          }
        }
      }
    }
    
    # Regularization
    penalty <- penalty + 1e-4 * sum(params^2)
    
    -ll + penalty
  }

  fit <- optim(
    par = c(initial_weights_raw, unlist(initial_dist_params_list)),
    fn = neg_log_likelihood,
    method = "Nelder-Mead",
    control = list(maxit = 5000)
  )
  
  final <- unpack(fit$par)
  result <- list(
    weights = final$weights, components = final$dist_fits,
    log_likelihood = -fit$value, n = length(data),
    distribution = "mixture", dist_names = dist_names,
    convergence = fit$convergence,
    empirical_mean = empirical_mean,
    params_internal = fit$par
  )
  class(result) <- "mixture"
  result
}
