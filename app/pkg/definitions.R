# rfit - Distribution Parameter Definitions

# -------------------------------
# Filtering Utility
# -------------------------------
#' @title filter_distributions
#' @description Returns a subset of available distributions matching specified criteria.
#' @param np Number of parameters (match exactly).
#' @param lb Lower bound of the distribution's support.
#' @param ub Upper bound of the distribution's support.
#' @param f_0 Boolean; if TRUE, filters for distributions where f(0) = 0 is possible.
#' @param f_0_strict Boolean; if TRUE, filters for distributions where f(0) = 0 is guaranteed.
#' @param mean_def Boolean; if TRUE, filters for distributions with well-defined means for all parameters.
#' @return A named list containing distribution metadata for matches.
#' @export
filter_distributions <- function(np = NULL, lb = NULL, ub = NULL, 
                                 f_0 = NULL, f_0_strict = NULL, 
                                 mean_def = NULL) {
  Filter(function(x) {
    keep <- TRUE
    if (!is.null(np) && x$np != np) keep <- FALSE
    if (!is.null(lb) && x$domain[1] != lb) keep <- FALSE
    if (!is.null(ub) && x$domain[2] != ub) keep <- FALSE
    if (!is.null(f_0) && x$f_0 != f_0) keep <- FALSE
    if (!is.null(f_0_strict) && x$f_0_strict != f_0_strict) keep <- FALSE
    if (!is.null(mean_def) && x$mean_def != mean_def) keep <- FALSE
    return(keep)
  }, .DIST_REGISTRY)
}

# -------------------------------
# Parameter Definitions
# -------------------------------
# This list defines the parameter names and transformation functions
# used by mixture fitting and other cross-distribution utilities.
# 'np' is the number of parameters fitted by the model.
# 'domain' is the support of the distribution c(lower, upper).
# 'f_0' is TRUE if f(0) = 0 under SOME valid parameter regimes.
# 'f_0_strict' is TRUE if f(0) = 0 under ALL valid parameter regimes.
# 'mean_def' is TRUE if the mean is well-defined for all valid parameter choices.
.DIST_REGISTRY <- list(
  lognormal_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = TRUE,
    mean_def = TRUE,
    names = c("mu", "sigma"),
    to_internal = function(p) c(p$mu, log(p$sigma)),
    from_internal = function(v) list(mu = v[1], sigma = exp(v[2]))
  ),
  normal_2p = list(
    np = 2,
    domain = c(-Inf, Inf),
    f_0 = FALSE,
    f_0_strict = FALSE,
    mean_def = TRUE,
    names = c("mu", "sigma"),
    to_internal = function(p) c(p$mu, log(p$sigma)),
    from_internal = function(v) list(mu = v[1], sigma = exp(v[2]))
  ),
  weibull_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires shape > 1
    mean_def = TRUE,
    names = c("shape", "scale"),
    to_internal = function(p) c(log(p$shape), log(p$scale)),
    from_internal = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  pareto_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = TRUE,
    mean_def = FALSE, # shape > 1
    names = c("shape", "scale"),
    to_internal = function(p) c(log(p$shape), log(p$scale)),
    from_internal = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  dpln_4p = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires beta > 1
    mean_def = FALSE, # alpha > 1
    names = c("alpha", "beta", "nu", "tau"),
    to_internal = function(p) c(log(p$alpha), log(p$beta), p$nu, log(p$tau)),
    from_internal = function(v) list(alpha = exp(v[1]), beta = exp(v[2]), nu = v[3], tau = exp(v[4]))
  ),
  gb2_4p = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires a*p > 1
    mean_def = FALSE, # a * q > 1
    names = c("a", "b", "p", "q"),
    to_internal = function(p) c(log(p$a), log(p$b), log(p$p), log(p$q)),
    from_internal = function(v) list(a = exp(v[1]), b = exp(v[2]), p = exp(v[3]), q = exp(v[4]))
  ),
  kappa4_4p = list(
    np = 4,
    domain = c(-Inf, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires xi > 0
    mean_def = FALSE, # k > -1
    names = c("xi", "alpha", "k", "h"),
    to_internal = function(p) c(p$xi, log(p$alpha), p$k, log(p$h)),
    from_internal = function(v) list(xi = v[1], alpha = exp(v[2]), k = v[3], h = exp(v[4]))
  ),
  bradford_1p = list(
    np = 1,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires lower > 0
    mean_def = TRUE,
    names = c("shape", "lower", "upper"),
    to_internal = function(p) c(log(p$shape), p$lower, p$upper),
    from_internal = function(v) list(shape = exp(v[1]), lower = v[2], upper = v[3])
  ),
  fisk_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires shape > 1
    mean_def = FALSE, # shape > 1
    names = c("scale", "shape"),
    to_internal = function(p) c(log(p$scale), log(p$shape)),
    from_internal = function(v) list(scale = exp(v[1]), shape = exp(v[2]))
  ),
  johnsonsb_4p = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires xi > 0
    mean_def = TRUE,
    names = c("gamma", "delta", "xi", "lambda"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi, log(p$lambda)),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3], lambda = exp(v[4]))
  ),
  johnsonsl_3p = list(
    np = 3,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires xi >= 0
    mean_def = TRUE,
    names = c("gamma", "delta", "xi"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3])
  ),
  johnsonsu_4p = list(
    np = 4,
    domain = c(-Inf, Inf),
    f_0 = FALSE,
    f_0_strict = FALSE,
    mean_def = TRUE,
    names = c("gamma", "delta", "xi", "lambda"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi, log(p$lambda)),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3], lambda = exp(v[4]))
  ),
  rayleigh_1p = list(
    np = 1,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires mu >= 0 (here mu is fixed at 0)
    mean_def = TRUE,
    names = c("sigma"),
    to_internal = function(p) c(log(p$sigma)),
    from_internal = function(v) list(sigma = exp(v[1]))
  ),
  rayleigh_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires mu >= 0
    mean_def = TRUE,
    names = c("mu", "sigma"),
    to_internal = function(p) c(p$mu, log(p$sigma)),
    from_internal = function(v) list(mu = v[1], sigma = exp(v[2]))
  ),
  betarayleigh_4p = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires mu >= 0
    mean_def = TRUE,
    names = c("a", "b", "mu", "sigma"),
    to_internal = function(p) c(log(p$a), log(p$b), p$mu, log(p$sigma)),
    from_internal = function(v) list(a = exp(v[1]), b = exp(v[2]), mu = v[3], sigma = exp(v[4]))
  ),
  levy_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires mu >= 0
    mean_def = FALSE, # no mean function
    names = c("mu", "sigma"),
    to_internal = function(p) c(p$mu, log(p$sigma)),
    from_internal = function(v) list(mu = v[1], sigma = exp(v[2]))
  ),
  nakagami_2p = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = FALSE, # Requires m > 0.5
    mean_def = TRUE,
    names = c("m", "omega"),
    to_internal = function(p) c(log(p$m), log(p$omega)),
    from_internal = function(v) list(m = exp(v[1]), omega = exp(v[2]))
  ),
  birnbaumsaunders_3p = list(
    np = 3,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = TRUE,
    mean_def = TRUE,
    names = c("alpha", "beta", "mu"),
    to_internal = function(p) c(log(p$alpha), log(p$beta), p$mu),
    from_internal = function(v) list(alpha = exp(v[1]), beta = exp(v[2]), mu = v[3])
  ),
  gig_3p = list(
    np = 3,
    domain = c(0, Inf),
    f_0 = TRUE,
    f_0_strict = TRUE,
    mean_def = TRUE,
    names = c("lambda", "chi", "psi"),
    to_internal = function(p) c(p$lambda, log(p$chi), log(p$psi)),
    from_internal = function(v) list(lambda = v[1], chi = exp(v[2]), psi = exp(v[3]))
  )
)
