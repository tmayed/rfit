# rfit - Distribution Parameter Definitions

# -------------------------------
# Parameter Definitions
# -------------------------------
# This list defines the parameter names and transformation functions
# used by mixture fitting and other cross-distribution utilities.
# 'np' is the number of parameters fitted by the model.
# 'domain' is the support of the distribution c(lower, upper).
# 'f_0' is TRUE if the density f(0) = 0 for ALL valid parameters, FALSE otherwise.
.dist_params <- list(
  lognormal = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE,
    names = c("mu", "sigma"),
    to_internal = function(p) c(p$mu, log(p$sigma)),
    from_internal = function(v) list(mu = v[1], sigma = exp(v[2]))
  ),
  normal = list(
    np = 2,
    domain = c(-Inf, Inf),
    f_0 = FALSE,
    names = c("mean", "sd"),
    to_internal = function(p) c(p$mean, log(p$sd)),
    from_internal = function(v) list(mean = v[1], sd = exp(v[2]))
  ),
  weibull = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = FALSE, # f(0) > 0 if shape <= 1
    names = c("shape", "scale"),
    to_internal = function(p) c(log(p$shape), log(p$scale)),
    from_internal = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  pareto = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = TRUE, # Support is [scale, Inf) with scale > 0
    names = c("shape", "scale"),
    to_internal = function(p) c(log(p$shape), log(p$scale)),
    from_internal = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  dpln = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = TRUE,
    names = c("alpha", "beta", "nu", "tau"),
    to_internal = function(p) c(log(p$alpha), log(p$beta), p$nu, log(p$tau)),
    from_internal = function(v) list(alpha = exp(v[1]), beta = exp(v[2]), nu = v[3], tau = exp(v[4]))
  ),
  gb2 = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = FALSE, # f(0) can be non-zero depending on parameters
    names = c("a", "b", "p", "q"),
    to_internal = function(p) c(log(p$a), log(p$b), log(p$p), log(p$q)),
    from_internal = function(v) list(a = exp(v[1]), b = exp(v[2]), p = exp(v[3]), q = exp(v[4]))
  ),
  kappa4 = list(
    np = 4,
    domain = c(-Inf, Inf),
    f_0 = FALSE,
    names = c("xi", "alpha", "k", "h"),
    to_internal = function(p) c(p$xi, log(p$alpha), p$k, log(p$h)),
    from_internal = function(v) list(xi = v[1], alpha = exp(v[2]), k = v[3], h = exp(v[4]))
  ),
  bradford = list(
    np = 1,
    domain = c(0, Inf),
    f_0 = FALSE,
    names = c("shape", "lower", "upper"),
    to_internal = function(p) c(log(p$shape), p$lower, p$upper),
    from_internal = function(v) list(shape = exp(v[1]), lower = v[2], upper = v[3])
  ),
  fisk = list(
    np = 2,
    domain = c(0, Inf),
    f_0 = FALSE, # f(0) > 0 if shape <= 1
    names = c("scale", "shape"),
    to_internal = function(p) c(log(p$scale), log(p$shape)),
    from_internal = function(v) list(scale = exp(v[1]), shape = exp(v[2]))
  ),
  johnsonsb = list(
    np = 4,
    domain = c(0, Inf),
    f_0 = FALSE, # Depends on xi and shapes
    names = c("gamma", "delta", "xi", "lambda"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi, log(p$lambda)),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3], lambda = exp(v[4]))
  ),
  johnsonsl = list(
    np = 3,
    domain = c(0, Inf),
    f_0 = FALSE, # f(0) > 0 if xi < 0
    names = c("gamma", "delta", "xi"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3])
  ),
  johnsonsu = list(
    np = 4,
    domain = c(-Inf, Inf),
    f_0 = FALSE,
    names = c("gamma", "delta", "xi", "lambda"),
    to_internal = function(p) c(p$gamma, log(p$delta), p$xi, log(p$lambda)),
    from_internal = function(v) list(gamma = v[1], delta = exp(v[2]), xi = v[3], lambda = exp(v[4]))
  )
)
