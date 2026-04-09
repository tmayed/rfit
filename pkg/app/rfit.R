# rfit - Probability Distribution Fitting Package
# Main module that loads all distribution modules

#' @title rfit
#' @description A package to perform probability distribution fitting on data
#' @details This package provides functions for fitting various probability
#'          distributions using maximum likelihood estimation, including
#'          support for truncated data and weighted mixtures of distributions.
#' @export
NULL

# Source all distribution modules
source("dist/lognormal.R")
source("dist/pareto.R")

#' @export
lognormal_log_likelihood <- lognormal_log_likelihood

#' @export
lognormal_fit <- lognormal_fit

#' @export
lognormal_fit_truncated <- lognormal_fit_truncated

#' @export
lognormal_pdf <- lognormal_pdf

#' @export
lognormal_cdf <- lognormal_cdf

#' @export
lognormal_quantile <- lognormal_quantile

#' @export
lognormal_rand <- lognormal_rand

#' @export
pareto_log_likelihood <- pareto_log_likelihood

#' @export
pareto_fit <- pareto_fit

#' @export
pareto_fit_truncated <- pareto_fit_truncated

#' @export
pareto_pdf <- pareto_pdf

#' @export
pareto_cdf <- pareto_cdf

#' @export
pareto_quantile <- pareto_quantile

#' @export
pareto_rand <- pareto_rand
