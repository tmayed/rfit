# rfit - Probability Distribution Fitting Package
# Main module that loads all distribution modules

#' @title rfit
#' @description A package to perform probability distribution fitting on data
#' @details This package provides functions for fitting various probability
#'          distributions using maximum likelihood estimation, including
#'          support for truncated data and weighted mixtures of distributions.
#' @export
NULL

# Source all distribution modules relative to this file
# This assumes the distribution files are in a 'dist' subdirectory relative to this file
current_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) ".")
base_dir <- dirname(current_file)

if (base_dir == ".") {
  # Fallback for interactive use or if sys.frame(1)$ofile is not available
  source("dist/lognormal.R")
  source("dist/normal.R")
  source("dist/weibull.R")
  source("dist/pareto.R")
  source("dist/gb2.R")
  source("dist/dpln.R")
  source("dist/kappa4.R")
  source("dist/bradford.R")
  source("dist/johnsonsu.R")
  source("dist/johnsonsb.R")
  source("dist/johnsonsl.R")
  source("dist/fisk.R")
  source("fit.R")
  source("fit_2m.R")
  source("fit_3m.R")
  source("mixture.R")
  source("plots/cdf_plot.R")
  source("plots/pdf_plot.R")
} else {
  source(file.path(base_dir, "dist/lognormal.R"))
  source(file.path(base_dir, "dist/normal.R"))
  source(file.path(base_dir, "dist/weibull.R"))
  source(file.path(base_dir, "dist/pareto.R"))
  source(file.path(base_dir, "dist/gb2.R"))
  source(file.path(base_dir, "dist/dpln.R"))
  source(file.path(base_dir, "dist/kappa4.R"))
  source(file.path(base_dir, "dist/bradford.R"))
  source(file.path(base_dir, "dist/johnsonsu.R"))
  source(file.path(base_dir, "dist/johnsonsb.R"))
  source(file.path(base_dir, "dist/johnsonsl.R"))
  source(file.path(base_dir, "dist/fisk.R"))
  source(file.path(base_dir, "fit.R"))
  source(file.path(base_dir, "fit_2m.R"))
  source(file.path(base_dir, "fit_3m.R"))
  source(file.path(base_dir, "mixture.R"))
  source(file.path(base_dir, "plots/cdf_plot.R"))
  source(file.path(base_dir, "plots/pdf_plot.R"))
}

# Plotting Exports
#' @export
plot_cdf_comparison <- plot_cdf_comparison
#' @export
plot_pdf_comparison <- plot_pdf_comparison

# Orchestrator Exports
#' @export
fit_all <- fit_all
#' @export
fit_all_2m <- fit_all_2m
#' @export
fit_all_3m <- fit_all_3m

# Fisk Exports
#' @export
fisk_log_likelihood <- fisk_log_likelihood
#' @export
fisk_fit <- fisk_fit
#' @export
fisk_fit_truncated <- fisk_fit_truncated
#' @export
fisk_pdf <- fisk_pdf
#' @export
fisk_cdf <- fisk_cdf
#' @export
fisk_quantile <- fisk_quantile
#' @export
fisk_rand <- fisk_rand
#' @export
fisk_mean <- fisk_mean
#' @export
fisk_std <- fisk_std
#' @export
fisk_sf <- fisk_sf
#' @export
fisk_isf <- fisk_isf
#' @export
fisk_logpdf <- fisk_logpdf
#' @export
fisk_logcdf <- fisk_logcdf
#' @export
fisk_logsf <- fisk_logsf
#' @export
fisk_var <- fisk_var
#' @export
fisk_moment <- fisk_moment
#' @export
fisk_skew <- fisk_skew
#' @export
fisk_kurtosis <- fisk_kurtosis
#' @export
fisk_median <- fisk_median
#' @export
fisk_interval <- fisk_interval
#' @export
fisk_entropy <- fisk_entropy
#' @export
fisk_expect <- fisk_expect

# Normal Exports
#' @export
normal_log_likelihood <- normal_log_likelihood
#' @export
normal_fit <- normal_fit
#' @export
normal_fit_truncated <- normal_fit_truncated
#' @export
normal_pdf <- normal_pdf
#' @export
normal_cdf <- normal_cdf
#' @export
normal_quantile <- normal_quantile
#' @export
normal_rand <- normal_rand
#' @export
normal_mean <- normal_mean
#' @export
normal_std <- normal_std
#' @export
normal_sf <- normal_sf
#' @export
normal_isf <- normal_isf
#' @export
normal_logpdf <- normal_logpdf
#' @export
normal_logcdf <- normal_logcdf
#' @export
normal_logsf <- normal_logsf
#' @export
normal_var <- normal_var
#' @export
normal_moment <- normal_moment
#' @export
normal_skew <- normal_skew
#' @export
normal_kurtosis <- normal_kurtosis
#' @export
normal_median <- normal_median
#' @export
normal_interval <- normal_interval
#' @export
normal_entropy <- normal_entropy
#' @export
normal_expect <- normal_expect

# Weibull Exports
#' @export
weibull_log_likelihood <- weibull_log_likelihood
#' @export
weibull_fit <- weibull_fit
#' @export
weibull_fit_truncated <- weibull_fit_truncated
#' @export
weibull_pdf <- weibull_pdf
#' @export
weibull_cdf <- weibull_cdf
#' @export
weibull_quantile <- weibull_quantile
#' @export
weibull_rand <- weibull_rand
#' @export
weibull_mean <- weibull_mean
#' @export
weibull_std <- weibull_std
#' @export
weibull_sf <- weibull_sf
#' @export
weibull_isf <- weibull_isf
#' @export
weibull_logpdf <- weibull_logpdf
#' @export
weibull_logcdf <- weibull_logcdf
#' @export
weibull_logsf <- weibull_logsf
#' @export
weibull_var <- weibull_var
#' @export
weibull_moment <- weibull_moment
#' @export
weibull_skew <- weibull_skew
#' @export
weibull_kurtosis <- weibull_kurtosis
#' @export
weibull_median <- weibull_median
#' @export
weibull_interval <- weibull_interval
#' @export
weibull_entropy <- weibull_entropy
#' @export
weibull_expect <- weibull_expect

# Kappa 4 Exports
#' @export
kappa4_log_likelihood <- kappa4_log_likelihood
#' @export
kappa4_fit <- kappa4_fit
#' @export
kappa4_fit_truncated <- kappa4_fit_truncated
#' @export
kappa4_pdf <- kappa4_pdf
#' @export
kappa4_cdf <- kappa4_cdf
#' @export
kappa4_quantile <- kappa4_quantile
#' @export
kappa4_rand <- kappa4_rand
#' @export
kappa4_mean <- kappa4_mean
#' @export
kappa4_std <- kappa4_std
#' @export
kappa4_sf <- kappa4_sf
#' @export
kappa4_isf <- kappa4_isf
#' @export
kappa4_logpdf <- kappa4_logpdf
#' @export
kappa4_logcdf <- kappa4_logcdf
#' @export
kappa4_logsf <- kappa4_logsf
#' @export
kappa4_var <- kappa4_var
#' @export
kappa4_moment <- kappa4_moment
#' @export
kappa4_skew <- kappa4_skew
#' @export
kappa4_kurtosis <- kappa4_kurtosis
#' @export
kappa4_median <- kappa4_median
#' @export
kappa4_interval <- kappa4_interval
#' @export
kappa4_entropy <- kappa4_entropy
#' @export
kappa4_expect <- kappa4_expect

# GB2 Exports
#' @export
gb2_log_likelihood <- gb2_log_likelihood
#' @export
gb2_fit <- gb2_fit
#' @export
gb2_fit_truncated <- gb2_fit_truncated
#' @export
gb2_pdf <- gb2_pdf
#' @export
gb2_cdf <- gb2_cdf
#' @export
gb2_quantile <- gb2_quantile
#' @export
gb2_rand <- gb2_rand
#' @export
gb2_mean <- gb2_mean
#' @export
gb2_std <- gb2_std
#' @export
gb2_sf <- gb2_sf
#' @export
gb2_isf <- gb2_isf
#' @export
gb2_logpdf <- gb2_logpdf
#' @export
gb2_logcdf <- gb2_logcdf
#' @export
gb2_logsf <- gb2_logsf
#' @export
gb2_var <- gb2_var
#' @export
gb2_moment <- gb2_moment
#' @export
gb2_skew <- gb2_skew
#' @export
gb2_kurtosis <- gb2_kurtosis
#' @export
gb2_median <- gb2_median
#' @export
gb2_interval <- gb2_interval
#' @export
gb2_entropy <- gb2_entropy
#' @export
gb2_expect <- gb2_expect

# Lognormal Exports
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
lognormal_mean <- lognormal_mean
#' @export
lognormal_std <- lognormal_std
#' @export
lognormal_sf <- lognormal_sf
#' @export
lognormal_isf <- lognormal_isf
#' @export
lognormal_logpdf <- lognormal_logpdf
#' @export
lognormal_logcdf <- lognormal_logcdf
#' @export
lognormal_logsf <- lognormal_logsf
#' @export
lognormal_var <- lognormal_var
#' @export
lognormal_moment <- lognormal_moment
#' @export
lognormal_skew <- lognormal_skew
#' @export
lognormal_kurtosis <- lognormal_kurtosis
#' @export
lognormal_median <- lognormal_median
#' @export
lognormal_interval <- lognormal_interval
#' @export
lognormal_entropy <- lognormal_entropy
#' @export
lognormal_expect <- lognormal_expect

# Pareto Exports
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
#' @export
pareto_mean <- pareto_mean
#' @export
pareto_std <- pareto_std
#' @export
pareto_sf <- pareto_sf
#' @export
pareto_isf <- pareto_isf
#' @export
pareto_logpdf <- pareto_logpdf
#' @export
pareto_logcdf <- pareto_logcdf
#' @export
pareto_logsf <- pareto_logsf
#' @export
pareto_var <- pareto_var
#' @export
pareto_moment <- pareto_moment
#' @export
pareto_skew <- pareto_skew
#' @export
pareto_kurtosis <- pareto_kurtosis
#' @export
pareto_median <- pareto_median
#' @export
pareto_interval <- pareto_interval
#' @export
pareto_entropy <- pareto_entropy
#' @export
pareto_expect <- pareto_expect

# dPLN Exports
#' @export
dpln_log_likelihood <- dpln_log_likelihood
#' @export
dpln_fit <- dpln_fit
#' @export
dpln_fit_truncated <- dpln_fit_truncated
#' @export
dpln_pdf <- dpln_pdf
#' @export
dpln_cdf <- dpln_cdf
#' @export
dpln_quantile <- dpln_quantile
#' @export
dpln_rand <- dpln_rand
#' @export
dpln_mean <- dpln_mean
#' @export
dpln_std <- dpln_std
#' @export
dpln_sf <- dpln_sf
#' @export
dpln_isf <- dpln_isf
#' @export
dpln_logpdf <- dpln_logpdf
#' @export
dpln_logcdf <- dpln_logcdf
#' @export
dpln_logsf <- dpln_logsf
#' @export
dpln_var <- dpln_var
#' @export
dpln_moment <- dpln_moment
#' @export
dpln_skew <- dpln_skew
#' @export
dpln_kurtosis <- dpln_kurtosis
#' @export
dpln_median <- dpln_median
#' @export
dpln_interval <- dpln_interval
#' @export
dpln_entropy <- dpln_entropy
#' @export
dpln_expect <- dpln_expect

# Bradford Exports
#' @export
bradford_log_likelihood <- bradford_log_likelihood
#' @export
bradford_fit <- bradford_fit
#' @export
bradford_fit_truncated <- bradford_fit_truncated
#' @export
bradford_pdf <- bradford_pdf
#' @export
bradford_cdf <- bradford_cdf
#' @export
bradford_quantile <- bradford_quantile
#' @export
bradford_rand <- bradford_rand
#' @export
bradford_mean <- bradford_mean
#' @export
bradford_std <- bradford_std
#' @export
bradford_sf <- bradford_sf
#' @export
bradford_isf <- bradford_isf
#' @export
bradford_logpdf <- bradford_logpdf
#' @export
bradford_logcdf <- bradford_logcdf
#' @export
bradford_logsf <- bradford_logsf
#' @export
bradford_var <- bradford_var
#' @export
bradford_moment <- bradford_moment
#' @export
bradford_skew <- bradford_skew
#' @export
bradford_kurtosis <- bradford_kurtosis
#' @export
bradford_median <- bradford_median
#' @export
bradford_interval <- bradford_interval
#' @export
bradford_entropy <- bradford_entropy
#' @export
bradford_expect <- bradford_expect

# Johnson SU Exports
#' @export
johnsonsu_log_likelihood <- johnsonsu_log_likelihood
#' @export
johnsonsu_fit <- johnsonsu_fit
#' @export
johnsonsu_fit_truncated <- johnsonsu_fit_truncated
#' @export
johnsonsu_pdf <- johnsonsu_pdf
#' @export
johnsonsu_cdf <- johnsonsu_cdf
#' @export
johnsonsu_quantile <- johnsonsu_quantile
#' @export
johnsonsu_rand <- johnsonsu_rand
#' @export
johnsonsu_mean <- johnsonsu_mean
#' @export
johnsonsu_std <- johnsonsu_std
#' @export
johnsonsu_sf <- johnsonsu_sf
#' @export
johnsonsu_isf <- johnsonsu_isf
#' @export
johnsonsu_logpdf <- johnsonsu_logpdf
#' @export
johnsonsu_logcdf <- johnsonsu_logcdf
#' @export
johnsonsu_logsf <- johnsonsu_logsf
#' @export
johnsonsu_var <- johnsonsu_var
#' @export
johnsonsu_moment <- johnsonsu_moment
#' @export
johnsonsu_skew <- johnsonsu_skew
#' @export
johnsonsu_kurtosis <- johnsonsu_kurtosis
#' @export
johnsonsu_median <- johnsonsu_median
#' @export
johnsonsu_interval <- johnsonsu_interval
#' @export
johnsonsu_entropy <- johnsonsu_entropy
#' @export
johnsonsu_expect <- johnsonsu_expect

# Johnson SB Exports
#' @export
johnsonsb_log_likelihood <- johnsonsb_log_likelihood
#' @export
johnsonsb_fit <- johnsonsb_fit
#' @export
johnsonsb_fit_truncated <- johnsonsb_fit_truncated
#' @export
johnsonsb_pdf <- johnsonsb_pdf
#' @export
johnsonsb_cdf <- johnsonsb_cdf
#' @export
johnsonsb_quantile <- johnsonsb_quantile
#' @export
johnsonsb_rand <- johnsonsb_rand
#' @export
johnsonsb_mean <- johnsonsb_mean
#' @export
johnsonsb_std <- johnsonsb_std
#' @export
johnsonsb_sf <- johnsonsb_sf
#' @export
johnsonsb_isf <- johnsonsb_isf
#' @export
johnsonsb_logpdf <- johnsonsb_logpdf
#' @export
johnsonsb_logcdf <- johnsonsb_logcdf
#' @export
johnsonsb_logsf <- johnsonsb_logsf
#' @export
johnsonsb_var <- johnsonsb_var
#' @export
johnsonsb_moment <- johnsonsb_moment
#' @export
johnsonsb_skew <- johnsonsb_skew
#' @export
johnsonsb_kurtosis <- johnsonsb_kurtosis
#' @export
johnsonsb_median <- johnsonsb_median
#' @export
johnsonsb_interval <- johnsonsb_interval
#' @export
johnsonsb_entropy <- johnsonsb_entropy
#' @export
johnsonsb_expect <- johnsonsb_expect

# Johnson SL Exports
#' @export
johnsonsl_log_likelihood <- johnsonsl_log_likelihood
#' @export
johnsonsl_fit <- johnsonsl_fit
#' @export
johnsonsl_fit_truncated <- johnsonsl_fit_truncated
#' @export
johnsonsl_pdf <- johnsonsl_pdf
#' @export
johnsonsl_cdf <- johnsonsl_cdf
#' @export
johnsonsl_quantile <- johnsonsl_quantile
#' @export
johnsonsl_rand <- johnsonsl_rand
#' @export
johnsonsl_mean <- johnsonsl_mean
#' @export
johnsonsl_std <- johnsonsl_std
#' @export
johnsonsl_sf <- johnsonsl_sf
#' @export
johnsonsl_isf <- johnsonsl_isf
#' @export
johnsonsl_logpdf <- johnsonsl_logpdf
#' @export
johnsonsl_logcdf <- johnsonsl_logcdf
#' @export
johnsonsl_logsf <- johnsonsl_logsf
#' @export
johnsonsl_var <- johnsonsl_var
#' @export
johnsonsl_moment <- johnsonsl_moment
#' @export
johnsonsl_skew <- johnsonsl_skew
#' @export
johnsonsl_kurtosis <- johnsonsl_kurtosis
#' @export
johnsonsl_median <- johnsonsl_median
#' @export
johnsonsl_interval <- johnsonsl_interval
#' @export
johnsonsl_entropy <- johnsonsl_entropy
#' @export
johnsonsl_expect <- johnsonsl_expect

# Mixture Exports
#' @export
mixture_fit <- mixture_fit
#' @export
mixture_pdf <- mixture_pdf
#' @export
mixture_cdf <- mixture_cdf
#' @export
mixture_rand <- mixture_rand
#' @export
mixture_mean <- mixture_mean
