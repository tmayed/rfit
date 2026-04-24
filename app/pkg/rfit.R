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

if (base_dir == "." || is.na(base_dir)) {
  # Fallback for interactive use
  source("dist/lognormal_2p.R")
  source("dist/normal_2p.R")
  source("dist/weibull_2p.R")
  source("dist/pareto_2p.R")
  source("dist/gb2_4p.R")
  source("dist/dpln_4p.R")
  source("dist/kappa4_4p.R")
  source("dist/bradford_1p.R")
  source("dist/johnsonsu_4p.R")
  source("dist/johnsonsb_4p.R")
  source("dist/johnsonsl_3p.R")
  source("dist/fisk_2p.R")
  source("dist/rayleigh_1p.R")
  source("dist/rayleigh_2p.R")
  source("dist/betarayleigh_4p.R")
  source("dist/levy_2p.R")
  source("dist/nakagami_2p.R")
  source("dist/gig_3p.R")
  source("dist/gamma_2p.R")
  source("definitions.R")
  source("fit.R")
  source("fit_2m.R")
  source("fit_3m.R")
  source("fit_ln_gamma_dpnl.R")
  source("fit_ln_gamma.R")
  source("mixture.R")
  source("mixture_fit_mean.R")
  source("segment_data.R")
  source("plots/cdf_plot.R")
  source("plots/pdf_plot.R")
  source("plots/diag_plot.R")
} else {
  source(file.path(base_dir, "dist/lognormal_2p.R"))
  source(file.path(base_dir, "dist/normal_2p.R"))
  source(file.path(base_dir, "dist/weibull_2p.R"))
  source(file.path(base_dir, "dist/pareto_2p.R"))
  source(file.path(base_dir, "dist/gb2_4p.R"))
  source(file.path(base_dir, "dist/dpln_4p.R"))
  source(file.path(base_dir, "dist/kappa4_4p.R"))
  source(file.path(base_dir, "dist/bradford_1p.R"))
  source(file.path(base_dir, "dist/johnsonsu_4p.R"))
  source(file.path(base_dir, "dist/johnsonsb_4p.R"))
  source(file.path(base_dir, "dist/johnsonsl_3p.R"))
  source(file.path(base_dir, "dist/fisk_2p.R"))
  source(file.path(base_dir, "dist/rayleigh_1p.R"))
  source(file.path(base_dir, "dist/rayleigh_2p.R"))
  source(file.path(base_dir, "dist/birnbaumsaunders_3p.R"))
  source(file.path(base_dir, "dist/betarayleigh_4p.R"))
  source(file.path(base_dir, "dist/levy_2p.R"))
  source(file.path(base_dir, "dist/nakagami_2p.R"))
  source(file.path(base_dir, "dist/gig_3p.R"))
  source(file.path(base_dir, "dist/gamma_2p.R"))
  source(file.path(base_dir, "definitions.R"))
  source(file.path(base_dir, "fit.R"))
  source(file.path(base_dir, "fit_2m.R"))
  source(file.path(base_dir, "fit_3m.R"))
  source(file.path(base_dir, "fit_ln_gamma_dpnl.R"))
  source(file.path(base_dir, "fit_ln_gamma.R"))
  source(file.path(base_dir, "mixture.R"))
  source(file.path(base_dir, "mixture_fit_mean.R"))
  source(file.path(base_dir, "segment_data.R"))
  source(file.path(base_dir, "plots/cdf_plot.R"))
  source(file.path(base_dir, "plots/pdf_plot.R"))
  source(file.path(base_dir, "plots/diag_plot.R"))
}

# Mixture Exports
#' @export
mixture_fit <- mixture_fit
#' @export
mixture_fit_mean <- mixture_fit_mean
#' @export
fit_ln_gamma_dpnl <- fit_ln_gamma_dpnl
#' @export
fit_ln_gamma <- fit_ln_gamma
#' @export
mixture_pdf <- mixture_pdf
#' @export
mixture_cdf <- mixture_cdf
#' @export
mixture_quantile <- mixture_quantile
#' @export
mixture_rand <- mixture_rand
#' @export
mixture_mean <- mixture_mean

# Plotting Exports
#' @export
plot_cdf_comparison <- plot_cdf_comparison
#' @export
plot_pdf_comparison <- plot_pdf_comparison
#' @export
plot_diag <- plot_diag
#' @export
segment_data <- segment_data
#' @export
logLik.ln_gamma_dpnl_mixture <- logLik.ln_gamma_dpnl_mixture

# Birnbaum-Saunders Exports
#' @export
birnbaumsaunders_3p_log_likelihood <- birnbaumsaunders_3p_log_likelihood
#' @export
birnbaumsaunders_3p_fit <- birnbaumsaunders_3p_fit
#' @export
birnbaumsaunders_3p_fit_truncated <- birnbaumsaunders_3p_fit_truncated
#' @export
birnbaumsaunders_3p_pdf <- birnbaumsaunders_3p_pdf
#' @export
birnbaumsaunders_3p_cdf <- birnbaumsaunders_3p_cdf
#' @export
birnbaumsaunders_3p_quantile <- birnbaumsaunders_3p_quantile
#' @export
birnbaumsaunders_3p_rand <- birnbaumsaunders_3p_rand
#' @export
birnbaumsaunders_3p_mean <- birnbaumsaunders_3p_mean
#' @export
birnbaumsaunders_3p_std <- birnbaumsaunders_3p_std
#' @export
birnbaumsaunders_3p_sf <- birnbaumsaunders_3p_sf
#' @export
birnbaumsaunders_3p_isf <- birnbaumsaunders_3p_isf
#' @export
birnbaumsaunders_3p_logpdf <- birnbaumsaunders_3p_logpdf
#' @export
birnbaumsaunders_3p_logcdf <- birnbaumsaunders_3p_logcdf
#' @export
birnbaumsaunders_3p_logsf <- birnbaumsaunders_3p_logsf
#' @export
birnbaumsaunders_3p_var <- birnbaumsaunders_3p_var
#' @export
birnbaumsaunders_3p_moment <- birnbaumsaunders_3p_moment
#' @export
birnbaumsaunders_3p_skew <- birnbaumsaunders_3p_skew
#' @export
birnbaumsaunders_3p_kurtosis <- birnbaumsaunders_3p_kurtosis
#' @export
birnbaumsaunders_3p_median <- birnbaumsaunders_3p_median
#' @export
birnbaumsaunders_3p_interval <- birnbaumsaunders_3p_interval
#' @export
birnbaumsaunders_3p_entropy <- birnbaumsaunders_3p_entropy
#' @export
birnbaumsaunders_3p_expect <- birnbaumsaunders_3p_expect

# GIG Exports
#' @export
gig_3p_log_likelihood <- gig_3p_log_likelihood
#' @export
gig_3p_fit <- gig_3p_fit
#' @export
gig_3p_fit_truncated <- gig_3p_fit_truncated
#' @export
gig_3p_pdf <- gig_3p_pdf
#' @export
gig_3p_cdf <- gig_3p_cdf
#' @export
gig_3p_quantile <- gig_3p_quantile
#' @export
gig_3p_rand <- gig_3p_rand
#' @export
gig_3p_mean <- gig_3p_mean
#' @export
gig_3p_std <- gig_3p_std
#' @export
gig_3p_sf <- gig_3p_sf
#' @export
gig_3p_isf <- gig_3p_isf
#' @export
gig_3p_logpdf <- gig_3p_logpdf
#' @export
gig_3p_logcdf <- gig_3p_logcdf
#' @export
gig_3p_logsf <- gig_3p_logsf
#' @export
gig_3p_var <- gig_3p_var
#' @export
gig_3p_moment <- gig_3p_moment
#' @export
gig_3p_skew <- gig_3p_skew
#' @export
gig_3p_kurtosis <- gig_3p_kurtosis
#' @export
gig_3p_median <- gig_3p_median
#' @export
gig_3p_interval <- gig_3p_interval
#' @export
gig_3p_entropy <- gig_3p_entropy
#' @export
gig_3p_expect <- gig_3p_expect

# Gamma Exports
#' @export
gamma_2p_log_likelihood <- gamma_2p_log_likelihood
#' @export
gamma_2p_fit <- gamma_2p_fit
#' @export
gamma_2p_fit_truncated <- gamma_2p_fit_truncated
#' @export
gamma_2p_pdf <- gamma_2p_pdf
#' @export
gamma_2p_cdf <- gamma_2p_cdf
#' @export
gamma_2p_quantile <- gamma_2p_quantile
#' @export
gamma_2p_rand <- gamma_2p_rand
#' @export
gamma_2p_mean <- gamma_2p_mean
#' @export
gamma_2p_std <- gamma_2p_std
#' @export
gamma_2p_sf <- gamma_2p_sf
#' @export
gamma_2p_isf <- gamma_2p_isf
#' @export
gamma_2p_logpdf <- gamma_2p_logpdf
#' @export
gamma_2p_logcdf <- gamma_2p_logcdf
#' @export
gamma_2p_logsf <- gamma_2p_logsf
#' @export
gamma_2p_var <- gamma_2p_var
#' @export
gamma_2p_moment <- gamma_2p_moment
#' @export
gamma_2p_skew <- gamma_2p_skew
#' @export
gamma_2p_kurtosis <- gamma_2p_kurtosis
#' @export
gamma_2p_median <- gamma_2p_median
#' @export
gamma_2p_interval <- gamma_2p_interval
#' @export
gamma_2p_entropy <- gamma_2p_entropy
#' @export
gamma_2p_expect <- gamma_2p_expect

# Nakagami Exports
#' @export
nakagami_2p_log_likelihood <- nakagami_2p_log_likelihood
#' @export
nakagami_2p_fit <- nakagami_2p_fit
#' @export
nakagami_2p_fit_truncated <- nakagami_2p_fit_truncated
#' @export
nakagami_2p_pdf <- nakagami_2p_pdf
#' @export
nakagami_2p_cdf <- nakagami_2p_cdf
#' @export
nakagami_2p_quantile <- nakagami_2p_quantile
#' @export
nakagami_2p_rand <- nakagami_2p_rand
#' @export
nakagami_2p_mean <- nakagami_2p_mean
#' @export
nakagami_2p_std <- nakagami_2p_std
#' @export
nakagami_2p_sf <- nakagami_2p_sf
#' @export
nakagami_2p_isf <- nakagami_2p_isf
#' @export
nakagami_2p_logpdf <- nakagami_2p_logpdf
#' @export
nakagami_2p_logcdf <- nakagami_2p_logcdf
#' @export
nakagami_2p_logsf <- nakagami_2p_logsf
#' @export
nakagami_2p_var <- nakagami_2p_var
#' @export
nakagami_2p_moment <- nakagami_2p_moment
#' @export
nakagami_2p_skew <- nakagami_2p_skew
#' @export
nakagami_2p_kurtosis <- nakagami_2p_kurtosis
#' @export
nakagami_2p_median <- nakagami_2p_median
#' @export
nakagami_2p_interval <- nakagami_2p_interval
#' @export
nakagami_2p_entropy <- nakagami_2p_entropy
#' @export
nakagami_2p_expect <- nakagami_2p_expect

# Levy Exports
#' @export
levy_2p_log_likelihood <- levy_2p_log_likelihood
#' @export
levy_2p_fit <- levy_2p_fit
#' @export
levy_2p_fit_truncated <- levy_2p_fit_truncated
#' @export
levy_2p_pdf <- levy_2p_pdf
#' @export
levy_2p_cdf <- levy_2p_cdf
#' @export
levy_2p_quantile <- levy_2p_quantile
#' @export
levy_2p_rand <- levy_2p_rand
#' @export
levy_2p_mean <- levy_2p_mean
#' @export
levy_2p_std <- levy_2p_std
#' @export
levy_2p_sf <- levy_2p_sf
#' @export
levy_2p_isf <- levy_2p_isf
#' @export
levy_2p_logpdf <- levy_2p_logpdf
#' @export
levy_2p_logcdf <- levy_2p_logcdf
#' @export
levy_2p_logsf <- levy_2p_logsf
#' @export
levy_2p_var <- levy_2p_var
#' @export
levy_2p_moment <- levy_2p_moment
#' @export
levy_2p_skew <- levy_2p_skew
#' @export
levy_2p_kurtosis <- levy_2p_kurtosis
#' @export
levy_2p_median <- levy_2p_median
#' @export
levy_2p_interval <- levy_2p_interval
#' @export
levy_2p_entropy <- levy_2p_entropy
#' @export
levy_2p_expect <- levy_2p_expect

# Beta-Rayleigh Exports
#' @export
betarayleigh_4p_log_likelihood <- betarayleigh_4p_log_likelihood
#' @export
betarayleigh_4p_fit <- betarayleigh_4p_fit
#' @export
betarayleigh_4p_fit_truncated <- betarayleigh_4p_fit_truncated
#' @export
betarayleigh_4p_pdf <- betarayleigh_4p_pdf
#' @export
betarayleigh_4p_cdf <- betarayleigh_4p_cdf
#' @export
betarayleigh_4p_quantile <- betarayleigh_4p_quantile
#' @export
betarayleigh_4p_rand <- betarayleigh_4p_rand
#' @export
betarayleigh_4p_mean <- betarayleigh_4p_mean
#' @export
betarayleigh_4p_std <- betarayleigh_4p_std
#' @export
betarayleigh_4p_sf <- betarayleigh_4p_sf
#' @export
betarayleigh_4p_isf <- betarayleigh_4p_isf
#' @export
betarayleigh_4p_logpdf <- betarayleigh_4p_logpdf
#' @export
betarayleigh_4p_logcdf <- betarayleigh_4p_logcdf
#' @export
betarayleigh_4p_logsf <- betarayleigh_4p_logsf
#' @export
betarayleigh_4p_var <- betarayleigh_4p_var
#' @export
betarayleigh_4p_moment <- betarayleigh_4p_moment
#' @export
betarayleigh_4p_skew <- betarayleigh_4p_skew
#' @export
betarayleigh_4p_kurtosis <- betarayleigh_4p_kurtosis
#' @export
betarayleigh_4p_median <- betarayleigh_4p_median
#' @export
betarayleigh_4p_interval <- betarayleigh_4p_interval
#' @export
betarayleigh_4p_entropy <- betarayleigh_4p_entropy
#' @export
betarayleigh_4p_expect <- betarayleigh_4p_expect

# Rayleigh 1P Exports
#' @export
rayleigh_1p_log_likelihood <- rayleigh_1p_log_likelihood
#' @export
rayleigh_1p_fit <- rayleigh_1p_fit
#' @export
rayleigh_1p_fit_truncated <- rayleigh_1p_fit_truncated
#' @export
rayleigh_1p_pdf <- rayleigh_1p_pdf
#' @export
rayleigh_1p_cdf <- rayleigh_1p_cdf
#' @export
rayleigh_1p_quantile <- rayleigh_1p_quantile
#' @export
rayleigh_1p_rand <- rayleigh_1p_rand
#' @export
rayleigh_1p_mean <- rayleigh_1p_mean
#' @export
rayleigh_1p_std <- rayleigh_1p_std
#' @export
rayleigh_1p_sf <- rayleigh_1p_sf
#' @export
rayleigh_1p_isf <- rayleigh_1p_isf
#' @export
rayleigh_1p_logpdf <- rayleigh_1p_logpdf
#' @export
rayleigh_1p_logcdf <- rayleigh_1p_logcdf
#' @export
rayleigh_1p_logsf <- rayleigh_1p_logsf
#' @export
rayleigh_1p_var <- rayleigh_1p_var
#' @export
rayleigh_1p_moment <- rayleigh_1p_moment
#' @export
rayleigh_1p_skew <- rayleigh_1p_skew
#' @export
rayleigh_1p_kurtosis <- rayleigh_1p_kurtosis
#' @export
rayleigh_1p_median <- rayleigh_1p_median
#' @export
rayleigh_1p_interval <- rayleigh_1p_interval
#' @export
rayleigh_1p_entropy <- rayleigh_1p_entropy
#' @export
rayleigh_1p_expect <- rayleigh_1p_expect

# Rayleigh 2P Exports
#' @export
rayleigh_2p_log_likelihood <- rayleigh_2p_log_likelihood
#' @export
rayleigh_2p_fit <- rayleigh_2p_fit
#' @export
rayleigh_2p_fit_truncated <- rayleigh_2p_fit_truncated
#' @export
rayleigh_2p_pdf <- rayleigh_2p_pdf
#' @export
rayleigh_2p_cdf <- rayleigh_2p_cdf
#' @export
rayleigh_2p_quantile <- rayleigh_2p_quantile
#' @export
rayleigh_2p_rand <- rayleigh_2p_rand
#' @export
rayleigh_2p_mean <- rayleigh_2p_mean
#' @export
rayleigh_2p_std <- rayleigh_2p_std
#' @export
rayleigh_2p_sf <- rayleigh_2p_sf
#' @export
rayleigh_2p_isf <- rayleigh_2p_isf
#' @export
rayleigh_2p_logpdf <- rayleigh_2p_logpdf
#' @export
rayleigh_2p_logcdf <- rayleigh_2p_logcdf
#' @export
rayleigh_2p_logsf <- rayleigh_2p_logsf
#' @export
rayleigh_2p_var <- rayleigh_2p_var
#' @export
rayleigh_2p_moment <- rayleigh_2p_moment
#' @export
rayleigh_2p_skew <- rayleigh_2p_skew
#' @export
rayleigh_2p_kurtosis <- rayleigh_2p_kurtosis
#' @export
rayleigh_2p_median <- rayleigh_2p_median
#' @export
rayleigh_2p_interval <- rayleigh_2p_interval
#' @export
rayleigh_2p_entropy <- rayleigh_2p_entropy
#' @export
rayleigh_2p_expect <- rayleigh_2p_expect

# Orchestrator Exports
#' @export
fit_all <- fit_all
#' @export
fit_all_2m <- fit_all_2m
#' @export
fit_all_3m <- fit_all_3m

# Fisk Exports
#' @export
fisk_2p_log_likelihood <- fisk_2p_log_likelihood
#' @export
fisk_2p_fit <- fisk_2p_fit
#' @export
fisk_2p_fit_truncated <- fisk_2p_fit_truncated
#' @export
fisk_2p_pdf <- fisk_2p_pdf
#' @export
fisk_2p_cdf <- fisk_2p_cdf
#' @export
fisk_2p_quantile <- fisk_2p_quantile
#' @export
fisk_2p_rand <- fisk_2p_rand
#' @export
fisk_2p_mean <- fisk_2p_mean
#' @export
fisk_2p_std <- fisk_2p_std
#' @export
fisk_2p_sf <- fisk_2p_sf
#' @export
fisk_2p_isf <- fisk_2p_isf
#' @export
fisk_2p_logpdf <- fisk_2p_logpdf
#' @export
fisk_2p_logcdf <- fisk_2p_logcdf
#' @export
fisk_2p_logsf <- fisk_2p_logsf
#' @export
fisk_2p_var <- fisk_2p_var
#' @export
fisk_2p_moment <- fisk_2p_moment
#' @export
fisk_2p_skew <- fisk_2p_skew
#' @export
fisk_2p_kurtosis <- fisk_2p_kurtosis
#' @export
fisk_2p_median <- fisk_2p_median
#' @export
fisk_2p_interval <- fisk_2p_interval
#' @export
fisk_2p_entropy <- fisk_2p_entropy
#' @export
fisk_2p_expect <- fisk_2p_expect

# Normal Exports
#' @export
normal_2p_log_likelihood <- normal_2p_log_likelihood
#' @export
normal_2p_fit <- normal_2p_fit
#' @export
normal_2p_fit_truncated <- normal_2p_fit_truncated
#' @export
normal_2p_pdf <- normal_2p_pdf
#' @export
normal_2p_cdf <- normal_2p_cdf
#' @export
normal_2p_quantile <- normal_2p_quantile
#' @export
normal_2p_rand <- normal_2p_rand
#' @export
normal_2p_mean <- normal_2p_mean
#' @export
normal_2p_std <- normal_2p_std
#' @export
normal_2p_sf <- normal_2p_sf
#' @export
normal_2p_isf <- normal_2p_isf
#' @export
normal_2p_logpdf <- normal_2p_logpdf
#' @export
normal_2p_logcdf <- normal_2p_logcdf
#' @export
normal_2p_logsf <- normal_2p_logsf
#' @export
normal_2p_var <- normal_2p_var
#' @export
normal_2p_moment <- normal_2p_moment
#' @export
normal_2p_skew <- normal_2p_skew
#' @export
normal_2p_kurtosis <- normal_2p_kurtosis
#' @export
normal_2p_median <- normal_2p_median
#' @export
normal_2p_interval <- normal_2p_interval
#' @export
normal_2p_entropy <- normal_2p_entropy
#' @export
normal_2p_expect <- normal_2p_expect

# Weibull Exports
#' @export
weibull_2p_log_likelihood <- weibull_2p_log_likelihood
#' @export
weibull_2p_fit <- weibull_2p_fit
#' @export
weibull_2p_fit_truncated <- weibull_2p_fit_truncated
#' @export
weibull_2p_pdf <- weibull_2p_pdf
#' @export
weibull_2p_cdf <- weibull_2p_cdf
#' @export
weibull_2p_quantile <- weibull_2p_quantile
#' @export
weibull_2p_rand <- weibull_2p_rand
#' @export
weibull_2p_mean <- weibull_2p_mean
#' @export
weibull_2p_std <- weibull_2p_std
#' @export
weibull_2p_sf <- weibull_2p_sf
#' @export
weibull_2p_isf <- weibull_2p_isf
#' @export
weibull_2p_logpdf <- weibull_2p_logpdf
#' @export
weibull_2p_logcdf <- weibull_2p_logcdf
#' @export
weibull_2p_logsf <- weibull_2p_logsf
#' @export
weibull_2p_var <- weibull_2p_var
#' @export
weibull_2p_moment <- weibull_2p_moment
#' @export
weibull_2p_skew <- weibull_2p_skew
#' @export
weibull_2p_kurtosis <- weibull_2p_kurtosis
#' @export
weibull_2p_median <- weibull_2p_median
#' @export
weibull_2p_interval <- weibull_2p_interval
#' @export
weibull_2p_entropy <- weibull_2p_entropy
#' @export
weibull_2p_expect <- weibull_2p_expect

# Kappa 4 Exports
#' @export
kappa4_4p_log_likelihood <- kappa4_4p_log_likelihood
#' @export
kappa4_4p_fit <- kappa4_4p_fit
#' @export
kappa4_4p_fit_truncated <- kappa4_4p_fit_truncated
#' @export
kappa4_4p_pdf <- kappa4_4p_pdf
#' @export
kappa4_4p_cdf <- kappa4_4p_cdf
#' @export
kappa4_4p_quantile <- kappa4_4p_quantile
#' @export
kappa4_4p_rand <- kappa4_4p_rand
#' @export
kappa4_4p_mean <- kappa4_4p_mean
#' @export
kappa4_4p_std <- kappa4_4p_std
#' @export
kappa4_4p_sf <- kappa4_4p_sf
#' @export
kappa4_4p_isf <- kappa4_4p_isf
#' @export
kappa4_4p_logpdf <- kappa4_4p_logpdf
#' @export
kappa4_4p_logcdf <- kappa4_4p_logcdf
#' @export
kappa4_4p_logsf <- kappa4_4p_logsf
#' @export
kappa4_4p_var <- kappa4_4p_var
#' @export
kappa4_4p_moment <- kappa4_4p_moment
#' @export
kappa4_4p_skew <- kappa4_4p_skew
#' @export
kappa4_4p_kurtosis <- kappa4_4p_kurtosis
#' @export
kappa4_4p_median <- kappa4_4p_median
#' @export
kappa4_4p_interval <- kappa4_4p_interval
#' @export
kappa4_4p_entropy <- kappa4_4p_entropy
#' @export
kappa4_4p_expect <- kappa4_4p_expect

# GB2 Exports
#' @export
gb2_4p_log_likelihood <- gb2_4p_log_likelihood
#' @export
gb2_4p_fit <- gb2_4p_fit
#' @export
gb2_4p_fit_truncated <- gb2_4p_fit_truncated
#' @export
gb2_4p_pdf <- gb2_4p_pdf
#' @export
gb2_4p_cdf <- gb2_4p_cdf
#' @export
gb2_4p_quantile <- gb2_4p_quantile
#' @export
gb2_4p_rand <- gb2_4p_rand
#' @export
gb2_4p_mean <- gb2_4p_mean
#' @export
gb2_4p_std <- gb2_4p_std
#' @export
gb2_4p_sf <- gb2_4p_sf
#' @export
gb2_4p_isf <- gb2_4p_isf
#' @export
gb2_4p_logpdf <- gb2_4p_logpdf
#' @export
gb2_4p_logcdf <- gb2_4p_logcdf
#' @export
gb2_4p_logsf <- gb2_4p_logsf
#' @export
gb2_4p_var <- gb2_4p_var
#' @export
gb2_4p_moment <- gb2_4p_moment
#' @export
gb2_4p_skew <- gb2_4p_skew
#' @export
gb2_4p_kurtosis <- gb2_4p_kurtosis
#' @export
gb2_4p_median <- gb2_4p_median
#' @export
gb2_4p_interval <- gb2_4p_interval
#' @export
gb2_4p_entropy <- gb2_4p_entropy
#' @export
gb2_4p_expect <- gb2_4p_expect

# Lognormal Exports
#' @export
lognormal_2p_log_likelihood <- lognormal_2p_log_likelihood
#' @export
lognormal_2p_fit <- lognormal_2p_fit
#' @export
lognormal_2p_fit_truncated <- lognormal_2p_fit_truncated
#' @export
lognormal_2p_pdf <- lognormal_2p_pdf
#' @export
lognormal_2p_cdf <- lognormal_2p_cdf
#' @export
lognormal_2p_quantile <- lognormal_2p_quantile
#' @export
lognormal_2p_rand <- lognormal_2p_rand
#' @export
lognormal_2p_mean <- lognormal_2p_mean
#' @export
lognormal_2p_std <- lognormal_2p_std
#' @export
lognormal_2p_sf <- lognormal_2p_sf
#' @export
lognormal_2p_isf <- lognormal_2p_isf
#' @export
lognormal_2p_logpdf <- lognormal_2p_logpdf
#' @export
lognormal_2p_logcdf <- lognormal_2p_logcdf
#' @export
lognormal_2p_logsf <- lognormal_2p_logsf
#' @export
lognormal_2p_var <- lognormal_2p_var
#' @export
lognormal_2p_moment <- lognormal_2p_moment
#' @export
lognormal_2p_skew <- lognormal_2p_skew
#' @export
lognormal_2p_kurtosis <- lognormal_2p_kurtosis
#' @export
lognormal_2p_median <- lognormal_2p_median
#' @export
lognormal_2p_interval <- lognormal_2p_interval
#' @export
lognormal_2p_entropy <- lognormal_2p_entropy
#' @export
lognormal_2p_expect <- lognormal_2p_expect

# Pareto Exports
#' @export
pareto_2p_log_likelihood <- pareto_2p_log_likelihood
#' @export
pareto_2p_fit <- pareto_2p_fit
#' @export
pareto_2p_fit_truncated <- pareto_2p_fit_truncated
#' @export
pareto_2p_pdf <- pareto_2p_pdf
#' @export
pareto_2p_cdf <- pareto_2p_cdf
#' @export
pareto_2p_quantile <- pareto_2p_quantile
#' @export
pareto_2p_rand <- pareto_2p_rand
#' @export
pareto_2p_mean <- pareto_2p_mean
#' @export
pareto_2p_std <- pareto_2p_std
#' @export
pareto_2p_sf <- pareto_2p_sf
#' @export
pareto_2p_isf <- pareto_2p_isf
#' @export
pareto_2p_logpdf <- pareto_2p_logpdf
#' @export
pareto_2p_logcdf <- pareto_2p_logcdf
#' @export
pareto_2p_logsf <- pareto_2p_logsf
#' @export
pareto_2p_var <- pareto_2p_var
#' @export
pareto_2p_moment <- pareto_2p_moment
#' @export
pareto_2p_skew <- pareto_2p_skew
#' @export
pareto_2p_kurtosis <- pareto_2p_kurtosis
#' @export
pareto_2p_median <- pareto_2p_median
#' @export
pareto_2p_interval <- pareto_2p_interval
#' @export
pareto_2p_entropy <- pareto_2p_entropy
#' @export
pareto_2p_expect <- pareto_2p_expect

# dPLN Exports
#' @export
dpln_4p_log_likelihood <- dpln_4p_log_likelihood
#' @export
dpln_4p_fit <- dpln_4p_fit
#' @export
dpln_4p_fit_truncated <- dpln_4p_fit_truncated
#' @export
dpln_4p_pdf <- dpln_4p_pdf
#' @export
dpln_4p_cdf <- dpln_4p_cdf
#' @export
dpln_4p_quantile <- dpln_4p_quantile
#' @export
dpln_4p_rand <- dpln_4p_rand
#' @export
dpln_4p_mean <- dpln_4p_mean
#' @export
dpln_4p_std <- dpln_4p_std
#' @export
dpln_4p_sf <- dpln_4p_sf
#' @export
dpln_4p_isf <- dpln_4p_isf
#' @export
dpln_4p_logpdf <- dpln_4p_logpdf
#' @export
dpln_4p_logcdf <- dpln_4p_logcdf
#' @export
dpln_4p_logsf <- dpln_4p_logsf
#' @export
dpln_4p_var <- dpln_4p_var
#' @export
dpln_4p_moment <- dpln_4p_moment
#' @export
dpln_4p_skew <- dpln_4p_skew
#' @export
dpln_4p_kurtosis <- dpln_4p_kurtosis
#' @export
dpln_4p_median <- dpln_4p_median
#' @export
dpln_4p_interval <- dpln_4p_interval
#' @export
dpln_4p_entropy <- dpln_4p_entropy
#' @export
dpln_4p_expect <- dpln_4p_expect

# Bradford Exports
#' @export
bradford_1p_log_likelihood <- bradford_1p_log_likelihood
#' @export
bradford_1p_fit <- bradford_1p_fit
#' @export
bradford_1p_fit_truncated <- bradford_1p_fit_truncated
#' @export
bradford_1p_pdf <- bradford_1p_pdf
#' @export
bradford_1p_cdf <- bradford_1p_cdf
#' @export
bradford_1p_quantile <- bradford_1p_quantile
#' @export
bradford_1p_rand <- bradford_1p_rand
#' @export
bradford_1p_mean <- bradford_1p_mean
#' @export
bradford_1p_std <- bradford_1p_std
#' @export
bradford_1p_sf <- bradford_1p_sf
#' @export
bradford_1p_isf <- bradford_1p_isf
#' @export
bradford_1p_logpdf <- bradford_1p_logpdf
#' @export
bradford_1p_logcdf <- bradford_1p_logcdf
#' @export
bradford_1p_logsf <- bradford_1p_logsf
#' @export
bradford_1p_var <- bradford_1p_var
#' @export
bradford_1p_moment <- bradford_1p_moment
#' @export
bradford_1p_skew <- bradford_1p_skew
#' @export
bradford_1p_kurtosis <- bradford_1p_kurtosis
#' @export
bradford_1p_median <- bradford_1p_median
#' @export
bradford_1p_interval <- bradford_1p_interval
#' @export
bradford_1p_entropy <- bradford_1p_entropy
#' @export
bradford_1p_expect <- bradford_1p_expect

# Johnson SU Exports
#' @export
johnsonsu_4p_log_likelihood <- johnsonsu_4p_log_likelihood
#' @export
johnsonsu_4p_fit <- johnsonsu_4p_fit
#' @export
johnsonsu_4p_fit_truncated <- johnsonsu_4p_fit_truncated
#' @export
johnsonsu_4p_pdf <- johnsonsu_4p_pdf
#' @export
johnsonsu_4p_cdf <- johnsonsu_4p_cdf
#' @export
johnsonsu_4p_quantile <- johnsonsu_4p_quantile
#' @export
johnsonsu_4p_rand <- johnsonsu_4p_rand
#' @export
johnsonsu_4p_mean <- johnsonsu_4p_mean
#' @export
johnsonsu_4p_std <- johnsonsu_4p_std
#' @export
johnsonsu_4p_sf <- johnsonsu_4p_sf
#' @export
johnsonsu_4p_isf <- johnsonsu_4p_isf
#' @export
johnsonsu_4p_logpdf <- johnsonsu_4p_logpdf
#' @export
johnsonsu_4p_logcdf <- johnsonsu_4p_logcdf
#' @export
johnsonsu_4p_logsf <- johnsonsu_4p_logsf
#' @export
johnsonsu_4p_var <- johnsonsu_4p_var
#' @export
johnsonsu_4p_moment <- johnsonsu_4p_moment
#' @export
johnsonsu_4p_skew <- johnsonsu_4p_skew
#' @export
johnsonsu_4p_kurtosis <- johnsonsu_4p_kurtosis
#' @export
johnsonsu_4p_median <- johnsonsu_4p_median
#' @export
johnsonsu_4p_interval <- johnsonsu_4p_interval
#' @export
johnsonsu_4p_entropy <- johnsonsu_4p_entropy
#' @export
johnsonsu_4p_expect <- johnsonsu_4p_expect

# Johnson SB Exports
#' @export
johnsonsb_4p_log_likelihood <- johnsonsb_4p_log_likelihood
#' @export
johnsonsb_4p_fit <- johnsonsb_4p_fit
#' @export
johnsonsb_4p_fit_truncated <- johnsonsb_4p_fit_truncated
#' @export
johnsonsb_4p_pdf <- johnsonsb_4p_pdf
#' @export
johnsonsb_4p_cdf <- johnsonsb_4p_cdf
#' @export
johnsonsb_4p_quantile <- johnsonsb_4p_quantile
#' @export
johnsonsb_4p_rand <- johnsonsb_4p_rand
#' @export
johnsonsb_4p_mean <- johnsonsb_4p_mean
#' @export
johnsonsb_4p_std <- johnsonsb_4p_std
#' @export
johnsonsb_4p_sf <- johnsonsb_4p_sf
#' @export
johnsonsb_4p_isf <- johnsonsb_4p_isf
#' @export
johnsonsb_4p_logpdf <- johnsonsb_4p_logpdf
#' @export
johnsonsb_4p_logcdf <- johnsonsb_4p_logcdf
#' @export
johnsonsb_4p_logsf <- johnsonsb_4p_logsf
#' @export
johnsonsb_4p_var <- johnsonsb_4p_var
#' @export
johnsonsb_4p_moment <- johnsonsb_4p_moment
#' @export
johnsonsb_4p_skew <- johnsonsb_4p_skew
#' @export
johnsonsb_4p_kurtosis <- johnsonsb_4p_kurtosis
#' @export
johnsonsb_4p_median <- johnsonsb_4p_median
#' @export
johnsonsb_4p_interval <- johnsonsb_4p_interval
#' @export
johnsonsb_4p_entropy <- johnsonsb_4p_entropy
#' @export
johnsonsb_4p_expect <- johnsonsb_4p_expect

# Johnson SL Exports
#' @export
johnsonsl_3p_log_likelihood <- johnsonsl_3p_log_likelihood
#' @export
johnsonsl_3p_fit <- johnsonsl_3p_fit
#' @export
johnsonsl_3p_fit_truncated <- johnsonsl_3p_fit_truncated
#' @export
johnsonsl_3p_pdf <- johnsonsl_3p_pdf
#' @export
johnsonsl_3p_cdf <- johnsonsl_3p_cdf
#' @export
johnsonsl_3p_quantile <- johnsonsl_3p_quantile
#' @export
johnsonsl_3p_rand <- johnsonsl_3p_rand
#' @export
johnsonsl_3p_mean <- johnsonsl_3p_mean
#' @export
johnsonsl_3p_std <- johnsonsl_3p_std
#' @export
johnsonsl_3p_sf <- johnsonsl_3p_sf
#' @export
johnsonsl_3p_isf <- johnsonsl_3p_isf
#' @export
johnsonsl_3p_logpdf <- johnsonsl_3p_logpdf
#' @export
johnsonsl_3p_logcdf <- johnsonsl_3p_logcdf
#' @export
johnsonsl_3p_logsf <- johnsonsl_3p_logsf
#' @export
johnsonsl_3p_var <- johnsonsl_3p_var
#' @export
johnsonsl_3p_moment <- johnsonsl_3p_moment
#' @export
johnsonsl_3p_skew <- johnsonsl_3p_skew
#' @export
johnsonsl_3p_kurtosis <- johnsonsl_3p_kurtosis
#' @export
johnsonsl_3p_median <- johnsonsl_3p_median
#' @export
johnsonsl_3p_interval <- johnsonsl_3p_interval
#' @export
johnsonsl_3p_entropy <- johnsonsl_3p_entropy
#' @export
johnsonsl_3p_expect <- johnsonsl_3p_expect
 johnsonsl_3p_expect
