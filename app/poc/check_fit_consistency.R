# Check consistency between single and mixture fits

source("../pkg/rfit.R")

# Load data
input_file <- "inputs/whole_dataset.csv"
data_raw <- read.csv(input_file)$traffic
data <- data_raw[!is.na(data_raw) & data_raw > 0]
set.seed(42)
data <- sample(data, 1000) # smaller subset for speed

# 1. Fit single Pareto
cat("Fitting single pareto_2p...\n")
fit_s <- pareto_2p_fit(data)
cat(sprintf("Single Pareto Log-Lik: %.4f\n", fit_s$log_likelihood))

# 2. Fit mixture of Pareto + something else
cat("\nFitting mixture pareto_2p + normal_2p...\n")
fit_m <- mixture_fit(data, c("pareto_2p", "normal_2p"))
cat(sprintf("Mixture Pareto+Normal Reported Log-Lik: %.4f\n", fit_m$log_likelihood))
cat(sprintf("Mixture Weights: %.4f, %.4f\n", fit_m$weights[1], fit_m$weights[2]))

# Calculate TRUE log-likelihood for the mixture
mixture_ll_true <- function(data, fit) {
  total_dens <- numeric(length(data))
  for (i in seq_along(fit$components)) {
    pdf_func <- get(paste0(fit$dist_names[i], "_pdf"))
    dens <- pdf_func(data, fit$components[[i]])
    total_dens <- total_dens + fit$weights[i] * dens
  }
  sum(log(total_dens))
}

cat(sprintf("Mixture True Log-Lik: %.4f\n", mixture_ll_true(data, fit_m)))

# 3. Check if we can get back to single Pareto LL by setting weight to 1
comp_pareto <- fit_s
comp_pareto$distribution <- "pareto_2p"
pseudo_fit <- list(
  weights = c(1, 0),
  components = list(comp_pareto, normal_2p_fit(data)),
  dist_names = c("pareto_2p", "normal_2p")
)
cat(sprintf("\nPseudo-mixture (weight 1 on Pareto) True Log-Lik: %.4f\n", 
            mixture_ll_true(data, pseudo_fit)))
