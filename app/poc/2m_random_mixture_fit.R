# POC: Random Mixture Generation and Fitting

# Initialize renv
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../pkg/rfit.R")
source("../pkg/plots/cdf_plot.R")

# Prevent Rplots.pdf from being created
if (!interactive()) {
  pdf(NULL)
}

set.seed(Sys.time())

# 1. Define a list of distributions for generation
gen_dists <- c(
  "lognormal_2p", "normal_2p", "weibull_2p", "pareto_2p", "gb2_4p", "dpln_4p",
  "bradford_1p", "fisk_2p", "johnsonsb_4p", "johnsonsl_3p", "johnsonsu_4p", "kappa4_4p"
)

# 2. Randomly pick two distributions
picked <- sample(gen_dists, 2, replace = TRUE)
cat(sprintf("Generating random mixture of: %s and %s\n", picked[1], picked[2]))

# 3. Randomly pick weights
w1 <- runif(1, 0.2, 0.8)
w2 <- 1 - w1
weights <- c(w1, w2)
cat(sprintf("Weights: %.2f (for %s), %.2f (for %s)\n", weights[1], picked[1], weights[2], picked[2]))

# 4. Generate random parameters for each component
# We'll use some reasonable ranges to ensure they are well-behaved
gen_params <- function(name) {
  if (name == "lognormal_2p") {
    list(mu = runif(1, 0, 2), sigma = runif(1, 0.2, 0.8), distribution = "lognormal_2p")
  } else if (name == "normal_2p") {
    list(mean = runif(1, 5, 15), sd = runif(1, 1, 3), distribution = "normal_2p")
  } else if (name == "weibull_2p") {
    list(shape = runif(1, 1, 3), scale = runif(1, 1, 10), distribution = "weibull_2p")
  } else if (name == "pareto_2p") {
    list(shape = runif(1, 2, 5), scale = runif(1, 1, 2), distribution = "pareto_2p")
  } else if (name == "fisk_2p") {
    list(scale = runif(1, 2, 8), shape = runif(1, 2, 5), distribution = "fisk_2p")
  } else if (name == "gb2_4p") {
    list(a = runif(1, 1, 3), b = runif(1, 5, 15), p = runif(1, 1, 3), q = runif(1, 1, 3), distribution = "gb2_4p")
  } else if (name == "dpln_4p") {
    list(alpha = runif(1, 1, 3), beta = runif(1, 1, 3), nu = runif(1, -1, 1), tau = runif(1, 0.5, 1.5), distribution = "dpln_4p")
  } else if (name == "bradford_1p") {
    upper <- runif(1, 10, 50)
    list(shape = runif(1, 2, 10), lower = 0, upper = upper, distribution = "bradford_1p")
  } else if (name == "johnsonsb_4p") {
    list(gamma = runif(1, -1, 1), delta = runif(1, 1, 3), xi = 0, lambda = runif(1, 10, 50), distribution = "johnsonsb_4p")
  } else if (name == "johnsonsl_3p") {
    list(gamma = runif(1, -1, 1), delta = runif(1, 1, 3), xi = 0, distribution = "johnsonsl_3p")
  } else if (name == "johnsonsu_4p") {
    list(gamma = runif(1, -1, 1), delta = runif(1, 1, 3), xi = runif(1, 0, 10), lambda = runif(1, 1, 5), distribution = "johnsonsu_4p")
  } else if (name == "kappa4_4p") {
    list(xi = runif(1, 0, 5), alpha = runif(1, 5, 15), k = runif(1, -0.5, 0.5), h = runif(1, 0.1, 2), distribution = "kappa4_4p")
  }
}

comp1 <- gen_params(picked[1])
comp2 <- gen_params(picked[2])

# 5. Create a mixture object for generation
true_mixture <- list(
  weights = weights,
  components = list(comp1, comp2),
  dist_names = picked
)
class(true_mixture) <- "mixture"

# 6. Generate sample data
n_samples <- 1000
cat(sprintf("Generating %d samples from the mixture...\n", n_samples))
sample_data <- mixture_rand(n_samples, true_mixture)

# 7. Fit all 2-component mixtures
cat("\nFitting all 2-component mixtures to generated data...\n")
all_fits <- fit_all_2m(sample_data, criterion = "AIC")

# 8. Display top results
cat("\nTop 25 Results (ordered by AIC):\n")
for (i in 1:min(25, length(all_fits))) {
  fit <- all_fits[[i]]
  cat(sprintf("%d. %s (AIC: %.2f)\n", i, fit$distribution_mixture, fit$fit_criterion_value))
}

# 9. Plot the best fit
best_fit <- all_fits[[1]]
output_dir <- "outputs"
output_file <- "2m_random_mixture_best_fit"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

plot_cdf_comparison(
  sample_data = sample_data,
  fit = best_fit,
  dist_cdf = mixture_cdf,
  output_dir = output_dir,
  output_file = output_file,
  title = sprintf("Best Fit to Random Mixture (Truth: %s+%s, Fitted: %s)", 
                  picked[1], picked[2], best_fit$distribution_mixture),
  x_label = "Value"
)

cat(sprintf("\nProcessing complete. Best fit plot saved to %s/%s.png\n", output_dir, output_file))
