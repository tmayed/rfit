# POC: 3-Component Random Mixture Generation and Fitting

# Initialize renv
if (file.exists("../renv/activate.R")) {
  Sys.setenv(RENV_PROJECT = normalizePath(".."))
  source("../renv/activate.R")
}

source("../app/rfit.R")
source("../app/plots/cdf_plot.R")

set.seed(Sys.time())

# 1. Define a list of well-behaved distributions for generation
gen_dists <- c("lognormal", "weibull", "pareto", "fisk")

# 2. Randomly pick three distributions
picked <- sample(gen_dists, 3, replace = TRUE)
cat(sprintf("Generating random mixture of: %s, %s, and %s\n", picked[1], picked[2], picked[3]))

# 3. Randomly pick weights (softmax style)
w_raw <- runif(3)
weights <- exp(w_raw) / sum(exp(w_raw))
cat(sprintf("Weights: %.2f (%s), %.2f (%s), %.2f (%s)\n", 
            weights[1], picked[1], weights[2], picked[2], weights[3], picked[3]))

# 4. Generate random parameters for each component
gen_params <- function(name) {
  if (name == "lognormal") {
    list(mu = runif(1, 0, 2), sigma = runif(1, 0.2, 0.8), distribution = "lognormal")
  } else if (name == "weibull") {
    list(shape = runif(1, 1, 3), scale = runif(1, 1, 10), distribution = "weibull")
  } else if (name == "pareto") {
    list(shape = runif(1, 2, 5), scale = runif(1, 1, 2), distribution = "pareto")
  } else if (name == "fisk") {
    list(scale = runif(1, 2, 8), shape = runif(1, 2, 5), distribution = "fisk")
  }
}

comp1 <- gen_params(picked[1])
comp2 <- gen_params(picked[2])
comp3 <- gen_params(picked[3])

# 5. Create a mixture object for generation
true_mixture <- list(
  weights = weights,
  components = list(comp1, comp2, comp3),
  dist_names = picked
)
class(true_mixture) <- "mixture"

# 6. Generate sample data
n_samples <- 1500
cat(sprintf("Generating %d samples from the mixture...\n", n_samples))
sample_data <- mixture_rand(n_samples, true_mixture)

# 7. Fit all 3-component mixtures
cat("\nFitting all 3-component mixtures to generated data...\n")
all_fits <- fit_all_3m(sample_data, criterion = "AIC")

# 8. Display top results
cat("\nTop 5 Results (ordered by AIC):\n")
for (i in 1:min(5, length(all_fits))) {
  fit <- all_fits[[i]]
  cat(sprintf("%d. %s (AIC: %.2f)\n", i, fit$distribution_mixture, fit$fit_criterion_value))
}

# 9. Plot the best fit
best_fit <- all_fits[[1]]
cat("\nBest Fit Details:\n")
print(best_fit)

output_dir <- "outputs"
output_file <- "3m_random_mixture_best_fit"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

plot_cdf_comparison(
  sample_data = sample_data,
  fit = best_fit,
  dist_cdf = mixture_cdf,
  output_dir = output_dir,
  output_file = output_file,
  title = sprintf("Best Fit to 3-Mixture (Truth: %s+%s+%s, Fitted: %s)", 
                  picked[1], picked[2], picked[3], best_fit$distribution_mixture),
  x_label = "Value"
)

cat(sprintf("\nProcessing complete. Best fit plot saved to %s/%s.png\n", output_dir, output_file))
