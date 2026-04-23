# rfit Package

## Overview

rfit is an R package for probability distribution fitting using maximum likelihood estimation (MLE).

## Architecture

- Each distribution is stored in a separate file under `pkg/dist/`
- All functions for a distribution reside in its dedicated file
- **Implementation Standard:** Distributions MUST be implemented using standard R functions only. External libraries (e.g., `fitdistrplus`, `actuar`) are NOT permitted for distribution-specific logic (PDF, CDF, fitting, etc.).
- The package runs in a Docker container for reproducible execution

## Current Distributions

- Lognormal (`pkg/dist/lognormal_2p.R`)
- Normal (`pkg/dist/normal_2p.R`)
- Weibull (`pkg/dist/weibull_2p.R`)
- Pareto (`pkg/dist/pareto_2p.R`)
- GB2 (`pkg/dist/gb2_4p.R`)
- dPLN (`pkg/dist/dpln_4p.R`)
- Bradford (`pkg/dist/bradford_1p.R`)
- Fisk (`pkg/dist/fisk_2p.R`)
- Johnson SB (`pkg/dist/johnsonsb_4p.R`)
- Johnson SL (`pkg/dist/johnsonsl_3p.R`)
- Johnson SU (`pkg/dist/johnsonsu_4p.R`)
- Kappa 4 (`pkg/dist/kappa4_4p.R`)
- Rayleigh 1P (`pkg/dist/rayleigh_1p.R`)
- Rayleigh (`pkg/dist/rayleigh_2p.R`)
- Birnbaum-Saunders (`pkg/dist/birnbaumsaunders_3p.R`)
- Beta-Rayleigh (`pkg/dist/betarayleigh_4p.R`)
- Lévy (`pkg/dist/levy_2p.R`)
- Nakagami (`pkg/dist/nakagami_2p.R`)
- Generalized Inverse Gaussian (`pkg/dist/gig_3p.R`)
- Mixture (`pkg/mixture.R`)

## Function Naming Convention

Each distribution implements the following functions using the pattern `{dist}_{suffix}`:

| Function | Description |
|----------|-------------|
| `{dist}_log_likelihood` | Compute log-likelihood for given parameters |
| `{dist}_fit` | Fit distribution using MLE |
| `{dist}_fit_truncated` | Fit truncated distribution using MLE |
| `{dist}_pdf` | Probability density function |
| `{dist}_cdf` | Cumulative distribution function |
| `{dist}_quantile` | Quantile function |
| `{dist}_rand` | Generate random samples |
| `{dist}_mean` | Compute mean |
| `{dist}_std` | Compute standard deviation |
| `{dist}_sf` | Survival function (1 − CDF), useful for tail probabilities |
| `{dist}_isf` | Inverse survival function (inverse of SF) |
| `{dist}_logpdf` | Log of the probability density function |
| `{dist}_logcdf` | Log of the cumulative distribution function |
| `{dist}_logsf` | Log of the survival function |
| `{dist}_var` | Compute variance of the distribution |
| `{dist}_moment` | Compute non-central moment of specified order |
| `{dist}_skew` | Compute skewness (third standardized moment) |
| `{dist}_kurtosis` | Compute kurtosis (fourth standardized moment) |
| `{dist}_median` | Compute median of the distribution |
| `{dist}_interval` | Compute central interval with given confidence level |
| `{dist}_entropy` | Compute (differential) entropy of the distribution |
| `{dist}_expect` | Compute expected value of a function under the distribution |

## Definitions

The `pkg/definitions.R` file serves as the central registry for all probability distributions in the package. It defines the metadata and transformation logic required for high-level orchestrators (like `fit_all`, `mixture_fit`) to handle distributions dynamically and consistently.

### Distribution Registry (`.DIST_REGISTRY`)

Each distribution is defined by an entry in the `.DIST_REGISTRY` list with the following fields:

| Field | Description |
|-------|-------------|
| `np` | **Number of Parameters:** The count of parameters fitted by the model (e.g., 2 for Normal, 4 for dPLN). |
| `domain` | **Support Bounds:** A numeric vector `c(lower, upper)` defining the theoretical range of the distribution. |
| `f_0` | **Zero Density at Zero (Conditional):** A boolean indicating if the PDF $f(0) = 0$ is possible under certain parameter regimes (e.g., Weibull with shape > 1). |
| `f_0_strict` | **Zero Density at Zero (Strict):** A boolean indicating if the PDF $f(0) = 0$ is guaranteed for *all* valid parameter regimes (e.g., Lognormal). If `FALSE`, a `#` comment must be provided next to the value explaining the condition where it is true (e.g., `# shape > 1`). |
| `mean_def` | **Mean Existence:** A boolean indicating if the distribution's mean is well-defined for all valid parameter choices. If `FALSE`, a `#` comment must be provided next to the value explaining the condition where it is true (e.g., `# alpha > 1`) or if it never exists (e.g., `# no mean function`). |
| `names` | **Parameter Names:** A character vector of the human-readable names for the distribution parameters (e.g., `c("mu", "sigma")`). |
| `to_internal` | **Parameter Transformation:** A function that maps human-readable parameters (in a list) to the "unconstrained" internal scale used by the optimizer (typically log-transform for positive parameters). |
| `from_internal` | **Inverse Transformation:** A function that maps the optimizer's unconstrained values back to the human-readable parameter list. |

### Filtering Utility

The `filter_distributions()` function allows for querying the registry based on statistical or architectural constraints.

**Parameters:**
- `np`: Filter by exact number of parameters.
- `lb` / `ub`: Filter by lower or upper bounds of the support domain.
- `f_0` / `f_0_strict`: Filter by density properties at zero.
- `mean_def`: Filter for distributions with guaranteed finite means.

## File Structure

```
pkg/
├── dist/
├── plots/
├── mixture.R
├── rfit.R
├── fit.R
├── fit_2m.R
└── fit_3m.R
tests/
└── dist/
poc/
```

# Environment Execution Rules
This codebase runs inside an isolated, containerized R environment. 
You cannot run R code or tests directly on the local machine. 

To test code, execute scripts, or inspect the environment, you MUST use the global `run-env` command. This securely routes your commands into the execution container.

**Key Execution Context:**
* **Working Directory:** Commands run via `run-env` execute from the `/workspace` root inside the container.
* **Path Resolution:** If a script uses relative paths (e.g., R `source("../../path/to/file")`), you MUST `cd` into the script's directory within the `run-env` command to ensure paths resolve correctly.

**Examples:**
* To run an R script with relative sources: `run-env "cd pkg/tests/dist && Rscript test_file.R"`
* To run a simple command: `run-env Rscript my_file.R`
* To install a new package via renv: `run-env R -e 'renv::install("package_name")'`

When I ask you to "test the code" or "run tests", first search for existing tests in the codebase. If relevant tests exist, run them via `run-env`. If tests do not exist or are insufficient for the current task, write or update them as needed. After ensuring the tests are ready, execute them via `run-env`. If the tests fail, read the error output, fix the code, and run them again.
