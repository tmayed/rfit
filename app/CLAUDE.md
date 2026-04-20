# rfit Package

## Overview

rfit is an R package for probability distribution fitting using maximum likelihood estimation (MLE).

## Architecture

- Each distribution is stored in a separate file under `pkg/dist/`
- All functions for a distribution reside in its dedicated file
- **Implementation Standard:** Distributions MUST be implemented using standard R functions only. External libraries (e.g., `fitdistrplus`, `actuar`) are NOT permitted for distribution-specific logic (PDF, CDF, fitting, etc.).
- The package runs in a Docker container for reproducible execution

## Current Distributions

- Lognormal (`pkg/dist/lognormal.R`)
- Normal (`pkg/dist/normal.R`)
- Weibull (`pkg/dist/weibull.R`)
- Pareto (`pkg/dist/pareto.R`)
- GB2 (`pkg/dist/gb2.R`)
- dPLN (`pkg/dist/dpln.R`)
- Bradford (`pkg/dist/bradford.R`)
- Fisk (`pkg/dist/fisk.R`)
- Johnson SB (`pkg/dist/johnsonsb.R`)
- Johnson SL (`pkg/dist/johnsonsl.R`)
- Johnson SU (`pkg/dist/johnsonsu.R`)
- Kappa 4 (`pkg/dist/kappa4.R`)
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
