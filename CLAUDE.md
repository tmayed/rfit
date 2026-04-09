# rfit Package

## Overview

rfit is an R package for probability distribution fitting using maximum likelihood estimation (MLE).

## Architecture

- Each distribution is stored in a separate file under `app/dist/`
- All functions for a distribution reside in its dedicated file
- The package runs in a Docker container for reproducible execution

## Current Distributions

- Lognormal (`app/dist/lognormal.R`)
- Pareto (`app/dist/pareto.R`)

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

## File Structure

```
app/
├── dist/
│   ├── lognormal.R
│   └── pareto.R
├── rfit.R
└── Dockerfile
tests/
└── test.R
poc/
└── poc.R
```

## Usage

Shell into the Docker container and load the modules:

```r
source("app/dist/lognormal.R")
source("app/dist/pareto.R")

# Fit a distribution
fit <- lognormal_fit(data)

# Evaluate PDF
densities <- lognormal_pdf(x, fit)
```
