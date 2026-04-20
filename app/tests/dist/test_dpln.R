# Test suite for dPLN Distribution

source("../../pkg/dist/dpln.R")

# Test parameters
fit <- list(alpha = 2.0, beta = 2.0, nu = 0, tau = 1.0)

# 1. Test PDF
x <- c(0.1, 1, 10)
pdf_vals <- dpln_pdf(x, fit)
print("PDF values:")
print(pdf_vals)

# 2. Test CDF
cdf_vals <- dpln_cdf(x, fit)
print("CDF values:")
print(cdf_vals)

# 3. Test Random Generation and Fitting
set.seed(123)
data <- dpln_rand(500, fit)
print("Mean of generated data:")
print(mean(data))
print("Theoretical mean:")
print(dpln_mean(fit))

# Try fitting
print("Attempting to fit dPLN to generated data...")
fitted <- dpln_fit(data)
print("Fitted parameters:")
print(fitted)

# 4. Test Quantile
p <- c(0.25, 0.5, 0.75)
q <- dpln_quantile(p, fit)
print("Quantiles:")
print(q)
print("CDF at quantiles (should be close to p):")
print(dpln_cdf(q, fit))

# ---------------------------------------------------------
# Extended tests for dPLN conventions
# ---------------------------------------------------------

# Test parameters for moments (ensure moments exist)
fit2 <- list(alpha = 5.0, beta = 5.0, nu = 0, tau = 0.5)

# 1. Test SF / ISF
x2 <- 1.0
sf <- dpln_sf(x2, fit2)
cdf2 <- dpln_cdf(x2, fit2)
print(paste("SF + CDF at 1.0 (should be 1):", sf + cdf2))

p2 <- 0.3
isf <- dpln_isf(p2, fit2)
q2 <- dpln_quantile(1 - p2, fit2)
print(paste("ISF(0.3) == Quantile(0.7):", isf == q2))

# 2. Test Log functions
lx <- log(dpln_pdf(x2, fit2))
lpdf <- dpln_logpdf(x2, fit2)
print(paste("log(PDF) vs logpdf:", lx, lpdf))

lcdf <- dpln_logcdf(x2, fit2)
print(paste("log(CDF) vs logcdf:", log(cdf2), lcdf))

# 3. Test Moments
print(paste("Mean:", dpln_mean(fit2)))
print(paste("Variance:", dpln_var(fit2)))
print(paste("Std:", dpln_std(fit2)))
print(paste("Skewness:", dpln_skew(fit2)))
print(paste("Kurtosis:", dpln_kurtosis(fit2)))

# 4. Test Median and Interval
print(paste("Median:", dpln_median(fit2)))
print("95% Interval:")
print(dpln_interval(0.95, fit2))

# 5. Test Entropy and Expect
print(paste("Entropy:", dpln_entropy(fit2)))
print(paste("Expect(x^2):", dpln_expect(function(x) x^2, fit2)))
print(paste("Moment(2):", dpln_moment(2, fit2)))
