\name{bspline.curve.admm_lambdas}
\alias{bspline.curve.admm_lambdas}
\title{B-spline Curve Fitting using ADMM}
\usage{
bspline.curve.admm_lambdas(y, D, B, lambdas = NULL, lam_max = 100,
lam_min = 1e-05, n_lambda = 100, max_iter = 500, epsilon = 1e-5, eta_c = 0.1)
}
\arguments{
\item{y}{A numeric vector or matrix representing the observed data.}
\item{D}{A difference matrix used for regularization in the ADMM algorithm. It is typically generated using \code{bspline_jump(knots, order)}.}
\item{B}{A B-spline basis matrix corresponding to the design matrix. It is typically generated using \code{bsplines(t, knots, order)}.}
\item{lambdas}{A vector of regularization parameters. If \code{NULL}, a sequence of lambda values is automatically generated.}
\item{lam_max}{The maximum lambda value for automatic generation (only used if \code{lambdas} is \code{NULL}).}
\item{lam_min}{The minimum lambda value for automatic generation (only used if \code{lambdas} is \code{NULL}).}
\item{n_lambda}{The number of lambda values to generate (only used if \code{lambdas} is \code{NULL}).}
\item{max_iter}{The maximum number of iterations for the ADMM algorithm.}
\item{epsilon}{The convergence threshold for the iterative algorithm. The algorithm stops when the absolute difference between the penalized residual sum of squares (RSS) in consecutive iterations is less than \code{epsilon}.}
\item{eta_c}{A constant factor to compute the ADMM parameter \code{eta}.}
}

\value{
A list containing results for all lambda values and model selection criteria. The list includes:

\describe{
\item{xi}{Estimated coefficient matrix for each lambda value.}
\item{Bb}{Fitted values computed as \code{B \%*\% xi} for each lambda.}
\item{alpha}{Auxiliary variable estimates from the ADMM algorithm.}
\item{eta}{The computed ADMM parameter \code{eta}.}
\item{lambda}{The sequence of regularization parameters used in fitting.}
\item{u}{Dual variables from the ADMM algorithm.}
\item{iteration}{Number of iterations performed for convergence.}
\item{obj_f}{Final penalized residual sum of squares (RSS) for each lambda.}
\item{aic}{A vector of Akaike Information Criterion (AIC) values for all lambda values.}
\item{bic}{A vector of Bayesian Information Criterion (BIC) values for all lambda values.}
\item{r}{The total number of lambda values considered.}
\item{best_aic}{The model results corresponding to the lambda with the lowest AIC.}
\item{best_bic}{The model results corresponding to the lambda with the lowest BIC.}
}
}

\description{
This function performs B-spline curve fitting using the Alternating Direction Method of Multipliers (ADMM). It solves for the B-spline coefficients under a regularization scheme that encourages smoothness by applying a group lasso penalty. The function returns fitted coefficients and model selection criteria across a range of regularization parameters.
}
\examples{
set.seed(123)
n <- 100
order <- 3
dimension <- 10

# Generate time points

t <- seq(0, 1, length.out = n)

# Compute B-spline basis and difference matrix

knots <- knots_quantile(t, dimension, order)
xi <- matrix(rnorm(2 * dimension), ncol = 2)
B <- bsplines(t, knots, order)
D <- bspline_jump(knots, order)

# Generate synthetic response data

f <- matrix(0, n, 2)
for(i in 1:n)
  f[i, ] = crossprod(B[i, ], xi)
e <- matrix(rnorm(2 * n, sd = 0.1), ncol = 2)
y <- f + e

# Perform B-spline curve fitting using ADMM

fit <- bspline.curve.admm_lambdas(y, D, B)

# Plot observed data and fitted curve

best_index = which.min(fit$bic)
plot(y, col = "gray")
# true control points and curve
text(xi, col = "black", cex = 1)
lines(f, col = "black")
# fitted control points and curve
text(fit[[best_index]]$xi, col = "red", cex = 1)
lines(fit[[best_index]]$Bb, col = "red", lty = 2)
title("B-spline Curve Fitting via ADMM")
}

