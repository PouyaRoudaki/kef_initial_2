% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kef.R
\name{kef}
\alias{kef}
\title{Kernel Density Estimation Using the KEF Method}
\usage{
kef(samples, grids, lambda, tau, boundaries = NULL)
}
\arguments{
\item{samples}{A numeric vector representing the observed samples points.}

\item{grids}{A numeric vector representing the evaluation grids.}

\item{lambda}{A numeric scalar for the lambda parameter (regularization term).}

\item{tau}{A numeric scalar for the tau parameter.}

\item{boundaries}{A numeric vector (for 1D) or matrix (for multi-dimensional input) specifying the domain boundaries. Optional; if not provided, the domain boundaries are extended by 10\% around the samples range.}
}
\value{
A list containing:
\itemize{
\item \code{weights} - The estimated weight vector.
\item \code{probs_samples} - Density estimates at samples points.
\item \code{probs_grids} - Density estimates at grids points (only if \code{length(samples) < 1000}).
\item \code{time} - The function's execution time.
}
}
\description{
This function estimates density using kernel exponential family (KEF) methods.
It constructs a kernel matrix and computes weights using the Barzilai-Borwein
optimization method. If the samples size is large, it avoids computational complexity
by returning density estimates only at samples points.
}
\details{
If the samples size is \code{>= 1000}, the function \strong{only returns densities at samples points}
to avoid computational complexity. A warning is issued in this case.
}
\examples{
samples <- rnorm(500)
grids <- seq(-3, 3, length.out = 100)
lambda <- 0.1
tau <- 0.5
result <- kef(samples = samples, grids = grids,lambda = lambda, tau = tau)
plot(grids, result$grids, type = "l", main = "Estimated Density", ylab = "Density")
}
