% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maximum_marginal_likelihood.R
\name{maximize_marginal_likelihood}
\alias{maximize_marginal_likelihood}
\title{Maximize Marginal Log-Likelihood with Convergence Check}
\usage{
maximize_marginal_likelihood(
  samples,
  grids,
  initial_lambda = 1,
  initial_tau = 1,
  initial_weights = rep(0, length(samples)),
  MC_iterations = 10000,
  max.iterations = 5,
  tol = 1e-04,
  parallel_computing = TRUE,
  seed = 1,
  boundaries = NULL
)
}
\arguments{
\item{samples}{Vector of sampled points.}

\item{grids}{Vector of centering grids.}

\item{initial_lambda}{Initial value for \eqn{\lambda} (default: 1).}

\item{initial_tau}{Initial value for \eqn{\tau} (default: 1).}

\item{initial_weights}{Initial weight vector (default: zeros of length \code{samples}).}

\item{MC_iterations}{Number of Monte Carlo iterations.}

\item{max.iterations}{Maximum number of iterations (default: 5).}

\item{tol}{Convergence tolerance (default: \code{1e-4}). If the relative change in
\eqn{\lambda} and \eqn{\tau} falls below this threshold, optimization stops early.}

\item{parallel_computing}{Boolean flag indicating whether parallelization should be
used for optimization (default: \code{TRUE}).}

\item{seed}{An integer that controls the randomness.}

\item{boundaries}{A vector including the boundaries}
}
\value{
A list containing:
\describe{
\item{lambda}{Optimized value of \eqn{\lambda}.}
\item{tau}{Optimized value of \eqn{\tau}.}
\item{max_marginal_likelihood}{Maximum marginal log-likelihood value.}
\item{converged}{Logical value indicating whether the optimization converged
before reaching \code{max.iterations}.}
}
}
\description{
This function optimizes the hyper-parameters \eqn{\lambda} and \eqn{\tau} by maximizing
the marginal log-likelihood using L-BFGS-B optimization. Instead of a grid search,
it efficiently finds the best parameters while checking for convergence.
}
\details{
\itemize{
\item Uses the \strong{L-BFGS-B} optimization method to efficiently find the best
\eqn{\lambda} and \eqn{\tau} values.
\item The optimization stops early if the relative change in \eqn{\lambda} and
\eqn{\tau} is smaller than \code{tol}.
\item If convergence is not achieved within \code{max.iterations}, the function
prints a warning and returns the best estimate.
}
}
\examples{
\dontrun{
result <- optimize_marginal_log_likelihood(
  centered_kernel_mat_samples, min_x, max_x, samples,
  MC_iterations = 10000, max.iterations = 50, tol = 1e-4
)
print(result)
}

}
