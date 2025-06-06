% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centered_kernel_matrix.R
\name{centered_kernel_matrix}
\alias{centered_kernel_matrix}
\title{Compute the Centered Kernel Matrix}
\usage{
centered_kernel_matrix(
  dimension,
  eval_points_1,
  eval_points_2,
  centering_grid,
  hurst_coef
)
}
\arguments{
\item{dimension}{Integer. Dimensionality of the input space.}

\item{eval_points_1}{Matrix or vector. First set of evaluation points (e.g., \eqn{\bx}).}

\item{eval_points_2}{Matrix or vector. Second set of evaluation points (e.g., \eqn{\bx'}).}

\item{hurst_coef}{Numeric scalar. The Hurst coefficient.}

\item{centering_points}{Matrix or vector. Grid used for kernel centering (e.g., \eqn{\bz_1, \ldots, \bz_n}).}
}
\value{
A numeric matrix of kernel evaluations.
}
\description{
This function computes the centered fractional Brownian motion (fBM) kernel matrix.
}
