#' Compute the RKHS Squared Error (MSE) in an RKHS Given the True Weights
#'
#' This function calculates the MSE using a Reproducing Kernel Hilbert Space (RKHS) formulation.
#'
#' @param w_hat_vec Numeric vector. Estimated weights.
#' @param w_vec Numeric vector. True weights.
#' @param kernel_matrix_at_samples Numeric matrix. The kernel matrix evaluated at the sample points.
#'
#' @return A numeric scalar representing the computed MSE.
#' @export
rkhs_se <- function(w_hat_vec, w_vec, kernel_matrix_at_samples){
  se <- t(w_hat_vec - w_vec) %*%
    kernel_matrix_at_samples %*%
    (w_hat_vec - w_vec)
  return(as.numeric(se))
}


#' Compute the Integrated Squared Error (ISE) for L2 Norm Given the True Densities
#'
#' This function calculates the Integrated Squared Error (ISE) between the true and estimated density functions using numerical integration.
#'
#' @param true_density_grid Numeric vector. The true density values evaluated at `grid` points.
#' @param estimated_density_grid Numeric vector. The estimated density values evaluated at `grid` points.
#'
#' @return A numeric scalar representing the computed ISE.
#' @export
l2_ise <- function(true_density_grid, estimated_density_grid){
  n <- length(true_density_grid)
  ise <- n^(-1)*sum((true_density_grid - estimated_density_grid)^2)
  return(as.numeric(ise))
}


#' Compute the Maximum Mean Discrepancy (MMD)
#'
#' This function calculates the Maximum Mean Discrepancy (MMD)
#'
#' @param kernel_matrix_at_grids Numeric matrix The kernel matrix evaluated at the grid points.
#' @param true_density_grids Numeric vector. The true density values evaluated at `grid` points.
#' @param estimated_density_grids Numeric vector. The estimated density values evaluated at `grid` points.
#'
#' @return A numeric scalar representing the estimated MMD on grids.
#' @export
mmd_grids <- function(kernel_matrix_at_grids, true_density_grids, estimated_density_grids){
  n <- dim(kernel_matrix_at_grids)[1]
  mmd_grids <- (n^(-2))* (t(true_density_grids - estimated_density_grids) %*% kernel_matrix_at_grids %*% (true_density_grids - estimated_density_grids))
  return(as.numeric(mmd_grids))
}

#' Compute the Maximum Mean Discrepancy (MMD) using samples
#'
#' This function calculates the Maximum Mean Discrepancy (MMD) using samples.
#'
#' @param kernel_matrix_at_samples Numeric matrix The kernel matrix evaluated at the sample points.
#' @param true_density_samples Numeric vector. The true density values evaluated at `sample` points.
#' @param estimated_density_samples Numeric vector. The estimated density values evaluated at `sample` points.
#' @param base_measure_weights Numeric vector. The base measure weights.
#'
#' @return A numeric scalar representing the estimated MMD on samples.
#' @export
mmd_samples <- function(kernel_matrix_at_samples, true_density_samples, estimated_density_samples, base_measure_weights){
  n <- dim(kernel_matrix_at_grids)[1]
  mmd_samples <- t(base_measure_weights * (true_density_samples - estimated_density_samples)) %*%
    kernel_matrix_at_samples %*%
    (base_measure_weights * (true_density_samples - estimated_density_samples))
  return(as.numeric(mmd_samples))
}


