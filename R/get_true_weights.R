#' Compute Approximation of True Weights for Density Estimation
#'
#' This function estimates the true weights (`w`) used in a density estimation model.
#' It computes the centered kernel matrices, evaluates the density function, and
#' optimizes the weights using least squares error minimization.
#'
#' @param sample A numeric vector of sampled points.
#' @param grid A numeric vector representing the grid over which the density is estimated.
#' @param true_density_at_grid A numeric vector of true density values at each grid point.
#'
#' @return A numeric vector of estimated weights (`w`) corresponding to the sampled points.
#' @importFrom pracma trapz Norm
#' @importFrom stats optim
#'
#' @examples
#' sample <- seq(-3, 3, length.out = 10)
#' grid <- seq(-3, 3, length.out = 100)
#' true_density_at_grid <- dnorm(grid)  # Example: Gaussian density
#' w_estimates <- get_true_weights(sample, grid, true_density_at_grid)
#' print(w_estimates)
#'
#' @export
get_true_weights <- function(samples, grids, true_density_grids) {

  # Check samples and grids dimensionality
  if (is.vector(samples)) {
    dimension <- 1
    if (!is.vector(grids)) stop("If samples is a vector, grids must also be a vector.")
    if (length(samples) > length(grids)) {
      warning("grids has fewer points than the samples. Consider increasing grids size for better resolution.")
    }

  } else if (is.matrix(samples)) {
    dimension <- ncol(samples)
    if (!is.matrix(grids) || ncol(grids) != dimension)
      stop("If samples is a matrix, grids must also be a matrix with the same number of columns.")
    if (nrow(samples) > nrow(grids)) {
      warning("grids has fewer points than the samples. Consider increasing grids size for better resolution.")
    }
  } else {
    stop("samples must be either a numeric vector or a numeric matrix.")
  }

  # Compute the centered kernel matrix at sampled points
  centered_kernel_mat_samples <- centered_kernel_matrix(
    dimension = dimension,
    eval_points_1 = samples,
    eval_points_2 = samples,
    centering_grid = grids,
    hurst_coef = 0.5
  )

  # Compute the centered kernel matrix at grid points
  centered_kernel_mat_grids <- centered_kernel_matrix(
    dimension = dimension,
    eval_points_1 = samples,
    eval_points_2 = grids,
    centering_grid = grids,
    hurst_coef = 0.5
  )

  # Compute the diagonal of the centered kernel matrix at grid points
  centered_kernel_self_grid <- diag(centered_kernel_matrix(
    dimension = dimension,
    eval_points_1 = grids,
    eval_points_2 = grids,
    centering_grid = grids,
    hurst_coef = 0.5
  ))

  # Given true f(x) for all grid points
  f_x <- true_density_grids

  # Estimating the base measure
  base_measure_weights <- get_base_measures(samples,
                                            c(min(grids), max(grids)),
                                            dimension = dimension)

  # Function to compute predicted f(x) given weights w
  predicted_f <- function(w) {
    density <- get_dens(centered_kernel_mat_samples,
                        centered_kernel_mat_grids,
                        centered_kernel_self_grids,
                        samples,
                        grids,
                        base_measure_weights,
                        dimension,
                        lambda = 1,
                        w)$grids

    return(density)
  }

  # Define loss function for optimization
  loss_function <- function(w) {
    return((pracma::Norm(predicted_f(w) - f_x, 2))^2)  # Least squares error
  }

  # Choose method for solving the nonlinear system of equations
  initial_w <- rep(0, length(samples))  # Initial guess for weights


  result <- optim(initial_w, loss_function, method = "L-BFGS-B")
  approx_true_w <- result$par


  return(approx_true_w)
}
