#' Compute Marginal Log Likelihood Using Rcpp
#'
#' This function computes the marginal log likelihood using a Monte Carlo approach
#' with optional parallel computing.
#'
#' @param centered_kernel_mat_at_samples A matrix of centered kernel values at sampled points.
#' @param samples A numeric vector of sampled points.
#' @param min_x The minimum x value.
#' @param max_x The maximum x value.
#' @param p_vec A probability vector (default: uniform distribution).
#' @param lambda A scalar for the lambda hyperparameter.
#' @param tau A scalar for the tau hyperparameter.
#' @param std_rnorm_matrix A matrix of standard normal random values for Monte Carlo sampling.
#' @param MC_iterations The number of Monte Carlo iterations.
#' @param parallel_computing Logical; if TRUE, enables parallel computing.
#'
#' @return The computed log of the marginal likelihood.
#' @export
marginal_log_likelihood <- function(centered_kernel_mat_samples,
                                    samples,
                                    min_x,
                                    max_x,
                                    p_vec = rep(1, nrow(centered_kernel_mat_at_samples)),
                                    lambda,
                                    tau,
                                    std_rnorm_matrix,
                                    MC_iterations,
                                    parallel_computing = TRUE) {

  # Call the C++ function using `.Call()`
  .Call("_kefV1_marginal_log_likelihood",
        centered_kernel_mat_samples,
        samples,
        min_x,
        max_x,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
}




#' Optimize Marginal Log-Likelihood with Convergence Check
#'
#' This function optimizes the hyper-parameters \eqn{\lambda} and \eqn{\tau} by maximizing
#' the marginal log-likelihood using L-BFGS-B optimization. Instead of a grid search,
#' it efficiently finds the best parameters while checking for convergence.
#'
#' @param centered_kernel_mat_samples The kernel matrix centered at sampled points.
#' @param min_x Minimum value of the sampled domain.
#' @param max_x Maximum value of the sampled domain.
#' @param samples Vector of sampled points.
#' @param initial_lambda Initial value for \eqn{\lambda} (default: 1).
#' @param initial_tau Initial value for \eqn{\tau} (default: 1).
#' @param initial_w Initial weight vector (default: zeros of length `sampled_x`).
#' @param MC_iterations Number of Monte Carlo iterations.
#' @param max.iterations Maximum number of iterations (default: 5).
#' @param tol Convergence tolerance (default: `1e-4`). If the relative change in
#' \eqn{\lambda} and \eqn{\tau} falls below this threshold, optimization stops early.
#' @param parallel_computing Boolean flag indicating whether parallelization should be
#' used for optimization (default: `TRUE`).
#' @param seed An integer that controls the randomness.
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda}{Optimized value of \eqn{\lambda}.}
#'   \item{tau}{Optimized value of \eqn{\tau}.}
#'   \item{max_marginal_likelihood}{Maximum marginal log-likelihood value.}
#'   \item{converged}{Logical value indicating whether the optimization converged
#'   before reaching `max.iterations`.}
#' }
#'
#' @details
#' - Uses the **L-BFGS-B** optimization method to efficiently find the best
#'   \eqn{\lambda} and \eqn{\tau} values.
#' - The optimization stops early if the relative change in \eqn{\lambda} and
#'   \eqn{\tau} is smaller than `tol`.
#' - If convergence is not achieved within `max.iterations`, the function
#'   prints a warning and returns the best estimate.
#'
#' @examples
#' \dontrun{
#' result <- optimize_marginal_log_likelihood(
#'   centered_kernel_mat_samples, min_x, max_x, samples,
#'   MC_iterations = 10000, max.iterations = 50, tol = 1e-4
#' )
#' print(result)
#' }
#'
#' @export
maximise_marginal_likelihood <- function(samples,
                                         boundaries,
                                         initial_lambda = 1,
                                         initial_tau = 1,
                                         initial_weights = rep(0, length(samples)),
                                         MC_iterations = 10000,
                                         max.iterations = 5,
                                         tol = 1e-4,  # Convergence tolerance
                                         parallel_computing = TRUE,
                                         seed = 1) {

  # Setting default boundaries with 10% padding if the boundaries are not provided.
  if (is.null(boundaries)) {
    if (is.vector(samples)) {
      n <- length(samples)
      min_x <- min(samples)
      max_x <- max(samples)
      padding <- 0.1 * (max_x - min_x)
      boundaries <- c(min_x - padding, max_x + padding)
    } else if (is.matrix(samples)) {
      n <- dim(samples)[1]
      d <- ncol(samples)
      boundaries <- matrix(NA, nrow = d, ncol = 2)
      for (j in seq_len(d)) {
        min_j <- min(samples[, j])
        max_j <- max(samples[, j])
        padding_j <- 0.1 * (max_j - min_j)
        boundaries[j, ] <- c(min_j - padding_j, max_j + padding_j)
      }
    } else {
      stop("Unsupported samples type: expected vector or matrix.")
    }
  }

  # Check samples and grids dimensionality
  if (is.vector(samples)) {
    dimension <- 1

    sort_index <- order(samples)
    unsort_index <- order(sort_index)  # this gives the positions to map back

    # Sorted samples for computing base measures
    # This is required because further we used trapz for 1 dimensional normalizing
    # constant and trapz function needs ordered input. Instead of doing this
    # inside our functions it's better to do it here.
    # Get the sort index and reverse mapping
    samples <- samples[sort_index]

    if (!is.vector(grids)) stop("If samples is a vector, grids must also be a vector.")
    if (!is.vector(boundaries) || length(boundaries) != 2)
      stop("For 1D input, boundaries must be a vector of length 2 specifying the min and max of the domain.")
    if (length(samples) > length(grids)) {
      warning("grids has fewer points than the samples. Consider increasing grids size for better resolution.")
    }

  } else if (is.matrix(samples)) {
    dimension <- ncol(samples)
    if (!is.matrix(grids) || ncol(grids) != dimension)
      stop("If samples is a matrix, grids must also be a matrix with the same number of columns.")
    if (!is.matrix(boundaries) || nrow(boundaries) != dimension || ncol(boundaries) != 2)
      stop("For d-dimensional input, boundaries must be a d × 2 matrix with each row specifying the min and max of the domain for one dimension.")
    if (nrow(samples) > nrow(grids)) {
      warning("grids has fewer points than the samples. Consider increasing grids size for better resolution.")
    }

  } else {
    stop("samples must be either a numeric vector or a numeric matrix.")
  }


  # Generate matrix with each row independently sampled from Normal(0,1)
  set.seed(seed)
  std_rnorm_matrix <- matrix(rnorm(MC_iterations * n, mean = 0, sd = 1),
                             nrow = MC_iterations,
                             ncol = n)

  # Initializations
  t <- 1
  lambda <- initial_lambda
  tau <- initial_tau
  weight_vec <- initial_weights
  converged <- FALSE


  # Estimating the base measure
  base_measure_weights <- get_base_measures(samples, boundaries, dimension = dimension)

  # Compute the centered kernel matrix at samples
  centered_kernel_mat_samples <- centered_kernel_matrix(
    dimension = dimension,
    eval_points_1 = samples,
    eval_points_2 = samples,
    centering_grid = grids,
    hurst_coef = 0.5
  )

  # Initial density
  dens_vec <- get_dens_wo_grid(centered_kernel_mat_samples,
                               samples,
                               base_measure_weights,
                               dimension,
                               lambda,
                               weight_vec)

  p_vec <- (dens_vec * base_measure_weights) / sum(dens_vec * base_measure_weights)

  p_vec_init <- as.numeric(p_vec)

  repeat {

    cat(paste0("Iteration: ", t, "\n"))

    # Store old values for convergence check
    lambda_old <- lambda
    tau_old <- tau

    # Define objective function for optimization
    objective_function <- function(params) {
      log_lambda <- params[1]
      log_tau <- params[2]

      lambda <- exp(log_lambda)
      tau <- exp(log_tau)

      - marginal_log_likelihood(
        centered_kernel_mat_samples,
        samples,
        p_vec,
        lambda,
        tau,
        std_rnorm_matrix,
        MC_iterations,
        parallel_computing)
    }

    cat(paste0("Initial lambda: ", 0.1,", Initial tau: ", 1e-4, "\n"))

    # Optimization using L-BFGS-B (bounded optimization)
    opt_result <- optim(
      par = c(log(0.1), log(1e-4)),  # Start close to expected values
      fn = objective_function,
      method = "L-BFGS-B",
      lower = c(log(1e-2), log(1e-6)),  # Lower bounds for lambda and tau
      upper = c(log(1e2), log(1e2))     # Upper bounds
    )

    # Retrieve optimal lambda and tau
    lambda <- exp(opt_result$par[1])
    tau <- exp(opt_result$par[2])


    cat(paste0("Optimized lambda: ", lambda, ", tau: ", tau, ", MLL: ", -opt_result$value,
               ", The ratio: ", lambda^2/tau ,"\n"))

    #cat(paste0("lambda_old: ",lambda_old,"tau_old: ",tau_old,".\n"))
    # Convergence check: Stop if parameters don't change significantly
    delta_lambda <- abs(lambda - lambda_old)
    delta_tau <- abs(tau - tau_old)

    #cat(paste0("delta_lambda: ",delta_lambda,"delta_tau: ",delta_tau,".\n"))

    if ( max(delta_lambda,delta_tau) < tol) {
      converged <- TRUE
      cat("✔ Convergence reached: Estimated Lambda and Tau by Maximum Marginal Likelihood are stable.\n")
      break
    }

    # Update weights
    weight_vec <- get_weights(lambda = lambda,
                          tau = tau,
                          centered_kernel_mat_samples = centered_kernel_mat_samples,
                          samples = samples,
                          base_measure_weights = base_measure_weights,
                          dimension = dimension,
                          print_trace = FALSE)

    # Update density and p_vec
    dens_vec <- get_dens_wo_grid(centered_kernel_mat_samples,
                                 samples,
                                 base_measure_weights,
                                 dimension,
                                 lambda,
                                 weight_vec)

    p_vec <- as.numeric((dens_vec * base_measure_weights) / sum(dens_vec * base_measure_weights))

    t <- t + 1

    if (t > max.iterations) {
      cat("Error: Max iterations reached without full convergence.\n")
      break
    }
  }

  return(list(lambda = lambda,
              tau = tau,
              max_marginal_log_likelihood = -opt_result$value,
              converged = converged,
              std_rnorm_matrix = std_rnorm_matrix,
              p_vec = p_vec_init))
}
