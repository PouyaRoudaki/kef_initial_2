#' Estimate Weights Using the Barzilai-Borwein Method (Optimized with Rcpp)
#'
#' This function estimates the weight vector using the Barzilai-Borwein method,
#' a numerical approach to solve nonlinear systems of equations. It optimizes
#' performance by leveraging Rcpp for computation and adaptive sub-sampling for
#' better initial values in the Barzilai-Borwein method.
#'
#' @param lambdas A scalar representing the lambda parameter.
#' @param tau A scalar representing the tau parameter.
#' @param centered_kernel_mat_samples A matrix representing the centered kernel at sampled points.
#' @param samples A numeric vector of sampled points.
#' @param base_measure_weights A numeric vector of base measures for sample points.
#' @param prior_variance_p_vector A numeric vector for prior variance probabilities. Default is NULL.
#' @param print_trace Logical; if TRUE, prints progress updates.
#'
#' @return A numeric vector of estimated weights.
#' @importFrom BB BBsolve
#' @export
get_weights <- function(lambda,
                        tau,
                        centered_kernel_mat_samples,
                        samples,
                        base_measure_weights,
                        dimension,
                        print_trace = FALSE) {


  # Number of sampled points
  n <- nrow(centered_kernel_mat_samples)

  # Wrapper for BBsolve using the Rcpp function
  s_function <- function(weight_vec) {
    result <- get_s_function(weight_vec,
                             lambda,
                             tau,
                             centered_kernel_mat_samples,
                             samples,
                             base_measure_weights,
                             dimension)

    # Ensure it's a standard numeric vector
    return(as.numeric(result))
  }

  # Default initial weights for BBsolve.
  initial_weights <- rep(0, n)

  # Solve using BBsolve for full dataset
  result <- BB::BBsolve(par = as.numeric(initial_weights),
                    fn = s_function,
                    control = list(maxit = 10000,
                                   tol = 1e-4,
                                   trace = print_trace))

  # Print trace
  if (print_trace) {
    print(result$message)
  }

  # Return optimized weight vector
  return(result$par)
}
