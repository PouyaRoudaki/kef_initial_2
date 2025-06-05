#' Estimate Weights Using the Barzilai-Borwein Method (Optimized with Rcpp)
#'
#' This function estimates the weight vector using the Barzilai-Borwein method,
#' a numerical approach to solve nonlinear systems of equations. It optimizes
#' performance by leveraging Rcpp for computation and adaptive sub-sampling for
#' better initial values in the Barzilai-Borwein method.
#'
#' @param lambda_hat A scalar representing the estimated lambda parameter.
#' @param tau_hat A scalar representing the estimated tau parameter.
#' @param centered_kernel_mat_at_sampled A matrix representing the centered kernel at sampled points.
#' @param sampled_x A numeric vector of sampled points.
#' @param min_x The minimum domain value.
#' @param max_x The maximum domain value.
#' @param prior_variance_p_vector A numeric vector for prior variance probabilities. Default is NULL.
#' @param print_trace Logical; if TRUE, prints progress updates.
#' @param init Logical; if TRUE, the weights calculation has a sub-sampling method for initialization.
#'
#' @return A numeric vector of estimated weights.
#' @export
get_weights <- function(lambda_hat, tau_hat, centered_kernel_mat_at_sampled,
                                        sampled_x, min_x, max_x,
                                        prior_variance_p_vector = NULL,
                                        print_trace = FALSE,
                                        init = FALSE,
                                        q_with_base = TRUE,
                                        with_prob_in_var = TRUE,
                                        normalised_q = TRUE,
                                        normalised_p = TRUE,
                                        p_with_base = FALSE) {

  n <- nrow(centered_kernel_mat_at_sampled)  # Number of sampled points

  # Wrapper for BBsolve using the Rcpp function
  s_function <- function(weight_hat_vec) {
    result <- get_s_function(weight_hat_vec, lambda_hat, tau_hat,
                             centered_kernel_mat_at_sampled, sampled_x, min_x, max_x,
                             prior_variance_p_vector,
                             q_with_base,
                             with_prob_in_var,
                             normalised_q,
                             normalised_p,
                             p_with_base)  # Pass prior_variance_p_vector

    return(as.numeric(result))  # Ensure it's a standard numeric vector
  }

  # Default initial weights for BBsolve.
  initial_weights <- rep(0, n)

  # Solve using BBsolve for full dataset
  result <- BB::BBsolve(par = as.numeric(initial_weights),
                    fn = s_function,
                    control = list(maxit = 10000,
                                   tol = 1e-4,
                                   trace = print_trace))
  if (print_trace) {
    print(result$message)
  }

  return(result$par)  # Return optimized weight vector
}
