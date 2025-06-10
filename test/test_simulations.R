# Load required packages
packages <- c("ks", "quantreg", "BB", "pracma",
              "tidyverse", "dplyr", "parallel", "doParallel", "foreach", "kefV1")


library("quantreg")
library("kefV1")
library("BB")
library("ks")
library("pracma")
library("parallelly")
library("doParallel")
library("foreach")

# Set up parallel backend
num_cores <- detectCores()-1  # Use all cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the weights for the mixture distribution
mixture_weights = c(1/2, 1/2)

# Define the parameters for the normal distributions
# First distribution: N(0, 1)
means = c(0, 1.5)
sds = c(0.5, 0.1)


# Define parameters
min_x <- -3
max_x <- 3

n_grid <- 1000
grids <- seq(min_x, max_x, length.out = n_grid)

# Mixture Normal Example
# Define a matrix of normal densities for each mean and standard deviation
density_matrix_grids <- sapply(seq_along(means), function(i) {
  dnorm(grids, mean = means[i], sd = sds[i])
})

# Calculate the true density by taking the weighted sum of the columns
true_density_grids <- as.numeric(density_matrix_grids %*% mixture_weights)

true_density_df_grids <- data.frame(grid = grids, true_pdf = true_density_grids)


centered_kernel_self_grids_mat <- centered_kernel_matrix(dimension = 1, grids, grids, grids, 0.5)

centered_kernel_self_grids <- diag(centered_kernel_self_grids_mat)

# Uniform Example
#params_uniform <- list(min = min_x, max = max_x)
#density_characterization_uniform <- list(type = "uniform", parameters = params_uniform)
#true_density_grid <- true_density_function(grid, density_characterization_uniform)

n_iter <- 10
result_list <- vector("list", n_iter)

# Parallel execution of iterations
result_list <- foreach(i = 1:n_iter, .packages = c("quantreg","ks", "pracma", "foreach", "kefV1")) %dopar% {
  cat("Iteration:", i, "\n")

  n_sample <- 100
  set.seed(i)
  samples <- sort(rnorm_mixture(n_sample, means, sds, mixture_weights))


  # Compute kernel matrices
  centered_kernel_mat_samples <- centered_kernel_matrix(dimension = 1, samples, samples, grids, 0.5)
  centered_kernel_mat_grids <- centered_kernel_matrix(dimension = 1, samples, grids, grids, 0.5)


  # Compute true weights
  w_true <- get_true_weights(samples, grids, true_density_grids)

  # Kernel Density Estimation (KDE)
  start_time <- Sys.time()
  estimated_density_grids <- kde(x = samples, h = hpi(samples), eval.points = grids)$estimate
  pi_time <- difftime(Sys.time(), start_time, units = "secs")
  pi_ise <- l2_ise(true_density_grids, estimated_density_grids)
  pi_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  # Biased Cross Validation KDE
  start_time <- Sys.time()
  estimated_density_grids <- kde(x = samples, h = hscv(samples), eval.points = grids)$estimate
  bcv_time <- difftime(Sys.time(), start_time, units = "secs")
  bcv_ise <- l2_ise(true_density_grids, estimated_density_grids)
  bcv_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  # Least Square Cross Validation KDE
  start_time <- Sys.time()
  estimated_density_grids <- kde(x = samples, h = hlscv(samples), eval.points = grids)$estimate
  lscv_time <- difftime(Sys.time(), start_time, units = "secs")
  lscv_ise <- l2_ise(true_density_grids, estimated_density_grids)
  lscv_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  # Normal Scale KDE
  start_time <- Sys.time()
  estimated_density_grids <- kde(x = samples, h = hns(samples), eval.points = grids)$estimate
  ns_time <- difftime(Sys.time(), start_time, units = "secs")
  ns_ise <- l2_ise( true_density_grids, estimated_density_grids)
  ns_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  # Adaptive KDE (Adhoc)
  start_time <- Sys.time()
  estimated_density_grids <- akj(x = samples, z = grids)$dens
  adaptive_default_time <- difftime(Sys.time(), start_time, units = "secs")
  adaptive_default_ise <- l2_ise(true_density_grids, estimated_density_grids)
  adaptive_default_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)
  # Adaptive bandwidth: CV optim

  adaptive_cv_optim <- function(samples, grids, alpha_init, kappa_init, lower_bounds, upper_bounds) {

    # Define the objective function to minimize (cv_error)
    objective_function <- function(params) {
      alpha <- params[1]
      kappa <- params[2]

      # Adaptive kernel density estimation
      adaptive_result <- akj(samples, grids, alpha = alpha, kappa = kappa)$dens

      # Roughness (L2 norm)
      roughness <- pracma::trapz(grids, (adaptive_result)^2)

      # Optimized Cross Product Term using proper indexing
      leave_one_out_densities <- vapply(seq_along(samples), function(i) {
        # Remove only the i-th occurrence, not all instances of sample[i]
        sample_wo_i <- samples[-i]

        # Compute Adaptive kernel density estimation without the i-th sample point
        adaptive_result_wo_i <- akj(x = sample_wo_i, z = samples[i], alpha = alpha, kappa = kappa)$dens

        # Extract density values
        densities_wo_i <- adaptive_result_wo_i

        return(densities_wo_i)
      }, numeric(1))  # Ensures output is numeric

      # Cross-validation error calculation
      cross_product_sum <- sum(leave_one_out_densities, na.rm = TRUE)
      cv_error <- roughness - (2 / length(samples)) * cross_product_sum

      return(cv_error)
    }

    # Optimization using Nelder-Mead or BFGS
    result <- optim(
      par = c(alpha_init, kappa_init),   # Initial values
      fn = objective_function,           # Function to minimize
      method = "L-BFGS-B",               # Bounded optimization
      lower = lower_bounds,              # Lower bounds for alpha and kappa
      upper = upper_bounds               # Upper bounds for alpha and kappa
    )

    # Extract best parameters
    best_alpha <- result$par[1]
    best_kappa <- result$par[2]
    min_cv_error <- result$value

    # Return results
    return(list(best_alpha = best_alpha, best_kappa = best_kappa, min_cv_error = min_cv_error))
  }



  start_time_adap <- Sys.time()

  alpha_init <- 0.5  # Initial guess for alpha
  kappa_init <- 0.5  # Initial guess for kappa
  lower_bounds <- c(0.01, 0.01)  # Lower bounds for alpha and kappa
  upper_bounds <- c(5, 5)  # Upper bounds for alpha and kappa

  adaptive_cv_optim <- adaptive_cv_optim(samples, grids, alpha_init, kappa_init, lower_bounds, upper_bounds)

  estimated_density_grids <- akj(x = samples, z = grids, alpha = adaptive_cv_optim$best_alpha, kappa = adaptive_cv_optim$best_kappa)$dens
  end_time_adap <- Sys.time()
  adaptive_cv_time <- difftime(end_time_adap, start_time_adap, units = "secs")
  adaptive_cv_optim_ise <- l2_ise(true_density_grids,estimated_density_grids)
  adaptive_cv_optim_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  # KEF: Rule of Thumb
  kef_rot_lambda <- 1
  kef_rot_tau <- 1 / 1350
  kef_rot_ratio <- (kef_rot_lambda)^2 / kef_rot_tau
  kef_rot <- kef(samples, grids, lambda = kef_rot_lambda, tau = kef_rot_tau)
  kef_rot_time <- kef_rot$time
  estimated_density_grids <- kef_rot$dens_grids
  kef_rot_ise <- l2_ise( true_density_grids, estimated_density_grids)
  kef_rot_se <- rkhs_se(w_hat_vec = kef_rot$weights, w_vec = w_true, kernel_matrix_at_samples = centered_kernel_mat_samples)
  kef_rot_mmd_grids <- mmd_grids(centered_kernel_self_grids_mat, true_density_grids, estimated_density_grids)

  # KEF: Marginal Log Likelihood (MML) Optimization
  start_time_mml <- Sys.time()
  optimized_mml <- optimize_marginal_log_likelihood(
    centered_kernel_mat_samples, min_x, max_x, samples,
    initial_lambda = 1, initial_tau = 1, initial_w = rep(0, length(samples)),
    MC_iterations = 100000, max.iterations = 10, tol = 1e-1,
    parallel_computing = TRUE, seed = 4
  )
  kef_mml <- kef(samples, grids, lambda = optimized_mml$lambda, tau = optimized_mml$tau)
  kef_mml_time <- as.numeric(kef_mml$time) + as.numeric(difftime(Sys.time(), start_time_mml, units = "secs"))
  estimated_density_grids <- kef_mml$probs_grids
  kef_mml_ise <- l2_ise(grids, true_density_grids, estimated_density_grids)
  kef_mml_se <- rkhs_se(w_hat_vec = kef_mml$weights, w_vec = w_true, kernel_matrix_at_samples = centered_kernel_mat_samples)
  kef_mml_mmd_grids <- mmd_grids(centered_kernel_mat_grids, true_density_grids, estimated_density_grids)

  kef_mml_lambda <- optimized_mml$lambda
  kef_mml_tau <- optimized_mml$tau
  kef_mml_ratio <- (kef_mml_lambda)^2 / kef_mml_tau

  # Store results in a dataframe
  result_iter_i <- data.frame(
    method = c("fixed_pi", "fixed_bcv", "fixed_lscv", "fixed_ns",
               "adaptive_default","adaptive_cv_optim", "kef_rot"),
    time = c(pi_time, bcv_time, lscv_time, ns_time, adaptive_default_time, adaptive_cv_time, kef_rot_time),
    ISE = c(pi_ise, bcv_ise, lscv_ise, ns_ise, adaptive_default_ise, adaptive_cv_optim_ise , kef_rot_ise),
    MMD_error_on_grid = c(pi_mmd_grids, bcv_mmd_grids, lscv_mmd_grids, ns_mmd_grids, adaptive_default_mmd_grids, adaptive_cv_optim_mmd_grids , kef_rot_mmd_grids),
    RKHS_SE = c(NA, NA, NA, NA, NA, NA, kef_rot_se),
    lambda = c(NA, NA, NA, NA, NA, NA, kef_rot_lambda),
    tau = c(NA, NA, NA, NA, NA, NA, kef_rot_tau),
    ratio = c(NA, NA, NA, NA, NA, NA, kef_rot_ratio)
  )

  return(result_iter_i)
}

# Stop parallel cluster
stopCluster(cl)

#saveRDS(result_list, "bimodal.rds")

# Combine results into one data frame
df_combined <- bind_rows(result_list, .id = "iteration")

# Compute summary statistics
result_summary <- df_combined %>%
  group_by(method) %>%
  summarise(
    time_mean = mean(time, na.rm = TRUE),
    time_se = sd(time, na.rm = TRUE),
    MISE = mean(ISE, na.rm = TRUE),
    MedISE = median(ISE, na.rm = TRUE),
    ISE_se = sd(ISE, na.rm = TRUE),
    RKHS_MSE = mean(RKHS_SE, na.rm = TRUE),
    RKHS_MedSE = median(RKHS_SE, na.rm = TRUE),
    RKHS_SE_se = sd(RKHS_SE, na.rm = TRUE),
    MMD_Mean = mean(MMD_error_on_grid, na.rm = TRUE),
    MMD_Med = median(MMD_error_on_grid, na.rm = TRUE),
    MMD_se = sd(MMD_error_on_grid, na.rm = TRUE),
    lambda_log_mean = mean(log10(lambda), na.rm = TRUE),
    lambda_log_se = sd(log10(lambda), na.rm = TRUE),
    tau_log_mean = mean(log10(tau), na.rm = TRUE),
    tau_log_se = sd(log10(tau), na.rm = TRUE),
    ratio_log_mean = mean(log10(ratio), na.rm = TRUE),
    ratio_log_se = sd(log10(ratio), na.rm = TRUE)
  )

# Save results
#write.csv(result_summary, file = "/results/summary_bimodal.csv")
