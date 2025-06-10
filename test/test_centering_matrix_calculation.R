centered_kernel_matrix_R <- function(first_mat_kernel, second_mat_kernel, centering_grid, hurst_coef) {
  n0 <- nrow(first_mat_kernel)
  n1 <- nrow(second_mat_kernel)
  n2 <- nrow(centering_grid)

  # Initialize all required matrices
  term1_matrix <- matrix(0, nrow = n0, ncol = n1)
  term2_matrix <- matrix(0, nrow = n0, ncol = n2)
  term3_matrix <- matrix(0, nrow = n1, ncol = n2)
  term4_matrix <- matrix(0, nrow = n2, ncol = n2)

  # Compute term1_matrix: ((||x - y||^2)^H)
  for (i in 1:n0) {
    for (j in 1:n1) {
      diff <- first_mat_kernel[i, ] - second_mat_kernel[j, ]
      norm2_sq <- sqrt(sum(diff^2))
      term1_matrix[i, j] <- norm2_sq^(2*hurst_coef)
    }
  }

  # Compute term2_matrix: ((||x - z||^2)^H)
  for (i in 1:n0) {
    for (j in 1:n2) {
      diff <- first_mat_kernel[i, ] - centering_grid[j, ]
      norm2_sq <- sqrt(sum(diff^2))
      term2_matrix[i, j] <- norm2_sq^(2*hurst_coef)
    }
  }

  # Compute term3_matrix: ((||y - z||^2)^H)
  for (i in 1:n1) {
    for (j in 1:n2) {
      diff <- second_mat_kernel[i, ] - centering_grid[j, ]
      norm2_sq <- sqrt(sum(diff^2))
      term3_matrix[i, j] <- norm2_sq^(2*hurst_coef)
    }
  }

  # Compute term4_matrix: ((||z_i - z_j||^2)^H)
  for (i in 1:n2) {
    for (j in 1:n2) {
      diff <- centering_grid[i, ] - centering_grid[j, ]
      norm2_sq <- sqrt(sum(diff^2))
      term4_matrix[i, j] <- norm2_sq^(2*hurst_coef)
    }
  }

  # Compute row means of term2_matrix
  row_means_term2 <- numeric(n0)
  for (i in 1:n0) {
    row_means_term2[i] <- mean(term2_matrix[i, ])
  }

  # Compute row means of term3_matrix
  row_means_term3 <- numeric(n1)
  for (j in 1:n1) {
    row_means_term3[j] <- mean(term3_matrix[j, ])
  }

  # Compute mean of term4_matrix
  mean_term4 <- mean(term4_matrix)

  # Compute final result matrix
  result_matrix <- matrix(0, nrow = n0, ncol = n1)
  for (i in 1:n0) {
    for (j in 1:n1) {
      result_matrix[i, j] <- -0.5 * (
        term1_matrix[i, j] - row_means_term2[i] - row_means_term3[j] + mean_term4
      )
    }
  }

  return(result_matrix)
}

x_grid <- seq(-2.5, 2.5, length.out = 40)
y_grid <- seq(-2.5, 2.5, length.out = 40)
grid_mat <- as.matrix(expand.grid(x_grid, y_grid))

set.seed(7)

# === Generate samples from a standard 2D normal ===
n <- 100
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)  # Correlated normal
samples <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

R_result <- centered_kernel_matrix_R(first_mat_kernel = samples,
                         second_mat_kernel = grid_mat,
                         centering_grid = grid_mat,
                         hurst_coef = 0.5)

CPP_result <- centered_kernel_matrix(dimension = 2,
                                     eval_points_1 = samples,
                                     eval_points_2 = grid_mat,
                                     centering_grid = grid_mat,
                                     hurst_coef = 0.5)


head(R_result[,1:6])
head(CPP_result[,1:6])


grid_mat[1561,]

dim(CPP_result)

CPP_result[,1561] %*% kef_res$weights
