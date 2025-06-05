#' Compute Base Measure Weights
#'
#' Computes the base measure weights for either 1D (using grid midpoints) or 2D (using Voronoi tessellation).
#'
#' @param sample A numeric vector or matrix of sample points. For 1D, a vector; for 2D, a matrix/data.frame with two columns.
#' @param boundries A numeric vector (length 2) for 1D, or a 2x2 matrix for 2D defining min/max bounds of each axis.
#' @param dim Integer. Dimensionality of the sample: either 1 or 2.
#'
#' @return A numeric vector of base measure weights.
#' @export
#'
#' @examples
#' get_base_measures(sample = runif(100), boundries = c(0,1), dim = 1)
#' get_base_measures(sample = matrix(runif(200), ncol=2), boundries = matrix(c(0,1,0,1), nrow=2), dim = 2)
get_base_measures <- function(sample, boundries, dim = 1) {

  if (dim == 1) {
    # Validate input
    if (length(boundries) != 2) stop("For 1D, 'boundries' must be a length-2 vector.")
    sampled_x <- sort(sample)

    # Compute middle points between sorted x
    min_x <- boundries[1]
    max_x <- boundries[2]
    middle_points <- get_middle_points_grid(min_x, sampled_x, max_x)

    # Base measure is length of sub-intervals
    base_measure_weights <- diff(middle_points)

  } else if (dim == 2) {
    # Validate input
    if (!is.matrix(boundries) || dim(boundries)[1] != 2 || dim(boundries)[2] != 2)
      stop("For 2D, 'boundries' must be a 2x2 matrix.")
    if (!is.matrix(sample) || ncol(sample) != 2)
      stop("For 2D, 'sample' must be a matrix with 2 columns.")

    x1 <- sample[,1]
    x2 <- sample[,2]
    min_x1 <- boundries[1,1]
    max_x1 <- boundries[1,2]
    min_x2 <- boundries[2,1]
    max_x2 <- boundries[2,2]

    # Use Voronoi diagram from deldir
    base_measure_weights <- deldir::deldir(x1, x2,
                                           rw = c(min_x1, max_x1, min_x2, max_x2))$summary$dir.area
  } else {
    stop("For the moment function only supports 1D or 2D input.")
  }

  return(base_measure_weights)
}

