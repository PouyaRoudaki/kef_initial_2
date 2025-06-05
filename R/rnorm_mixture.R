#' Generate Random Samples from a Normal Mixture Model
#'
#' This function generates random samples from a mixture of normal distributions.
#' The user specifies the number of samples, the means, standard deviations,
#' and mixture weights for each component of the mixture.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param means A numeric vector of length 4 specifying the means of the normal components.
#' @param sds A numeric vector of length 4 specifying the standard deviations of the normal components.
#' @param mixture_weights A numeric vector of length 4 specifying the mixture weights for the components.
#'        The weights should sum to 1.
#'
#' @return A numeric vector of length `n` containing the generated samples from the normal mixture model.
#' @export
#'
#' @examples
#' # Example usage:
#' n <- 1000
#' means <- c(0, 5, 10, 15)
#' sds <- c(1, 1, 2, 2)
#' mixture_weights <- c(0.25, 0.25, 0.25, 0.25)
#' samples <- normal_mixture(n, means, sds, mixture_weights)
#' hist(samples, breaks = 30, main = "Histogram of Normal Mixture Samples")
rnorm_mixture <- function(n, means, sds, mixture_weights) {
  # Ensure that the mixture weights sum to 1 (optional check)
  if (sum(mixture_weights) != 1) {
    stop("The mixture weights must sum to 1.")
  }

  if(length(means) != length(sds) || length(means) != length(probabilities)) {
    stop("Lengths of means, sds, and probabilities must be equal")
  }

  # Generate component assignments based on mixture probabilities
  components <- sample(1:length(probabilities), size = n, replace = TRUE, prob = probabilities)

  # Generate samples
  samples <- mapply(function(component, n) {
    rnorm(n, mean = means[component], sd = sds[component])
  }, component = components, n = rep(1, n))

  return(samples)
}

