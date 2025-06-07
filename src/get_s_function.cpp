#include <RcppArmadillo.h>
#include "density.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Compute the function s(weight_vec)
// [[Rcpp::export]]
arma::vec get_s_function(const arma::vec& weight_vec,
                         double lambda,
                         double tau,
                         const arma::mat& centered_kernel_mat_samples,
                         const arma::vec& samples,
                         const arma::vec& base_measure_weights,
                         double dimension) {

  // Sample size
  double n = samples.n_elem;

  // Compute densities
  arma::vec dens = get_dens_wo_grid(centered_kernel_mat_samples,
                                    samples,
                                    base_measure_weights,
                                    dimension,
                                    lambda,
                                    weight_vec);

  // Compute density-based probabilities
  arma::vec dens_sample_via_base = dens % base_measure_weights;

  // Compute normalized probability vector
  arma::vec prob_sample_via_base = dens_sample_via_base / sum(dens_sample_via_base);

  // Compute the function s(weight_vec)
  arma::vec s = lambda * (sum(centered_kernel_mat_samples, 1) -
    centered_kernel_mat_samples * prob_sample_via_base * n) -
    tau * (weight_vec / prob_sample_via_base);

  return s;
}


