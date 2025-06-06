// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// centered_kernel_matrix
arma::mat centered_kernel_matrix(const arma::vec& eval_points_1, const arma::vec& eval_points_2, const arma::vec& centering_grid, double hurst_coef);
RcppExport SEXP _kefV1_centered_kernel_matrix(SEXP eval_points_1SEXP, SEXP eval_points_2SEXP, SEXP centering_gridSEXP, SEXP hurst_coefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eval_points_1(eval_points_1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eval_points_2(eval_points_2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type centering_grid(centering_gridSEXP);
    Rcpp::traits::input_parameter< double >::type hurst_coef(hurst_coefSEXP);
    rcpp_result_gen = Rcpp::wrap(centered_kernel_matrix(eval_points_1, eval_points_2, centering_grid, hurst_coef));
    return rcpp_result_gen;
END_RCPP
}
// centered_kernel_matrix_hd
arma::mat centered_kernel_matrix_hd(const arma::mat& eval_points_1, const arma::mat& eval_points_2, const arma::mat& centering_grid, double hurst_coef);
RcppExport SEXP _kefV1_centered_kernel_matrix_hd(SEXP eval_points_1SEXP, SEXP eval_points_2SEXP, SEXP centering_gridSEXP, SEXP hurst_coefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eval_points_1(eval_points_1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type eval_points_2(eval_points_2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type centering_grid(centering_gridSEXP);
    Rcpp::traits::input_parameter< double >::type hurst_coef(hurst_coefSEXP);
    rcpp_result_gen = Rcpp::wrap(centered_kernel_matrix_hd(eval_points_1, eval_points_2, centering_grid, hurst_coef));
    return rcpp_result_gen;
END_RCPP
}
// unnormalised_density_samples
arma::vec unnormalised_density_samples(const arma::mat& centered_kernel_mat_samples, double lambda, const arma::vec& weight_vec);
RcppExport SEXP _kefV1_unnormalised_density_samples(SEXP centered_kernel_mat_samplesSEXP, SEXP lambdaSEXP, SEXP weight_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type centered_kernel_mat_samples(centered_kernel_mat_samplesSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight_vec(weight_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(unnormalised_density_samples(centered_kernel_mat_samples, lambda, weight_vec));
    return rcpp_result_gen;
END_RCPP
}
// unnormalised_density_grids
arma::vec unnormalised_density_grids(const arma::mat& centered_kernel_mat_grids, const arma::vec& centered_kernel_self_grids, double lambda, const arma::vec& weight_vec);
RcppExport SEXP _kefV1_unnormalised_density_grids(SEXP centered_kernel_mat_gridsSEXP, SEXP centered_kernel_self_gridsSEXP, SEXP lambdaSEXP, SEXP weight_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type centered_kernel_mat_grids(centered_kernel_mat_gridsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type centered_kernel_self_grids(centered_kernel_self_gridsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight_vec(weight_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(unnormalised_density_grids(centered_kernel_mat_grids, centered_kernel_self_grids, lambda, weight_vec));
    return rcpp_result_gen;
END_RCPP
}
// get_dens_wo_grid
arma::vec get_dens_wo_grid(const arma::mat& centered_kernel_mat_samples, const arma::vec& samples, const arma::vec& base_measure_weights, double dimension, double lambda, const arma::vec& weight_vec);
RcppExport SEXP _kefV1_get_dens_wo_grid(SEXP centered_kernel_mat_samplesSEXP, SEXP samplesSEXP, SEXP base_measure_weightsSEXP, SEXP dimensionSEXP, SEXP lambdaSEXP, SEXP weight_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type centered_kernel_mat_samples(centered_kernel_mat_samplesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type base_measure_weights(base_measure_weightsSEXP);
    Rcpp::traits::input_parameter< double >::type dimension(dimensionSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight_vec(weight_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dens_wo_grid(centered_kernel_mat_samples, samples, base_measure_weights, dimension, lambda, weight_vec));
    return rcpp_result_gen;
END_RCPP
}
// get_dens
Rcpp::List get_dens(const arma::mat& centered_kernel_mat_samples, const arma::mat& centered_kernel_mat_grids, const arma::vec& centered_kernel_self_grids, const arma::vec& samples, const arma::vec& grids, const arma::vec& base_measure_weights, double dimension, double lambda, const arma::vec& weight_vec);
RcppExport SEXP _kefV1_get_dens(SEXP centered_kernel_mat_samplesSEXP, SEXP centered_kernel_mat_gridsSEXP, SEXP centered_kernel_self_gridsSEXP, SEXP samplesSEXP, SEXP gridsSEXP, SEXP base_measure_weightsSEXP, SEXP dimensionSEXP, SEXP lambdaSEXP, SEXP weight_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type centered_kernel_mat_samples(centered_kernel_mat_samplesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type centered_kernel_mat_grids(centered_kernel_mat_gridsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type centered_kernel_self_grids(centered_kernel_self_gridsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type grids(gridsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type base_measure_weights(base_measure_weightsSEXP);
    Rcpp::traits::input_parameter< double >::type dimension(dimensionSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight_vec(weight_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dens(centered_kernel_mat_samples, centered_kernel_mat_grids, centered_kernel_self_grids, samples, grids, base_measure_weights, dimension, lambda, weight_vec));
    return rcpp_result_gen;
END_RCPP
}
// get_middle_points_grid
arma::vec get_middle_points_grid(double min, const arma::vec& samples, double max);
RcppExport SEXP _kefV1_get_middle_points_grid(SEXP minSEXP, SEXP samplesSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type min(minSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(get_middle_points_grid(min, samples, max));
    return rcpp_result_gen;
END_RCPP
}
// get_s_function
arma::vec get_s_function(const arma::vec& weight_vec, double lambda, double tau, const arma::mat& centered_kernel_mat_samples, const arma::vec& samples, const arma::vec& base_measure_weights, double dimension);
RcppExport SEXP _kefV1_get_s_function(SEXP weight_vecSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP centered_kernel_mat_samplesSEXP, SEXP samplesSEXP, SEXP base_measure_weightsSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type weight_vec(weight_vecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type centered_kernel_mat_samples(centered_kernel_mat_samplesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type base_measure_weights(base_measure_weightsSEXP);
    Rcpp::traits::input_parameter< double >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(get_s_function(weight_vec, lambda, tau, centered_kernel_mat_samples, samples, base_measure_weights, dimension));
    return rcpp_result_gen;
END_RCPP
}
// generate_voronoi
Rcpp::List generate_voronoi(Rcpp::NumericMatrix points, double x_min, double x_max, double y_min, double y_max);
RcppExport SEXP _kefV1_generate_voronoi(SEXP pointsSEXP, SEXP x_minSEXP, SEXP x_maxSEXP, SEXP y_minSEXP, SEXP y_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< double >::type y_min(y_minSEXP);
    Rcpp::traits::input_parameter< double >::type y_max(y_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_voronoi(points, x_min, x_max, y_min, y_max));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kefV1_centered_kernel_matrix", (DL_FUNC) &_kefV1_centered_kernel_matrix, 4},
    {"_kefV1_centered_kernel_matrix_hd", (DL_FUNC) &_kefV1_centered_kernel_matrix_hd, 4},
    {"_kefV1_unnormalised_density_samples", (DL_FUNC) &_kefV1_unnormalised_density_samples, 3},
    {"_kefV1_unnormalised_density_grids", (DL_FUNC) &_kefV1_unnormalised_density_grids, 4},
    {"_kefV1_get_dens_wo_grid", (DL_FUNC) &_kefV1_get_dens_wo_grid, 6},
    {"_kefV1_get_dens", (DL_FUNC) &_kefV1_get_dens, 9},
    {"_kefV1_get_middle_points_grid", (DL_FUNC) &_kefV1_get_middle_points_grid, 3},
    {"_kefV1_get_s_function", (DL_FUNC) &_kefV1_get_s_function, 7},
    {"_kefV1_generate_voronoi", (DL_FUNC) &_kefV1_generate_voronoi, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_kefV1(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
