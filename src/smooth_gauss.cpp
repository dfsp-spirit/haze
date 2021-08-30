
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <cassert>
//#include <iostream>

/// Perform nearest-neighbor smoothing of the given data, based on mesh adjacency list representation.
/// @param _mesh the input mesh
/// @param _data an R numerical vector with one value per mesh vertex
/// @param _fwhm the FWHM for the Gaussian kernel
/// @param _truncfactor the cutoff factor after which to end the Gaussian neighborhood, in Gaussian standard deviations
RcppExport SEXP smooth_data_gaussian(SEXP _mesh, SEXP _data, SEXP _fwhm, SEXP _truncfactor) {
  return wrap(1);
}

/// Used in C++ code only, not exported.
std::vector<std::vector<float>> gauss_weights(const std::vector<std::vector<int>> geod_neigh_indices, const std::vector<std::vector<float>> geod_neigh_dists, const float gstd) {
  std::vector<std::vector<float>> weights;

  float gvar2 = 2 * (gstd * gstd);
  float f = 1.0 / (sqrt(2 * M_PI) * gstd);

  return(weights);
}

/// Used in C++ code only, not exported.
std::vector<float> spatial_filter(const std::vector<float> data, const std::vector<std::vector<int>> geod_neigh_indices, const std::vector<std::vector<float>> geod_neigh_gauss_weights) {
  std::vector<float> smoothed_data(data.size());
  float smoothed_val;
  for(size_t i=0; i<data.size(); i++) {
    smoothed_val = 0.0;
    for(size_t j=0; j<geod_neigh_indices[i].size(); j++) {
      smoothed_val += data[geod_neigh_indices[i][j]] * geod_neigh_gauss_weights[i][j];
    }
    smoothed_data[i] = smoothed_val;
  }
  return(smoothed_data);
}

/// Used in C++ code only, not exported.
inline float fhwm_to_gstd(const float fwhm) {
  return(fwhm / sqrt(log(256.0)));
}
