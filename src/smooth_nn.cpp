#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <cassert>
//#include <iostream>

/// Perform nearest-neighbor smoothing of the given data, based on mesh adjacency list representation.
/// @param _mesh_adj an R list of integer vectors, the adjacency list representation of the mesh. This must not be the 1-ring neighborhood, you can pass any neighborhood definition you want.
/// @param _data an R numerical vector with one value per mesh vertex
/// @param _num_iter an R scalar positive integer, the number of smoothing iterations to perform
RcppExport SEXP smooth_data(SEXP _mesh_adj, SEXP _data, SEXP _num_iter) {

  Rcpp::List mesh_adj(_mesh_adj);
  Rcpp::NumericVector data(_data);
  const int num_iter = Rcpp::as<int>(_num_iter);
  const int num_values = data.size();
  Rcpp::NumericVector source_data;
  Rcpp::NumericVector smoothed_data(num_values);
  Rcpp::IntegerVector vert_neighbors;
  float neigh_sum;
  int num_non_na_values;

  //std::cout << "Smoothing " << std::to_string(num_iter) << " iterations over the " << std::to_string(num_values) << " data values in C++.\n";

  assert(num_iter > 0);
  assert(mesh_adj.size() == num_values);

  int num_skip_na_self = 0;
  int num_skip_na_neighbors = 0;

  for (int i = 0; i < num_iter; i++){
    if(i == 0) {
      source_data = data;
    } else {
      source_data = smoothed_data;
    }
    for (int j = 0; j < num_values; j++){
      if (NumericVector::is_na(source_data[j])) {
        num_skip_na_self++;
        smoothed_data[j] = NA_REAL;
        continue;
      }
      vert_neighbors = mesh_adj[j];
      neigh_sum = 0.0;
      num_non_na_values = 0;
      for(int k=0; k<vert_neighbors.size(); k++) {
        if (NumericVector::is_na(source_data[vert_neighbors[k]])) {
          num_skip_na_neighbors++;
          continue; // Ignore NA values.
        }
        num_non_na_values++;
        neigh_sum += source_data[vert_neighbors[k]];
      }
      smoothed_data[j] = neigh_sum / (float)num_non_na_values;
    }
  }
  //std::cout << "Ignored " << std::to_string(num_skip_na_self) << " NA vertices. Ignored " << std::to_string(num_skip_na_neighbors) << " NA neighbors of non-NA vertices.\n";
  return(smoothed_data);
}
