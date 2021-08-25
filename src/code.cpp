#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <iostream>


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

  std::cout << "Smoothing " << std::to_string(num_iter) << " iterations over the " << std::to_string(num_values) << " data values in C++.\n";

  for (int i = 0; i < num_iter; i++){
    if(i == 0) {
      source_data = data;
    } else {
      source_data = smoothed_data;
    }
    for (int j = 0; j < num_values; j++){
      vert_neighbors = mesh_adj[j];
      neigh_sum = 0.0;
      num_non_na_values = 0;
      for(int k=0; k<vert_neighbors.size(); k++) {
        // Could check for NA here
        num_non_na_values++;
        neigh_sum += data[vert_neighbors[k]];
      }
      smoothed_data[j] = neigh_sum / (float)num_non_na_values;
    }
  }
  return(smoothed_data);
}
