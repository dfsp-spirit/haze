#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <vector>


RcppExport SEXP smooth_data(Rcpp::List mesh_adj, Rcpp::NumericVector data, SEXP _num_iter) {

  const int num_iter = Rcpp::as<int >(_num_iter);
  const int num_values = data.size();
  Rcpp::NumericVector source_data;
  Rcpp::NumericVector smoothed_data(num_values);
  Rcpp::IntegerVector vert_neighbors;
  float neigh_sum;
  int num_non_na_values;

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
  return(Rcpp::wrap(smoothed_data));
}
