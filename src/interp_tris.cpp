#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <cassert>
//#include <iostream>

/// Interpolate values within mesh triangles.
/// @param _query_coordinates nx3 matrix of query coordinates.
RcppExport SEXP interp_tris_c(SEXP _query_coordinates, SEXP _mesh_vertices, SEXP _nearest_face_vertices, SEXP _pervertex_data, SEXP _iwd_beta) {

  Rcpp::List query_coordinates(_query_coordinates);         // nx3 double matrix
  Rcpp::List mesh_vertices(_mesh_vertices);                 // mx3 double matrix
  Rcpp::List nearest_face_vertices(_nearest_face_vertices); // nx3 int matrix

  Rcpp::NumericVector pervertex_data(_pervertex_data);
  Rcpp::NumericVector interp_data(_pervertex_data);
  const double iwd_beta = Rcpp::as<double>(_iwd_beta);


  return(interp_data);
}
