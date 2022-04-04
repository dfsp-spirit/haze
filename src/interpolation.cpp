#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <cassert>
//#include <iostream>

/// Interpolate values within mesh triangles.
/// @param _query_coordinates nx3 double matrix of query coordinates (points for which you need the interpolated per-vertex data).
/// @param _mesh_vertices mx3 matrix of mesh vertex coordinates (the mesh that holds the _pervertex_data)
/// @param _nearest_face_vertices nx3 int matrix, for each of the n query_coordinates, the face (given as 3 vertex indices) into which the query coordinate falls on the mesh
/// @param _pervertex_data double vector with m values, the per-vertex data for the mesh, that should be interpolated at the query_coordinates
/// @param _iwd_beta, scalar double, the beta for the inverse distance weighting algorithm. typically something between 1.0 and 2.0
RcppExport SEXP interp_tris_c(SEXP _query_coordinates, SEXP _mesh_vertices, SEXP _nearest_face_vertices, SEXP _pervertex_data, SEXP _iwd_beta) {

  Rcpp::List query_coordinates(_query_coordinates);         // nx3 double matrix
  Rcpp::List mesh_vertices(_mesh_vertices);                 // mx3 double matrix
  Rcpp::List nearest_face_vertices(_nearest_face_vertices); // nx3 int matrix

  Rcpp::NumericVector pervertex_data(_pervertex_data);
  Rcpp::NumericVector interp_data(_pervertex_data);
  const double iwd_beta = Rcpp::as<double>(_iwd_beta);


  return(interp_data);
}
