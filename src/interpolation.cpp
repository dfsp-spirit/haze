#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <cassert>
//#include <iostream>

inline double euclid(Rcpp::NumericVector v1, Rcpp::NumericVector v2) {
  double dx = v1(0) - v2(0);
  double dy = v1(1) - v2(1);
  double dz = v1(2) - v2(2);
  return(sqrt(dx * dx + dy * dy + dz * dz));
}

/// Interpolate values within mesh triangles.
/// @param _query_coordinates nx3 double matrix of query coordinates (points for which you need the interpolated per-vertex data).
/// @param _mesh_vertices mx3 matrix of mesh vertex coordinates (the mesh that holds the _pervertex_data)
/// @param _nearest_face_vertices nx3 int matrix, for each of the n query_coordinates, the face (given as 3 vertex indices) into which the query coordinate falls on the mesh
/// @param _pervertex_data double vector with m values, the per-vertex data for the mesh, that should be interpolated at the query_coordinates
/// @param _iwd_beta, scalar double, the beta for the inverse distance weighting algorithm. typically something between 1.0 and 2.0
RcppExport SEXP interp_tris_c(SEXP _query_coordinates, SEXP _mesh_vertices, SEXP _nearest_face_vertices, SEXP _pervertex_data, SEXP _iwd_beta) {

  Rcpp::NumericMatrix query_coordinates(_query_coordinates);         // nx3 double matrix
  Rcpp::NumericMatrix mesh_vertices(_mesh_vertices);                 // mx3 double matrix
  Rcpp::IntegerMatrix nearest_face_vertices(_nearest_face_vertices); // nx3 int matrix

  int nq = query_coordinates.nrow();

  Rcpp::NumericVector pervertex_data(_pervertex_data);
  Rcpp::NumericVector interp_data(nq);
  const double iwd_beta = Rcpp::as<double>(_iwd_beta);

  Rcpp::IntegerVector nearest_face_vertex_indices;
  Rcpp::NumericVector qc, rel_dist, weights, v1_vertex_coords, v2_vertex_coords, v3_vertex_coords;
  double dist_query_to_v1, dist_query_to_v2, dist_query_to_v3, total_dist, total_weights, part_v1, part_v2, part_v3;
  for(int row_idx = 0; row_idx < nq; row_idx++) {
    qc = query_coordinates( row_idx , _ );

    nearest_face_vertex_indices = {nearest_face_vertices(row_idx,0), nearest_face_vertices(row_idx,1), nearest_face_vertices(row_idx,2)};
    v1_vertex_coords = {mesh_vertices(nearest_face_vertex_indices[0], 0), mesh_vertices(nearest_face_vertex_indices[0], 1), mesh_vertices(nearest_face_vertex_indices[0], 2)};
    v2_vertex_coords = {mesh_vertices(nearest_face_vertex_indices[1], 0), mesh_vertices(nearest_face_vertex_indices[1], 1), mesh_vertices(nearest_face_vertex_indices[1], 2)};
    v3_vertex_coords = {mesh_vertices(nearest_face_vertex_indices[2], 0), mesh_vertices(nearest_face_vertex_indices[2], 1), mesh_vertices(nearest_face_vertex_indices[2], 2)};
    dist_query_to_v1 = euclid(qc, v1_vertex_coords);
    dist_query_to_v2 = euclid(qc, v2_vertex_coords);
    dist_query_to_v3 = euclid(qc, v3_vertex_coords);

    total_dist = dist_query_to_v1 + dist_query_to_v2 + dist_query_to_v3;
    rel_dist = {dist_query_to_v1/total_dist, dist_query_to_v2/total_dist, dist_query_to_v3/total_dist};
    weights = { std::pow(rel_dist[0], - iwd_beta), std::pow(rel_dist[1], - iwd_beta), std::pow(rel_dist[2], - iwd_beta)  };
    total_weights = weights[0] + weights[1] + weights[2];
    part_v1 = weights[0] * pervertex_data[nearest_face_vertices(row_idx, 0)];
    part_v2 = weights[1] * pervertex_data[nearest_face_vertices(row_idx, 1)];
    part_v3 = weights[2] * pervertex_data[nearest_face_vertices(row_idx, 2)];
    interp_data[row_idx] = (part_v1 + part_v2 + part_v3) / total_weights;
  }
  return(interp_data);
}




