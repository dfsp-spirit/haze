


#' @title Get per-vertex data at vertices closest to the given query coordinates on the mesh.
#'
#' @description Return per-vertex data at the vertices closest to the given query points.
#'
#' @param query_coordinates nx3 numerical matrix of x,y,z coordinates. These are typically the vertex positions of a second (spherical!) mesh for that you need per-vertex data (e.g., the \code{fsaverage6} mesh).
#'
#' @param mesh fs.surface instance, see \code{\link[freesurferformats]{read.fs.surface}} or \code{\link[fsbrain]{subject.surface}} to get one, or turn an \code{rgl} \code{tmesh} into one with \code{\link[fsbrain]{tmesh3d.to.fs.surface}}.
#'
#' @param pervertex_data numerical vector, the continuous per-vertex data for the vertices of the mesh.
#'
#' @return the per-vertex data for the vertices closest to the query coordinates.
#'
#' @export
nn_interpolate_kdtree <- function(query_coordinates, mesh, pervertex_data) {
  mesh = ensure.fs.surface(mesh);
  if(length(pervertex_data) != nrow(mesh$vertices)) {
    warning(sprintf("The 'pervertex_data' is for %d vertices, but the mesh has %d. Expected identical values.\n",length(pervertex_data), nrow(mesh$vertices)));
  }
  res = find_nv_kdtree(query_coordinates, mesh);
  return(pervertex_data[res$index]);
}


#' @title Interpolate per-vertex data at the query points. Or map per-vertex data between subjects.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @description This method uses inverse distance weight interpolation within a triangle. First, the face of the \code{mesh} that the \code{query_coordinate} falls into is determined. Then results in 3 vertices with respective per-vertex data, and a query coordinate. We then compute the distance to all 3 vertices, and perform inverse distance weight interpolation with a beta setting defined by parameter \code{iwd_beta}.
#'
#' @param iwd_beta scalar double, the \code{beta} parameter for the inverse distance weight interpolation with the triangle. See details.
#'
#' @param ... ignore, passed on to internal function \code{interp_tris}.
#'
#' @note The mesh must be spherical, and the \code{query_coordinates} must also be located on the mesh sphere.
#'
#' @return named list with entries: 'interp_values', the numerical vector of interpolated data at the query_coordinates. 'nearest_vertex_in_face' the nearest vertex in the face that the respective query coordinate falls into, 'nearest_face' the index of the nearest face that the respective query coordinate falls into.
#'
#' @importFrom stats dist
#'
#' @export
linear_interpolate_kdtree <- function(query_coordinates, mesh, pervertex_data, iwd_beta = 2.0, ...) {
  mesh = ensure.fs.surface(mesh);
  if(length(pervertex_data) != nrow(mesh$vertices)) {
    warning(sprintf("The 'pervertex_data' is for %d vertices, but the mesh has %d. Expected identical values.\n",length(pervertex_data), nrow(mesh$vertices)));
  }
  if(ncol(query_coordinates) != 3L) {
    stop("Parameter query_coordinates must be nx3 matrix of Cartesian 3d coordinates.");
  }
  if(length(pervertex_data) != nrow(mesh$vertices)) {
    warning(sprintf("The 'pervertex_data' is for %d vertices, but the mesh has %d. Expected identical values.\n", length(pervertex_data), nrow(mesh$vertices)));
  }

  res = find_nv_kdtree(query_coordinates, mesh);
  query_coords_closest_vertex = res$index;

  if(length(query_coords_closest_vertex) != nrow(query_coordinates)) {
    stop(sprintf("Number of query_coordinates (%d) must match number of closest vertices for the query coordinates (%d).\n", nrow(query_coordinates), length(query_coords_closest_vertex)));
  }

  #tmesh = ensure.tmesh3d(mesh);
  #vertex_neighbors = Rvcg::vcgVertexNeighbors(tmesh); # Compute vertex neighborhood of vertices.
  #vertex_faces = Rvcg::vcgVFadj(tmesh);  # Compute all faces the vertices are part of.

  # TODO: decide whether vertex_neighbors and vertex_faces should be list of vectors or a matrix (with NA entries), and enforce/check it here.
  # Get the maximal neighbor count over all mesh vertices. typically 6 or 7 for triangular meshes.
  #max_num_neighbors = ncol(vertex_neighbors); # for matrix, vertices with less will have NA entries at the end of their row.
  #max_num_neighbors = max(unlist(lapply(vertex_neighbors, length))); # for list.
  #max_num_vertex_faces = ncol(vertex_faces); # for matrix, vertices which are part of less faces will have NA entries at the end of their row.
  #max_num_vertex_faces = max(unlist(lapply(vertex_faces, length))); # for list.
  #cat(sprintf("Maximal number of vertex neighbors per vertex (vertex degree) is %d. Maximal number of faces a vertex is part of is %d.\n", max_num_neighbors, max_num_vertex_faces));

  # now, for each query coordinate:
  # -find the face that the coordinate falls into.
  #  * for this, see https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findFaces.h
  # - find the vertex of the face that is closest to the query coordinate (passed in as parameter 'query_coords_closest_vertex')
  # -then we retrieve the 3 vertices of the face and their pervertex_data values.
  # -then we interpolate the value at the query_coordinate between the 3 known values/coordinates.
  tmesh = ensure.tmesh3d(mesh);
  clost = Rvcg::vcgClost(query_coordinates, tmesh); # or use Rvcg::vcgClostKD(), need to benchmark which is faster.
  nearest_face = clost$faceptr; # vector, for each query coordinate the (index of the) closest face.
  nearest_face_vertices = mesh$faces[nearest_face, ]; # nx3 int matrix, the vertex indices (of verts forming the closest face).

  # One could project the query coordinate onto the triangle plane, then interpolate within a 2D plane.
  # This may give better results if the query points are not exactly on the sphere, but in our case, that should not be needed.
  # See https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_linearInterp.h for that approach.
  # To project z exactly onto the xy-plane, see also https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane.

  nq = nrow(query_coordinates);
  if(nq < 1L) {
    stop("Parameter 'query_coordinates' must contain at least one x,y,z row.");
  }

  interp_values = interp_tris(query_coordinates, mesh$vertices, nearest_face_vertices, pervertex_data, iwd_beta = iwd_beta, ...);

  nearest_vertex_in_face = query_coords_closest_vertex; # this currently is the vertex index (global, in the mesh). Should we compute and return the index in the face (1,2, or 3L) instead?
  return(list("interp_values"=interp_values, "nearest_vertex_in_face"=nearest_vertex_in_face, "nearest_face"=nearest_face));
}


#' @title Interpolate values within mesh triangles in R.
#'
#' @inheritParams linear_interpolate_kdtree
#'
#' @param mesh_vertices \code{nx3} matrix of x,y,z cartesian coordinates for the n mesh vertices (of the mesh that has the 'pervertex_data').
#'
#' @return vector with one entry per query coordinate (i.e., \code{nrow(query_coordinates)}).
#'
#' @param cpp logical, whether to use the much faster \code{C++} version. Leave alone unless you know what you are doing. Changing this to \code{FALSE} will make this a lot slower.
#'
#' @keywords internal
interp_tris <- function(query_coordinates, mesh_vertices, nearest_face_vertices, pervertex_data, iwd_beta = 2.0, cpp=TRUE) {
  if(cpp) {
    return(interp_tris_cpp_wrapper(query_coordinates, mesh_vertices, nearest_face_vertices, pervertex_data, iwd_beta = iwd_beta));
  }
  nq = nrow(query_coordinates);
  interp_values = rep(0.0, nq); # Allocation, gets filled below.

  if(length(pervertex_data) != nrow(mesh_vertices)) {
    stop(sprintf("Mesh has %d vertices, but pervertex_data has %d entries, counts must match.\n", nrow(mesh_vertices), length(pervertex_data)));
  }

  if(ncol(query_coordinates) != 3L) {
    stop("Parameter 'query_coordinates' must be an nx3 matrix.");
  }

  # The current approach uses inverse distance weighted (IWD) interpolation.
  # beta parameter for IWD
  for(row_idx in seq.int(nq)) {
    qc = query_coordinates[row_idx, ];
    #closest_vertex_in_closest_face_local_idx = which(nearest_face_vertices[row_idx, ] == query_coords_closest_vertex[row_idx]); # 1,2 or 3
    dist_query_to_v1 = stats::dist(rbind(qc, mesh_vertices[nearest_face_vertices[row_idx,],1]));
    dist_query_to_v2 = stats::dist(rbind(qc, mesh_vertices[nearest_face_vertices[row_idx,],2]));
    dist_query_to_v3 = stats::dist(rbind(qc, mesh_vertices[nearest_face_vertices[row_idx,],3]));

    total_dist = sum(dist_query_to_v1, dist_query_to_v2, dist_query_to_v3);
    rel_dist = c(dist_query_to_v1, dist_query_to_v2, dist_query_to_v3) / total_dist;
    # see https://rspatial.org/raster/analysis/4-interpolation.html and
    # https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Inverse-Distance-Weighting/index.html
    weights = rel_dist ** - iwd_beta;
    interp_values[row_idx] = sum(weights*pervertex_data[nearest_face_vertices[row_idx,]]) / sum(weights);
  }
  return(interp_values);
}


#' @title Interpolate values within mesh triangles in C++.
#'
#' @inheritParams linear_interpolate_kdtree
#'
#' @param mesh_vertices \code{nx3} matrix of x,y,z Cartesian coordinates for the n mesh vertices (of the mesh that has the 'pervertex_data').
#'
#' @return vector with one entry per query coordinate (i.e., \code{nrow(query_coordinates)}).
#'
#' @note No validation of input data is done, do this yourself before passing values to this function to avoid hard crashes (which happen when junk is passed to \code{C++}).
#'
#' @keywords internal
interp_tris_cpp_wrapper <- function(query_coordinates, mesh_vertices, nearest_face_vertices, pervertex_data, iwd_beta = 2.0) {

  # adjust indexing for C++
  #nearest_face_vertices = nearest_face_vertices + 1L;

  return(.Call("interp_tris_c", query_coordinates, mesh_vertices, nearest_face_vertices, pervertex_data, iwd_beta));
}


#' @title Find nearest mesh vertex for query coordinates using kdtree.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @param threads integer, the number of threads to run in parallel.
#'
#' @return named list with keys 'index' and 'distance'. 'index': integer vector, the \code{n} vertex indices which are closest to the \code{nx3} matrix of query_coordinates. 1-based indices for R are returned. 'distance': double vector, the distances to the respective vertices in the 'index' key.
#'
#' @seealso  \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findNV_kdTree.m}
#'
#' @keywords internal
#'
#' @note @note The mesh must be spherical, and the query_coordinates must be located on the mesh sphere.
#'
#' @importFrom parallel detectCores
#' @importFrom Rvcg vcgCreateKDtree vcgSearchKDtree
#' @export
find_nv_kdtree <- function(query_coordinates, mesh, threads = parallel::detectCores()) {

  mesh = ensure.fs.surface(mesh);
  if(! freesurferformats::is.fs.surface(mesh)) {
    stop("Parameter 'mesh' must be an fs.surface instance.");
  }
  if(is.matrix(query_coordinates)) {
    if(ncol(query_coordinates) != 3L) {
      stop("The matrix 'query_coordinates' must have exactly 3 columns.");
    }
    if(! is.numeric(query_coordinates)) {
      stop("The matrix 'query_coordinates' must be numeric.");
    }
    tmesh = ensure.tmesh3d(mesh);
    kdtree = Rvcg::vcgCreateKDtree(tmesh);
    vcg_res = Rvcg::vcgSearchKDtree(kdtree, query_coordinates, k=1L, threads = threads);
    res = list('index'=vcg_res$index, 'distance'=vcg_res$distance);
    return(res);
  } else {
    stop("Parameter 'query_coordinates' must be a matrix.");
  }
}


#' @title Convert homogeneous coordinates to Cartesian coordinates.
#'
#' @param homog nx4 numeric matrix of input coordinates
#'
#' @return nx3 matrix of Cartesian coordinates
#'
#' @examples
#' \dontrun{
#' homog = matrix(c(1,2,3,1,1,2,3,2), ncol=4, byrow=TRUE);
#' haze:::homogeneous_to_cartesian(homog);
#' }
#'
#' @keywords internal
homogeneous_to_cartesian <- function(homog) {
  if(! is.matrix(homog)) {
    stop("Parameter 'homog' must be a matrix.");
  }
  if(ncol(homog) != 4L) {
    stop("Parameter 'homog' must be a matrix with 4 columns.");
  }
  return(homog[,1:3] / homog[,4]);
}

