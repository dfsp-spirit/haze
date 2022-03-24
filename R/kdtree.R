


#' @title Get per-vertex data at vertices closest to the given query coordinates on the mesh.
#'
#' @description Return interpolated values of the mesh data at the given query points. For this to make any sense, the meshes must be spherical with identical radius, and only the vertex positions on the sphere differ.
#'
#' @param query_coordinates nx3 numerical matrix of x,y,z coordinates. These are typically the vertex positions of a second (spherical!) mesh for that you need per-vertex data (e.g., the \code{fsaverage6} mesh).
#'
#' @param mesh a spherical fs.surface instance, see \code{\link[freesurferformats]{read.fs.surface}} to get one.
#'
#' @param pervertex_data numerical vector, the continuous per-vertex data for the vertices of the mesh.
#'
#' @return the per-vertex data for the vertices closest to the query coordinates.
#'
#' @seealso \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_NNInterpolate_kdTree.m}
#'
#' @note The mesh must be spherical, and the query_coordinates must be located on the mesh sphere.
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


#' @title Map input values at points in 3D space onto mesh vertices.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @seealso  \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_linearInterpolate_kdTree.m}
#' @seealso \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_linearInterpolateAux.c}
#' @seealso \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findFace.m}
#'
#' @note The mesh must be spherical, and the query_coordinates must be located on the mesh sphere.
#'
#' @return named list with entries: 'interp_values', the numerical vector of interpolated data at the query_coordinates. 'nearest_vertex_in_face' the nearest vertex in the face that the respective query coordinate falls into, 'nearest_face' the index of the nearest face that the respective query coordinate falls into.
#'
#' @importFrom stats dist
#'
#' @export
linear_interpolate_kdtree <- function(query_coordinates, mesh, pervertex_data) {
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

  tmesh = haze:::ensure.tmesh3d(mesh);
  vertex_neighbors = Rvcg::vcgVertexNeighbors(tmesh); # Compute vertex neighborhood of vertices.
  vertex_faces = Rvcg::vcgVFadj(tmesh);  # Compute all faces the vertices are part of.

  # TODO: decide whether vertex_neighbors and vertex_faces should be list of vectors or a matrix (with NA entries), and enforce/check it here.

  # Get the maximal neighbor count over all mesh vertices. typically 6 or 7 for triangular meshes.
  #max_num_neighbors = ncol(vertex_neighbors); # for matrix, vertices with less will have NA entries at the end of their row.
  max_num_neighbors = max(unlist(lapply(vertex_neighbors, length))); # for list.
  #max_num_vertex_faces = ncol(vertex_faces); # for matrix, vertices which are part of less faces will have NA entries at the end of their row.
  max_num_vertex_faces = max(unlist(lapply(vertex_faces, length))); # for list.
  cat(sprintf("Maximal number of vertex neighbors per vertex (vertex degree) is %d. Maximal number of faces a vertex is part of is %d.\n", max_num_neighbors, max_num_vertex_faces));

  # now, for each query coordinate:
  # -find the face that the coordinate falls into.
  #  * for this, see https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findFaces.h
  # - find the vertex of the face that is closest to the query coordinate (passed in as parameter 'query_coords_closest_vertex')
  # -then we retrieve the 3 vertices of the face and their pervertex_data values.
  # -then we interpolate the value at the query_coordinate between the 3 known values/coordinates.
  tmesh = haze:::ensure.tmesh3d(mesh);
  clost = Rvcg::vcgClost(query_coordinates, tmesh); # or use Rvcg::vcgClostKD(), need to benchmark which is faster.
  nearest_face = clost$faceptr; # vector, for each query coordinate the (index of the) closest face.
  nearest_face_vertices = mesh$faces[nearest_face, ]; # nx3 int matrix, the vertex indices (of verts forming the closest face).

  # one could project the query coordinate onto the triangle plane, then interpolate within a 2D plane.
  # See https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_linearInterp.h for that approach.
  # The akima package could be interesting for R: https://cran.r-project.org/web/packages/akima/akima.pdf


  # TODO: project z exactly onto the xy-plane and rotate coord system to make triangle lie in the xy-plane.
  # see https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane
  nq = nrow(query_coordinates);
  if(nq < 1L) {
    stop("Parameter 'query_coordinates' must contain at least one x,y,z row.");
  }

  interp_values = rep(0.0, nq); # Allocation, gets filled below.

  # The current approach uses inverse distance weighted (IWD) interpolation.
  iwd_beta = 2.0; # beta parameter for IWD
  for(row_idx in seq.int(nq)) {
    qc = query_coordinates[row_idx, ];
    #closest_vertex_in_closest_face_local_idx = which(nearest_face_vertices[row_idx, ] == query_coords_closest_vertex[row_idx]); # 1,2 or 3
    dist_query_to_v1 = stats::dist(rbind(qc, mesh$vertices[nearest_face_vertices[row_idx,],1]));
    dist_query_to_v2 = stats::dist(rbind(qc, mesh$vertices[nearest_face_vertices[row_idx,],2]));
    dist_query_to_v3 = stats::dist(rbind(qc, mesh$vertices[nearest_face_vertices[row_idx,],3]));

    total_dist = sum(dist_query_to_v1, dist_query_to_v2, dist_query_to_v3);
    rel_dist = c(dist_query_to_v1, dist_query_to_v2, dist_query_to_v3) / total_dist;
    # see https://rspatial.org/raster/analysis/4-interpolation.html and
    # https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Inverse-Distance-Weighting/index.html
    weights = rel_dist ** - iwd_beta;
    interp_values[row_idx] = sum(weights*pervertex_data[nearest_face_vertices[row_idx,]])/ sum(weights);
  }


  nearest_vertex_in_face = query_coords_closest_vertex; # this currently is the vertex index (global, in the mesh). Should we compute and return the index in the face (1,2, or 3L) instead?
  return(list("interp_values"=interp_values, "nearest_vertex_in_face"=nearest_vertex_in_face, "nearest_face"=nearest_face));
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
    tmesh = haze:::ensure.tmesh3d(mesh);
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
#' homog = matrix(c(1,2,3,1,1,2,3,2), ncol=4, byrow=TRUE);
#' haze:::homogeneous_to_cartesian(homog);
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

