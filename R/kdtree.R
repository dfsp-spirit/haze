


#' @title Get per-vertex data at vertices closest to the given query coordinates on the mesh.
#'
#' @description Return interpolated values of the mesh data at the given query points. For this to make any sense, the meshes must be spherical with identical radius, and only the vertex positions on the sphere differ.
#'
#' @param query_coordinates nx3 numerical matrix of x,y,z coordinates.
#'
#' @param mesh an fs.surface instance, see \code{\link[freesurferformats]{read.fs.surface}} to get one.
#'
#' @param pervertex_data numerical vector, the continuous per-vertex data for the vertices of the mesh.
#'
#' @return the per-vertex data for the vertices closest to the query coordinates.
#'
#' @seealso \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_NNInterpolate_kdTree.m}
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
#'
#' @export
linear_interpolate_kdtree <- function(query_coordinates, mesh, pervertex_data) {
  mesh = ensure.fs.surface(mesh);
  if(length(pervertex_data) != nrow(mesh$vertices)) {
    warning(sprintf("The 'pervertex_data' is for %d vertices, but the mesh has %d. Expected identical values.\n",length(pervertex_data), nrow(mesh$vertices)));
  }
  res = find_nv_kdtree(query_coordinates, mesh);

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
#' @importFrom parallel detectCores
#' @importFrom Rvcg vcgCreateKDtreeFromBarycenters vcgSearchKDtree
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

