


#' @title Return interpolated values of the mesh data at the given query points.
#'
#' @description Return interpolated values of the mesh data at the given query points. For this to make any sense, the meshes must be spherical with identical radius, and only the vertex positions on the sphere differ.
#'
#' @param input_coordinates nx3 numerical matrix of x,y,z coordinates.
#'
#' @param mesh an fs.surface instance, see \code{\link[freesurferformats]{read.fs.surface}} to get one.
#'
#' @param input_values numerical vector, the continuous per-vertex data for the point at the input_coordinates.
#'
#' @return the interpolated per-vertex data for the mesh given in the \code{mesh} parameter.
#'
#' @seealso \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_NNInterpolate_kdTree.m}
#'
#' @keywords internal
# @export
nn_interpolate_kdtree <- function(input_coordinates, mesh, input_values) {
  warning("nn_interpolate_kdtree: NOT IMPLEMENTED YET, returning fake data.");
  return(rnorm(nrow(mesh$vertices), 5.0, 1.0));
}


#' @title Map input values at points in 3D space onto mesh vertices.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @seealso  \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_linearInterpolate_kdTree.m}
#'
#' @keywords internal
# @export
linear_interpolate_kdtree <- function(input_coordinates, mesh, input_values) {
  warning("linear_interpolate_kdtree: NOT IMPLEMENTED YET, returning fake data.");
  return(rnorm(nrow(mesh$vertices), 5.0, 1.0));
}


#' @title Find nearest mesh vertex for a coordinate using kdtree.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @return integer vector, the \code{n} vertex indices which are closest to the \code{nx3} matrix of input_coordinates.
#'
#' @seealso  \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findNV_kdTree.m}
#'
#' @keywords internal
#'
#' @importFrom parallel detectCores
#' @importFrom Rvcg vcgCreateKDtreeFromBarycenters vcgClostOnKDtreeFromBarycenters
# export
find_nv_kdtree <- function(input_coordinates, mesh) {
  warning("find_nv_kdtree: NOT IMPLEMENTED YET, returning fake data.");
  return(seq.int(nrow(input_coordinates)));

  mesh = ensure.fs.surface(mesh);
  if(! freesurferformats::is.fs.surface(mesh)) {
    stop("Parameter 'mesh' must be an fs.surface instance.");
  }
  if(is.matrix(input_coordinates)) {
    if(ncol(input_coordinates) != 3L) {
      stop("The matrix 'input_coordinates' must have exactly 3 columns.");
    }
    if(! is.numeric(input_coordinates)) {
      stop("The matrix 'input_coordinates' must be numeric.");
    }
    tmesh = ensure.tmesh3d(mesh);
    kdtreeBary = Rvcg::vcgCreateKDtreeFromBarycenters(tmesh);
    #nv_mesh = Rvcg::vcgClostOnKDtreeFromBarycenters(kdtreeBary, input_coordinates, threads=parallel::detectCores())
    #closest_vertex_coords = nv_mesh$vb;
    #closest_vertex_distances = nv_mesh$quality;
    res = Rvcg::vcgSearchKDtree(kdtreeBary, input_coordinates, k=1L, threads = parallel::detectCores());
  } else {
    stop("Parameter 'input_coordinates' must be a matrix.");
  }
}



