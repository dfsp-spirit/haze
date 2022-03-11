


#' @title Map input values at points in 3D space onto mesh vertices.
#'
#' @description Map input values at points in 3D space onto mesh vertices. For this to make any sense, the meshes must be very similar in shape. Typically both meshes are spherical meshes with identical radius, and only the vertex positions on the sphere differ.
#'
#' @param input_coordinates nx3 numerical matrix of x,y,z coordinates.
#'
#' @param mesh an fs.surface instance, see \code{\link[freesurferformats]{read.fs.surface}} to get one.
#'
#' @param input_values numerical vector, the continuous per-vertex data for the point at the input_coordinates.
#'
#' @return the interpolated per-vertex data for the mesh given in the \code{mesh} parameter.
#'
#' @seealso \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_NNInterpolate_kdTree.m}
#'
#' @export
nn_interpolate_kdtree <- function(input_coordinates, mesh, input_values) {

}


#' @title Map input values at points in 3D space onto mesh vertices.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @seealso  \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_linearInterpolate_kdTree.m}
#'
#' @export
linear_interpolate_kdtree <- function(input_coordinates, mesh, input_values) {

}


#' @title Find nearest mesh vertex for a coordinate using kdtree.
#'
#' @inheritParams nn_interpolate_kdtree
#'
#' @return integer vector, the \code{n} vertex indices which are closest to the \code{nx3} matrix of input_coordinates.
#'
#' @seealso  \code{https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findNV_kdTree.m}
#'
#' @export
find_nv_kdtree <- function(input_coordinates, mesh) {

}



