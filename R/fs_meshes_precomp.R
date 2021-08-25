# Functions that return pre-computed neighborhood information for FreeSurfer brain templates.

#' @title Return pre-computed neighborhood data for specific meshes.
#'
#' @param meshname a text identifier specifying the mesh you want connectivty data for. Currently supported meshes are listed here. 'lh_fsaverage': the left hemisphere of the FreeSurfer 6 fsaverage template. 'rh_fsaverage': the right hemisphere of the FreeSurfer 6 fsaverage template. 'lh_fsaverage6': the left hemisphere of the FreeSurfer 6 fsaverage6 template. 'rh_fsaverage6': the right hemisphere of the FreeSurfer 6 fsaverage6 template.
#'
#' @return list of vectors, the connectivity data as an adjacency list. The outer list has length n, where n is the number of vertices in the graph. The inner lists represent, for each vertex, all of its neighbors.
#'
#' @export
mesh_neigh_pre <- function(meshname) {
}
