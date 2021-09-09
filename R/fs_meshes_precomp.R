# Functions that return pre-computed neighborhood information for FreeSurfer brain templates.

#' @title Return pre-computed neighborhood data for specific meshes.
#'
#' @param meshname a text identifier specifying the mesh you want connectivity data for. Currently supported meshes are listed here. 'lh_fsaverage': the left hemisphere of the FreeSurfer 6 fsaverage template. 'rh_fsaverage': the right hemisphere of the FreeSurfer 6 fsaverage template. 'lh_fsaverage6': the left hemisphere of the FreeSurfer 6 fsaverage6 template. 'rh_fsaverage6': the right hemisphere of the FreeSurfer 6 fsaverage6 template.
#'
#' @return list of vectors, the connectivity data as an adjacency list. The outer list has length n, where n is the number of vertices in the graph. The inner lists represent, for each vertex, all of its neighbors.
#'
#' @export
mesh.neigh.pre <- function(meshname) {
  vvfile = NULL;
  if(meshname == "lh_fsaverage") {
    vvfile = system.file("extdata", "fsaverage_lh_white_meshdist_edge_1.vv", package = "haze", mustWork = TRUE);
  } else if(meshname == "rh_fsaverage") {
    vvfile = system.file("extdata", "fsaverage_rh_white_meshdist_edge_1.vv", package = "haze", mustWork = TRUE);
  } else if(meshname == "lh_fsaverage6") {
    vvfile = system.file("extdata", "fsaverage6_lh_white_meshdist_edge_1.vv", package = "haze", mustWork = TRUE);
  } else if(meshname == "rh_fsaverage6") {
    vvfile = system.file("extdata", "fsaverage6_rh_white_meshdist_edge_1.vv", package = "haze", mustWork = TRUE);
  } else {
    stop("Invalid mesh name. See function help for supported ones.");
  }
  return(read.vv(vvfile));
}
