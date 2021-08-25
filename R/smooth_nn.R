
#' @title Smooth per-vertex data based on mesh.
#'
#' @param surface a mesh, represented as an \code{fs.surface} instance from the \code{freesurferformats} package or a \code{tmesh3d} instance from \code{rgl}.
#'
#' @inheritParams pervertexdata.smoothnn.adj
#'
#' @return numerical vector, the smoothed data.
#'
#' @seealso \code{\link{pervertexdata.smoothnn.adj}} if you already have pre-computed adjacency data for the mesh. Using that data can increase performance considerably, especially if you need to smooth many data sets.
#'
#' @export
pervertexdata.smoothnn <- function(surface, data, num_iter) {
  k = 1L;
  if(! freesurferformats::is.fs.surface(surface)) {
    stop("Parameter 'surface' must be an fs.surface instance.");
  }
  if(nrow(surface$vertices) != length(data)) {
    stop("Number of surface vertices must match data length.");
  }

  if(requireNamespace("Rvcg", quietly = TRUE)) {
    tmesh = ensure.tmesh3d(surface);
    adj = Rvcg::vcgVertexNeighbors(tmesh, vi = NULL, numstep = k, include_self = TRUE);
  } else {
    stop("The 'Rvcg' package must be installed to use this functionality.");
  }

  nv = nrow(surface$vertices);
  if(length(data) != nv) {
    stop("Data and vertex count mismatch");
  }
  return(pervertexdata.smoothnn.adj(adj, data, num_iter));
}


#' @title Smooth per-vertex data based on mesh adjacency information.
#'
#' @param mesh_adj list of vectors of integers, the adjacency list representation of the mesh. One can use the pre-computed adjacency for some special meshes, see \code{\link{mesh_neigh_pre}}. Data for vertices should include the vertex itself.
#'
#' @param data numerical vector of per-vertex-data for the mesh, one value per vertex. Data values of \code{NA} will be ignored, allowing you to mask parts of the data.
#'
#' @param num_iter positive integer, number of smoothing iterations.
#'
#' @return numerical vector, the smoothed data.
#'
#' @seealso \code{\link{pervertexdata.smoothnn}} if you have a mesh and still need the connectivity to be computed.
#'
#' @export
pervertexdata.smoothnn.adj <- function(mesh_adj, data, num_iter) {
  nv = length(data);
  data_smoothed = rep(NA, nv);

  for(iteration in seq.int(num_iter)) {
    if(iteration == 1L) {
      source_data = data;
    } else {
      source_data = data_smoothed;
    }
    for(vidx in seq.int(nv)) {
      if(is.na(data[vidx])) { next; }
      #neigh = c(mesh_adj[[vidx]], vidx);
      neigh = mesh_adj[[vidx]];
      data_smoothed[vidx] = mean(source_data[neigh], na.rm = TRUE);
    }
  }
  return(data_smoothed);
}
