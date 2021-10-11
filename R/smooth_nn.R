
#' @title Smooth per-vertex data based on mesh.
#'
#' @param surface a mesh, represented as an \code{fs.surface} instance from the \code{freesurferformats} package or a \code{tmesh3d} instance from \code{rgl}, or a character string representing the path of a mesh to load with \code{freesurferformats::read.fs.surface}.
#'
#' @inheritParams pervertexdata.smoothnn.adj
#'
#' @inheritParams mesh.adj
#'
#' @return numerical vector, the smoothed data.
#'
#' @seealso \code{\link{pervertexdata.smoothnn.adj}} if you already have pre-computed adjacency data for the mesh. Using that data can increase performance considerably, especially if you need to smooth many data sets.
#'
#' @examples
#' \dontrun{
#' mesh = rgl::tetrahedron3d();
#' pvd = rnorm(nrow(mes2$vb), mean = 5.0, sd = 1.0);
#' pvd_smoothed = pervertexdata.smoothnn(mesh, pvd, num_iter = 30L);
#' }
#'
#' @export
pervertexdata.smoothnn <- function(surface, pvdata, num_iter, k=1L, method="C++") {
  surface = ensure.fs.surface(surface);
  if(! freesurferformats::is.fs.surface(surface)) {
    stop("Parameter 'surface' must be an fs.surface instance.");
  }
  if(nrow(surface$vertices) != length(data)) {
    stop("Number of surface vertices must match data length.");
  }

  tmesh = ensure.tmesh3d(surface);
  adj = mesh.adj(tmesh, k = k);

  nv = nrow(surface$vertices);
  if(length(pvdata) != nv) {
    stop("Data and vertex count mismatch");
  }
  return(pervertexdata.smoothnn.adj(adj, pvdata, num_iter, method=method));
}


#' @title Compute vertex neighborhoods for a mesh.
#'
#' @inheritParams pervertexdata.smoothnn.adj
#'
#' @inheritParams pervertexdata.smoothnn
#'
#' @param k scalar positive integer, the k value for the k-ring neighborhood. For k=1, this function computes the adjacency list representation of the graph (where the neighbors include the vertex itself).
#'
#' @return list of integer vectors, the neighborhood data
#'
#' @examples
#' \dontrun{
#' mesh = rgl::tetrahedron3d();
#' mesh_adj = mesh.adj(mesh, k = 1L);
#' }
#'
#' @export
mesh.adj <- function(surface, k = 1L) {
  if(requireNamespace("Rvcg", quietly = TRUE)) {
    tmesh = ensure.tmesh3d(surface);
    adj = Rvcg::vcgVertexNeighbors(tmesh, vi = NULL, numstep = k, include_self = TRUE);
  } else {
    stop("The 'Rvcg' package must be installed to use this functionality.");
  }
}


#' @title Smooth per-vertex data using nearest-neighbor smoothing based on mesh adjacency information.
#'
#' @param mesh_adj list of vectors of integers, the adjacency list representation of the mesh. One can use the pre-computed adjacency for some special meshes, see \code{\link{mesh.neigh.pre}}. Data for vertices should include the vertex itself.
#'
#' @param pvdata numerical vector of per-vertex-data for the mesh, one value per vertex. Data values of \code{NA} will be ignored, allowing you to mask parts of the data.
#'
#' @param num_iter positive integer, number of smoothing iterations.
#'
#' @param method character string, one of 'C++' or 'R'. The C++ version is much faster (about 50 times faster on our test machine), and there is little reason to ever use the R version in production code, so just ignore this.
#'
#' @return numerical vector, the smoothed data.
#'
#' @seealso \code{\link{pervertexdata.smoothnn}} if you have a mesh and still need the connectivity to be computed.
#'
#' @examples
#' \dontrun{
#' mesh = rgl::tetrahedron3d();
#' mesh_adj = mesh.adj(mesh, k = 1L);
#' pvd = rnorm(nrow(mes2$vb), mean = 5.0, sd = 1.0);
#' pvd_smoothed = pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 30L);
#' }
#'
#' @export
pervertexdata.smoothnn.adj <- function(mesh_adj, pvdata, num_iter, method="C++") {
  if(method == "C++") {
    return(pervertexdata.smoothnn.adj.cpp(mesh_adj, pvdata, num_iter));
  }
  if(method != "R") {
    stop("Parameter 'method' must be one of 'C++' or 'R'.");
  }

  nv = length(pvdata);
  #cat(sprintf("Smoothing %d iterations over the %d pvdata values in R.\n", num_iter, nv));
  pvdata_smoothed = rep(NA, nv);

  for(iteration in seq.int(num_iter)) {
    if(iteration == 1L) {
      source_data = pvdata;
    } else {
      source_data = pvdata_smoothed;
    }
    for(vidx in seq.int(nv)) {
      if(is.na(pvdata[vidx])) { next; }
      neigh = mesh_adj[[vidx]];
      pvdata_smoothed[vidx] = mean(source_data[neigh], na.rm = TRUE);
    }
  }
  return(pvdata_smoothed);
}


#' @title Smooth data, C++ version.
#'
#' @inheritParams pervertexdata.smoothnn.adj
#'
#' @return numerical vector, the smoothed data.
#'
#' @keywords internal
pervertexdata.smoothnn.adj.cpp <- function(mesh_adj, pvdata, num_iter) {
  if(! is.list(mesh_adj)) {
    stop("Parameter 'mesh_adj' must be a list of integer vectors.");
  }
  if(is.matrix(pvdata)) {
    if(nrow(pvdata) == 1L) {
      pvdata = as.vector(pvdata);
    }
  }
  if(! is.vector(pvdata)) {
    stop("Parameter 'pvdata' must be a numeric vector.");
  }
  if(length(num_iter) != 1L) {
    stop("Parameter 'num_iter' must be a scalar, positive integer.");
  }

  # Adapt the vertex indices for C++ (from 1-based R indexing to 0-based C++ indexing).
  mesh_adj = lapply(mesh_adj, function(x){x-1L});

  return(.Call("smooth_data", mesh_adj, pvdata, num_iter));
}
