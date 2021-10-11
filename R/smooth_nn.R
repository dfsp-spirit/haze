
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
  if(nrow(surface$vertices) != length(pvdata)) {
    stop("Number of surface vertices must match pvdata length.");
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
#' @param pvdata numerical vector of per-vertex-data for the mesh, one value per vertex. Data values of \code{NA} will be ignored, allowing you to mask parts of the data. If you pass an nxm matrix, the n rows will be treated as (independent) overlays that should be smoothed in parallel. To set the number of cores to use for parallel processing, set the 'mc_cores' option like this: \code{options("mc.cores"=22L);} before calling this function.
#'
#' @param num_iter positive integer, number of smoothing iterations.
#'
#' @param method character string, one of 'C++' or 'R'. The C++ version is much faster (about 50 times faster on our test machine), and there is little reason to ever use the R version in production code, so just ignore this.
#'
#' @param silent logical, whether to suppress output messages.
#'
#' @return numerical vector, the smoothed data.
#'
#' @seealso \code{\link{pervertexdata.smoothnn}} if you have a mesh and still need the connectivity to be computed.
#'
#' @examples
#' \dontrun{
#' mesh = rgl::tetrahedron3d();
#' mesh_adj = mesh.adj(mesh, k = 1L);
#' pvd = rnorm(nrow(mesh$vb), mean = 5.0, sd = 1.0);
#' pvd_smoothed = pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 30L);
#' }
#'
#' @export
pervertexdata.smoothnn.adj <- function(mesh_adj, pvdata, num_iter, method="C++", silent = getOption("haze.silent", default = TRUE)) {

  if(is.matrix(pvdata)) {
    if(nrow(pvdata) == 1L | ncol(pvdata) == 1L) {
      pvdata = as.vector(pvdata);
    } else {
      return(pervertexdata.smoothnn.adj.mat(mesh_adj, pvdata, num_iter, method=method));
    }
  }


  if(method == "C++") {
    return(pervertexdata.smoothnn.adj.cpp(mesh_adj, pvdata, num_iter));
  }
  if(method != "R") {
    stop("Parameter 'method' must be one of 'C++' or 'R'.");
  }

  nv = length(pvdata);
  if(! silent) {
    cat(sprintf("Smoothing %d iterations over the %d pvdata values in R.\n", num_iter, nv));
  }
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
  if(length(num_iter) != 1L) {
    stop("Parameter 'num_iter' must be a scalar, positive integer.");
  }
  if(! is.numeric(num_iter)) {
    stop("Parameter 'num_iter' must be a scalar, positive integer.");
  }
  if(num_iter < 1L) {
    stop("Parameter 'num_iter' must be a scalar, positive integer.");
  }
  # Adapt the vertex indices for C++ (from 1-based R indexing to 0-based C++ indexing).
  mesh_adj = lapply(mesh_adj, function(x){x-1L});

  if(! is.vector(pvdata)) {
    stop("Parameter 'pvdata' must be a numeric vector.");
  }
  return(.Call("smooth_data", mesh_adj, pvdata, num_iter));
}


#' @title Smooth matrix in parallel
#'
#' @inheritParams pervertexdata.smoothnn.adj
#'
#' @param pvdata numerical matrix of per-vertex-data for the mesh. Each row is an overlay with one value per vertex. Rows are independent and smoothed in parallel. Data values of \code{NA} will be ignored, allowing you to mask parts of the data.
#'
#' @examples
#' \dontrun{
#' mesh = rgl::tetrahedron3d();
#' mesh_adj = mesh.adj(mesh, k = 1L);
#' pvd1 = rnorm(nrow(mesh$vb), mean = 1.0, sd = 1.0);
#' pvd2 = rnorm(nrow(mesh$vb), mean = 4.0, sd = 1.0);
#' pvd3 = rnorm(nrow(mesh$vb), mean = 7.0, sd = 1.0);
#' pvd_matrix = rbind(pvd1, pvd2, pvd3);
#' # or: pvd_matrix = matrix(data=rnorm(nrow(mesh$vb)*5, mean=5.0, sd=1.0), nrow=5);
#' pvd_smoothed = haze:::pervertexdata.smoothnn.adj.mat(mesh_adj, pvd_matrix, num_iter = 30L);
#' }
#'
#' @note There is no need to call this manually, it will be called by \code{pervertexdata.smooth.adj} when its input is a matrix.
#'
#' @note To set the number of cores to use, set the 'mc_cores' option like this: \code{options("mc.cores"=22L);} before calling this function.
#'
#' @keywords internal
pervertexdata.smoothnn.adj.mat <- function(mesh_adj, pvdata, num_iter, method = "C++", silent = getOption("haze.silent", default = FALSE)) {
  if(! is.matrix(pvdata)) {
    stop("Parameter 'pvdata' must be a matrix.");
  }
  num_cores = getOption("mc.cores", default = 2L);

  if(! silent) {
    cat(sprintf("Handling %d vertex overlays with %d values each in parallel on %d CPU cores. (Use 'options(\"mc.cores\"=N)' to request N cores).\n", nrow(pvdata), ncol(pvdata), num_cores));
  }

  res_list = parallel::mclapply( 1L:nrow(pvdata), mc.cores = num_cores, function(row_idx){ pervertexdata.smoothnn.adj(mesh_adj, pvdata[row_idx, ], num_iter=num_iter, method = method) } );
  return(t(as.matrix(data.frame(res_list))));
}

