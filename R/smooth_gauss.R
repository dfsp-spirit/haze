

#' @title Compute gstd from FWHM.
#' @keywords internal
fwhm.to.gstd <- function(fwhm) { fwhm / sqrt(log(256.0)); }

#' @title Perform Gaussian smoothing of per-vertex data on a mesh given geodesic distance data.
#'
#' @inheritParams surf.sphere.gaussianweights
#'
#' @param data numerical vector of per-vertex data for the surface
#'
#' @param fwhm positive double, the full width at half maximum for the Gaussian kernel
#'
#' @param truncfactor the factor after how many stddevs to truncate the Gaussian kernel
#'
#' @return the smoothed data
#'
#' @export
pervertexdata.smoothgaussian <- function(surface, fwhm, geodesic_neigh = NULL, truncfactor = 3.5) {
  if(is.null(geodesic_neigh)) {
    geodesic_neigh = geodesic.neighborhoods(surface, fwhm, truncfactor);
  }
  gaussian_weights = mesh.gaussianweights(geodesic_neigh, fwhm.to.gstd(fwhm));
  smoothed_data = mesh.spatialfilter(data, geodesic_neigh, gaussian_weights);
  return(smoothed_data);
}

#' @title Compute geodesic neighborhoods for all mesh vertices.
#'
#' @note This can get quite large, and may fill RAM, but allows re-use of the neighborhoods. The alternative would be to compute on the fly, but then one has to re-compute this very expensive operation each time.
#'
#' @export
geodesic.neighborhoods <- function(surface, fwhm, truncfactor = 3.5) {
  maxdist = truncfactor * fwhm.to.gstd(fwhm);
  res = list("neigh_idx"=list(), "neigh_dist"=list());
  num_verts = ...; # FIX THIS
  for(vidx in seq.int(num_verts)) {
    single_vert_neigh = geod.vert.neighborhood(surface, vertex, max_distance=maxdist);
    message("FIX THIS: the following var names are incorrect");
    res$neigh_idx[[length(res$neigh_idx)+1L]] = single_vert_neigh$neig_idx;
    res$neigh_dist[[length(res$neigh_dist)+1L]] = single_vert_neigh$neig_dist;
  }
  return(res);
}

#' @title Compute all vertices within given geodesic distance on the mesh.
#'
#' @param mesh an instance of \code{rgl::tmesh3d} or \code{freesurferformats::fs.surface}.
#'
#' @param vertex positive integer (or vector of the latter), the index of the source vertex in the mesh. If a vector, the neighborhoods for all vertices will be computed separately.
#'
#' @param max_distance double, the neighborhood size. All mesh vertices in geodesic distance smaller than / up to this distance will be returned.
#'
#' @param include_max logical, whether the max_distance value is inclusive.
#'
#' @param return_distances logical, whether to compute the 'distances' entry in the returned list. Doing so is a little bit slower, so it can be turned off if not needed.
#'
#' @return named list with the following entries: 'vertices': integer vector, the indices of all vertices in the neigborhood. 'distances': double vector, the distances to the respective vertices (unless 'return_distances' is FALSE).
#'
#' @note This function uses the pseudo-geodesic distance along the mesh edges.
#'
#' @examples
#' \dontrun{
#'   sjd = fsaverage.path(TRUE);
#'   surface = subject.surface(sjd, 'fsaverage', surface = "white", hemi = "lh");
#'   res = geod.vert.neighborhood(surface, 12345L, max_distance = 10.0);
#'   res$vertices;
#' }
#'
#' @export
geod.vert.neighborhood <- function(mesh, vertex, max_distance=5.0, include_max = TRUE, return_distances = TRUE) {
  mesh = ensure.tmesh3d(mesh);
  if(requireNamespace("Rvcg", quietly = TRUE)) {
    neighborhood = vertex;
    neighborhood_distances = rep(0.0, length(vertex));
    for (v in vertex) {
      geodesic_dists_to_vertex = geodesic.dists.to.vertex(mesh, v);
      if(include_max) {
        this_neighborhood_indices = which(geodesic_dists_to_vertex <= max_distance);
      } else {
        this_neighborhood_indices = which(geodesic_dists_to_vertex < max_distance);
      }
      neighborhood = c(neighborhood, this_neighborhood_indices);
      neighborhood_distances = c(neighborhood_distances, geodesic_dists_to_vertex[this_neighborhood_indices]);
    }

    if(! return_distances) {
      return(list("vertices" = unique(neighborhood)));
    } else {
      # Make neighborhood unique, and also remove the corresponding duplicated distances (so we cannot simply use base::unique).
      verts_unique = c();
      dists_unique = c();
      for(neigh_idx in seq_len(length(neighborhood))) {
        vert_idx = neighborhood[neigh_idx];
        if(! (vert_idx %in% verts_unique)) {
          verts_unique = c(verts_unique, vert_idx);
          dists_unique = c(dists_unique, neighborhood_distances[neigh_idx]);
        }
      }
      return(list('vertices' = verts_unique, 'distances' = dists_unique));
    }
  } else {
    stop("The 'Rvcg' package must be installed to use this functionality.");
  }
}


#' @title Simple internal wrapper around \code{Rvcg::vcgDijkstra} with function check.
#'
#' @param mesh a tmesh3d instance.
#'
#' @param v positive integer, a vertex index in the mesh.
#'
#' @return double vector with length equal to num vertices in the mesh, the geodesic distances from all other vertices to the query vertex \code{v}.
#'
#' @keywords internal
geodesic.dists.to.vertex <- function(mesh, v) {
  if(! exists('vcgDijkstra', where=asNamespace('Rvcg'), mode='function')) {
    stop("Your Rvcg version does not export the vcgDijkstra function, which means it is too old. You need to install Rvcg from GitHub for this this functionality to be available. Try 'devtools::install_github('zarquon42b/Rvcg')'.");
  } else {
    return(Rvcg::vcgDijkstra(mesh, v));
  }
}


#' @title Apply spatial filter to surface data.
#'
#' @param source_data numerical vector, per-vertex data for a surface.
#'
#' @param sphere_dists named list with 3 entries, as returned by \code{surf.sphere.dist}
#'
#' @param gaussian_weight list of double vectors, the Gaussian weights for all neighbors of the respective vertex. As returned by \code{surf.sphere.gaussianweights}.
#'
#' @return numerical vector, the spatially filtered per-vertex data.
#'
#' @keywords internal
mesh.spatialfilter <- function(source_data, sphere_dists, gaussian_weights) {
  smoothed_data = rep(NA, length(source_data));
  if(! is.list(sphere_dists$neigh)) {
    stop("Parameter 'sphere_dists' member 'neigh' must be a list.");
  }
  if(! is.list(gaussian_weights)) {
    stop("Parameter 'gaussian_weights' must be a list.");
  }
  for(vidx in seq_along(source_data)) {
    smoothed_data[vidx] = sum(source_data[sphere_dists$neigh[[vidx]]] * gaussian_weights[[vidx]]);
  }
  return(smoothed_data);
}


#' @title Compute Gaussian weights
#'
#' @inheritParams surf.sphere.dist
#'
#' @param sphere_dists list of vectors, as returned by surf.sphere.dist
#'
#' @param gstd double, Gaussian standard deviation, can be computed from the FWHM as \code{gstd = fwhm / sqrt(log(256.0))}.
#'
#' @keywords internal
mesh.gaussianweights <- function(sphere_dists, gstd) {

  gvar2 = 2 * (gstd * gstd); # twice the variance
  f = 1.0 / (sqrt(2 * pi) * gstd);

  num_neighbors = unlist(lapply(sphere_dists$neigh, length));
  weights = list();
  for(vidx in seq(nv)) {
    gsum = 0.0;
    vert_weights = rep(NA, num_neighbors[vidx]);
    local_idx = 1L;

    for(neigh_vidx in sphere_dists$neigh[[vidx]]) {
      d = sphere_dists$neigh_dist_surface[[vidx]][local_idx];
      g = f * exp(-(d * d) / (gvar2));
      vert_weights[local_idx] = g;
      gsum = gsum + g;
      local_idx = local_idx + 1L;
    }
    vert_weights = vert_weights / gsum; # Normalize
    weights[[vidx]] = vert_weights;
    #cat(sprintf("Vertex %d has %d neighbors: %s. weights=%s\n", vidx, num_neighbors[vidx], paste(sphere_dists$neigh[[vidx]], collapse = " "), paste(weights[[vidx]], collapse = " ")));
  }
  return(weights);
}
