# Convert between FreeSurfer Gaussian smoothing (FWHM) and nearest-neighbor
# smoothing iteration counts.
#
# FreeSurfer emulates Gaussian smoothing on the cortical surface by iterating
# a cheap local operation: averaging each vertex with its direct (1-ring) edge
# neighbors. By the central limit theorem on a graph, repeating this many
# times produces an effective kernel that converges to a Gaussian, so no
# expensive geodesic distance computation is ever needed. The number of
# iterations that corresponds to a given Gaussian FWHM was calibrated
# empirically by FreeSurfer (see MRISfwhm2niters() and MRISniters2fwhm() in
# utils/mrisutils.cpp of the FreeSurfer source).


#' @title Compute the mean vertex area of a triangular mesh.
#'
#' @param surface a mesh, represented as an \code{fs.surface} instance from the \code{freesurferformats} package, a \code{tmesh3d} instance from \code{rgl}, or a character string representing the path of a surface file to load.
#' @param avg_vertex_area numeric scalar, the mean vertex area in squared mesh units (e.g., mm^2). If \code{NULL} (the default), it is computed from \code{surface}.
#'
#' @return numeric scalar, the mean (average) vertex area of the mesh in squared mesh units.
#'
#' @keywords internal
get.avg.vertex.area <- function(surface = NULL, avg_vertex_area = NULL) {
  if(! is.null(avg_vertex_area)) {
    if(! is.numeric(avg_vertex_area) || length(avg_vertex_area) != 1L || avg_vertex_area <= 0) {
      stop("Parameter 'avg_vertex_area' must be a positive number.");
    }
    return(avg_vertex_area);
  }
  if(is.null(surface)) {
    stop("Either 'surface' or 'avg_vertex_area' must be given.");
  }
  tmesh = ensure.tmesh3d(surface);
  total_area = sum(Rvcg::vcgArea(tmesh));   # per-face areas, so the sum is the total surface area.
  return(total_area / ncol(tmesh$vb));
}


#' @title Convert FreeSurfer Gaussian smoothing FWHM to number of NN smoothing iterations.
#'
#' @description Convert a Gaussian smoothing kernel width, given as the full width at half maximum (FWHM) in millimeters as used by FreeSurfer, into the number of nearest-neighbor (NN) smoothing iterations to use with \code{\link{pervertexdata.smoothnn}} (or \code{\link{pervertexdata.smoothnn.adj}}) to achieve a similar amount of smoothing.
#'
#' @param fwhm numeric scalar, the desired full width at half maximum (FWHM) of the Gaussian smoothing kernel in millimeters, as used by FreeSurfer's \code{mris_fwhm --fwhm <fwhm>} or \code{mri_surf2surf --fwhm <fwhm>}.
#' @param surface a mesh, represented as an \code{fs.surface} instance from the \code{freesurferformats} package, a \code{tmesh3d} instance from \code{rgl}, or a character string representing the path of a surface file to load. Used to compute the mean vertex area. Ignored if \code{avg_vertex_area} is given.
#' @param avg_vertex_area numeric scalar, the mean vertex area of the mesh in square millimeters (squared mesh units). If \code{NULL} (the default), it is computed from \code{surface}. See details for how FreeSurfer's own numbers differ.
#' @param k positive integer, the k-ring neighborhood size to use for NN smoothing. The FreeSurfer-based calibration below is only valid for the 1-ring neighborhood (\code{k=1}); a warning is issued otherwise.
#'
#' @return integer scalar, the number of NN smoothing iterations (the \code{num_iter} parameter of \code{\link{pervertexdata.smoothnn}}) that approximates Gaussian smoothing with the requested FWHM.
#'
#' @details
#' True Gaussian smoothing on a triangular mesh would require computing geodesic distances between (many) vertex pairs, which is very slow. FreeSurfer avoids this by emulating the Gaussian with many iterations of nearest-neighbor (1-ring) averaging: the whole trick is that \emph{repeated local averaging converges to a Gaussian kernel} (central limit theorem on a graph), so no geodesic distances are needed at all. Each iteration is a cheap local operation (average a vertex with its direct edge neighbors), which is exactly what \code{\link{pervertexdata.smoothnn}} implements.
#'
#' The relationship between the effective Gaussian standard deviation \code{gstd} after \code{n} iterations and the mean vertex area \code{avg_vertex_area} was calibrated empirically by FreeSurfer in \code{MRISfwhm2niters()} (see utils/mrisutils.cpp in the FreeSurfer source):
#' \preformatted{
#'   gstd   = fwhm / sqrt(log(256))                     # FWHM -> Gaussian sigma (log(256) = 8*ln(2))
#'   niters = round(1.14 * (4 * pi * gstd^2) / (7 * avg_vertex_area))
#' }
#' where \code{7} approximates the average number of neighbors in the 1-ring of a triangulated cortical surface, and \code{1.14} is an empirical fudge factor that FreeSurfer fitted so that the measured FWHM of the smoothed output matches the requested one. This function is the inverse of \code{\link{niters2fwhm}}.
#'
#' \strong{Note on matching FreeSurfer's own outputs:} FreeSurfer normalizes the mean vertex area by the \emph{group-average} surface area stored in its template surface files (\code{group_avg_surface_area}, about 82,000 mm^2 per hemisphere for fsaverage), which is larger than the raw geometric area of the fsaverage meshes (about 65,000 mm^2 for the packaged fsaverage lh white surface). FreeSurfer therefore reports somewhat lower iteration counts than you get when \code{avg_vertex_area} is computed directly from these meshes (e.g., for FWHM = 10 mm on fsaverage: 74 iterations in FreeSurfer vs. 92 from the raw mesh area). To reproduce FreeSurfer's exact numbers, pass \code{avg_vertex_area = group_avg_surface_area / nvertices} (about 0.5015 mm^2 for the fsaverage templates).
#'
#' @examples
#' \dontrun{
#' mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
#' niters = fwhm2niters(10, surface = mesh);      # 92 iterations, approximating ~10 mm FWHM
#' smoothed = pervertexdata.smoothnn(mesh, data, num_iter = niters);
#' }
#'
#' @seealso \code{\link{niters2fwhm}} for the inverse conversion, \code{\link{pervertexdata.smoothnn}} and \code{\link{pervertexdata.smoothnn.adj}} for the actual NN smoothing.
#'
#' @export
fwhm2niters <- function(fwhm, surface = NULL, avg_vertex_area = NULL, k = 1L) {
  if(! is.numeric(fwhm) || length(fwhm) != 1L || fwhm <= 0) {
    stop("Parameter 'fwhm' must be a positive number.");
  }
  if(k != 1L) {
    warning("The FreeSurfer-based FWHM to iteration conversion is only calibrated for k=1 (1-ring) neighborhoods; results for other k are not meaningful.");
  }
  avg_vertex_area = get.avg.vertex.area(surface, avg_vertex_area);
  gstd = fwhm / sqrt(log(256.0));
  niters = floor((1.14 * 4.0 * pi * gstd * gstd) / (7.0 * avg_vertex_area) + 0.5);
  niters = max(niters, 1L);
  return(as.integer(niters));
}


#' @title Convert number of NN smoothing iterations to a FreeSurfer Gaussian smoothing FWHM.
#'
#' @description Inverse of \code{\link{fwhm2niters}}: convert a number of nearest-neighbor (NN) smoothing iterations into the equivalent Gaussian smoothing kernel width, given as the full width at half maximum (FWHM) in millimeters.
#'
#' @inheritParams fwhm2niters
#' @param niters positive integer, the number of NN smoothing iterations (the \code{num_iter} parameter of \code{\link{pervertexdata.smoothnn}}).
#'
#' @return numeric scalar, the FWHM in millimeters of the Gaussian smoothing kernel that the given number of NN smoothing iterations approximates.
#'
#' @details
#' This is the inverse of \code{\link{fwhm2niters}}, using the same empirically calibrated FreeSurfer formula (utils/mrisutils.cpp in the FreeSurfer source):
#' \preformatted{
#'   gstd = sqrt((7 * avg_vertex_area * niters) / (1.14 * 4 * pi))
#'   fwhm = gstd * sqrt(log(256))
#' }
#' It is useful to check what amount of Gaussian smoothing (FWHM) a given number of NN smoothing iterations corresponds to, e.g., to report the equivalent FWHM in a paper or to compare settings across methods. The calibration is only valid for the 1-ring neighborhood (\code{k=1}).
#'
#' @examples
#' \dontrun{
#' mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
#' fwhm = niters2fwhm(92, surface = mesh);       # ~10 mm FWHM
#' }
#'
#' @seealso \code{\link{fwhm2niters}} for the forward conversion, \code{\link{pervertexdata.smoothnn}} and \code{\link{pervertexdata.smoothnn.adj}} for the actual NN smoothing.
#'
#' @export
niters2fwhm <- function(niters, surface = NULL, avg_vertex_area = NULL, k = 1L) {
  if(! is.numeric(niters) || length(niters) != 1L || niters < 1L) {
    stop("Parameter 'niters' must be a positive integer.");
  }
  if(k != 1L) {
    warning("The FreeSurfer-based iteration to FWHM conversion is only calibrated for k=1 (1-ring) neighborhoods; results for other k are not meaningful.");
  }
  avg_vertex_area = get.avg.vertex.area(surface, avg_vertex_area);
  gstd = sqrt((7.0 * avg_vertex_area * niters) / (1.14 * 4.0 * pi));
  fwhm = gstd * sqrt(log(256.0));
  return(fwhm);
}
