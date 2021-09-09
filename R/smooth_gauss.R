


#' @title Smooth per-vertex data based on a mesh using true Gaussian smoothing.
#'
#' @description This function performs Gaussian smoothing of vertices in the geodesic neighborhood of all vertices. This is computationally quite expensive and requires large amounts of memory for large meshes.
#'
#' @param surface a mesh, represented as an \code{fs.surface} instance from the \code{freesurferformats} package or a \code{tmesh3d} instance from \code{rgl}, or a character string representing the path of a mesh to load with \code{freesurferformats::read.fs.surface}.
#'
#' @inheritParams pervertexdata.smoothnn
#'
#' @param fwhm scalar double smoothing kernel full width at half max
#'
#' @param trunc_factor scalar double, truncation factor for Gaussian neighborhood, in Gaussian standard deviations.
#'
#' @return numerical vector, the smoothed data.
#'
#' @seealso \code{\link{pervertexdata.smoothnn}} can be used to approximate Gaussian smoothing with several iterations of nearest neighbor smoothing, and is a lot faster for large meshes.
#'
#' @examples
#' \dontrun{
#' mesh = rgl::tetrahedron3d();
#' pvd = rnorm(nrow(mes2$vb), mean = 5.0, sd = 1.0);
#' pvd_smoothed = pervertexdata.smoothgauss(mesh, pvd, fwhm=5.0);
#' }
#'
#' @export
pervertexdata.smoothgauss <- function(surface, data, fwhm, trunc_factor=3.5) {
  if(requireNamespace("Rvcg", quietly = TRUE)) {
    return(Rvcg::vcgSmoothPVD(ensure.tmesh3d(surface), data, fwhm = fwhm, trunc_factor = trunc_factor));
  } else {
    stop("The 'Rvcg' package is required to use this functionality. Use 'remotes::install_github('dfsp-spirit/Rvcg', ref='smooth_pervertex_data')' to obtain the correct package version.");
  }
}

