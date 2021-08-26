#' @title Ensure the mesh is a tmesh3d instance. Will convert fs.surfaces to one automatically.
#'
#' @param mesh whatever, but hopefully an \code{rgl::tmesh3d} or \code{freesurferformats::fs.surface} instance. Can be a character string, which will be loaded as a surface file if it exists.
#'
#' @return tmesh3d instance, the input or converted from the input.
#'
#' @note This function will stop if the mesh cannot be converted to tmesh3d.
#'
#' @keywords internal
ensure.tmesh3d <- function(mesh) {
  if(is.character(mesh) && length(mesh) == 1L) { # treat as filename}
    if(file.exists(mesh)) {
      mesh = freesurferformats::read.fs.surface(mesh);
    }
  }

  if(freesurferformats::is.fs.surface(mesh)) {
    return(fs.surface.to.tmesh3d(mesh));
  } else if ("mesh3d" %in% class(mesh)) {
    return(mesh);
  } else {
    stop("Cannot convert value in parameter 'mesh' to tmesh3d instance, invalid mesh.");
  }
}

#' @title Get an rgl tmesh3d instance from a brain surface mesh.
#'
#' @param surface an fs.surface instance, as returned by \code{subject.surface} or \code{freesurferformats::read.fs.surface}.
#'
#' @return a tmesh3d instance, see \code{rgl::tmesh3d} for details.
#'
#' @keywords internal
fs.surface.to.tmesh3d <- function(surface) {
  if( ! freesurferformats::is.fs.surface(surface)) {
    stop("Parameter 'surface' must be an instance of freesurferformats::fs.surface.");
  }
  return(rgl::tmesh3d(c(t(surface$vertices)), c(t(surface$faces)), homogeneous=FALSE));
}


#' @title Get an fs.surface brain mesh from an rgl tmesh3d instance.
#'
#' @param tmesh a tmesh3d instance, see \code{rgl::tmesh3d} for details.
#'
#' @return an fs.surface instance, as returned by \code{subject.surface} or \code{freesurferformats::read.fs.surface}.
#'
#' @keywords internal
tmesh3d.to.fs.surface <- function(tmesh) {
  vertices = t(tmesh$vb[1:3,]);
  faces = t(tmesh$it);
  surface = list('vertices'=vertices, 'faces'=faces);
  class(surface) <- c(class(surface), 'fs.surface');
  return(surface);
}


#' @title Check whether parameter is an fs.surface instance.
#'
#' @param surface an fs.surface instance which will be returned as-is, a tmesh3d which will be converted to a surface, or a character string which will be interpreted as a file system path and loaded with \code{freesurferformats::read.fs.surface}. Anything else will stop with an error.
#'
#' @return an fs.surface instance, unless an error occurs.
#'
#' @keywords internal
ensure.fs.surface <- function(surface) {
  if(freesurferformats::is.fs.surface(surface)) {
    return(surface);
  } else if(is.character(surface)) {
    return(freesurferformats::read.fs.surface(surface));
  } else if('mesh3d' %in% class(surface)) {
    return(tmesh3d.to.fs.surface(surface));
  } else {
    stop("Parameter 'surface' must be an fs.surface instance or a character string.");
  }
}
