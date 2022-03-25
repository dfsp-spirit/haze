

#' @title Create a submesh including only the given vertices.
#'
#' @param surface_mesh an \code{fs.surface} instance, the original mesh. See \code{\link[freesurferformats]{read.fs.surface}} or \code{\link[fsbrain]{subject.surface}} to get one. Can also be an rgl tmesh, see \code{\link[rgl]{tmesh3d}}.
#'
#' @param old_vertex_indices_to_use integer vector, the vertex indices of the 'surface_mesh' that should be used to construct the new sub mesh.
#'
#' @param ret_mappings whether to return the vertex mappings. If \code{TRUE}, the return value becomes a list with entries 'submesh', 'vmap_full_to_submesh', and 'vmap_submesh_to_full'.
#'
#' @return the new mesh, made up of the given 'old_vertex_indices_to_use' and all (complete) faces that exist between the query vertices in the source mesh. But see 'ret_mapping' parameter.
#'
#' @examples
#' \dontrun{
#' if(requireNamespace("fsbrain, quietly=T")) {
#' sjd = fsbrain::fsaverage.path(T);
#' sj = "fsaverage";
#' mesh = fsbrain::subject.surface(sjd, sj, hemi="lh");
#' lab = fsbrain::subject.label(sjd, sj, "cortex", hemi = "lh");
#' sm = submesh.vertex(mesh, lab);
#' fsbrain::vis.fs.surface(mesh); # show the full mesh.
#' fsbrain::vis.fs.surface(sm);   # show only the cortex.
#' }}
#'
#' @export
#' @importFrom stats complete.cases
submesh.vertex <- function(surface_mesh, old_vertex_indices_to_use, ret_mappings=FALSE) {

  surface_mesh = ensure.fs.surface(surface_mesh);

  if(! is.vector(old_vertex_indices_to_use)) {
    stop("Argument 'old_vertex_indices_to_use' must be a vector.");
  }
  old_vertex_indices_to_use = sort(as.integer(old_vertex_indices_to_use)); # sort is essential! The vertex indices in 'old_vertex_indices_to_use' may not be sorted,
  # and the order of the 'new_vertices' will be wrong then (vertices will be ordered incorrectly, and thus faces will be broken).

  nv_old = nrow(surface_mesh$vertices);
  if(min(old_vertex_indices_to_use) < 1L | max(old_vertex_indices_to_use) > nv_old) {
    stop(sprintf("Invalid 'old_vertex_indices_to_use' parameter: must be integer vector containing values >=1 and <=num_verts(surface_mesh), which is %d.\n", nv_old));
  }

  #nv_new = length(old_vertex_indices_to_use);

  vert_mapping = rep(NA, nv_old); # position/index is old vertex, value is new vertex. old ones not in new mesh receive value of NA.

  # Create a map from the old vertex indices to the new ones. Needed to construct faces later.
  mapped_new_vertex_index = 1L;
  vertex_is_retained = rep(FALSE, nv_old);
  vertex_is_retained[old_vertex_indices_to_use] = TRUE;
  for(old_vert_idx in seq(nv_old)) {
    if(vertex_is_retained[old_vert_idx]) {
      vert_mapping[old_vert_idx] = mapped_new_vertex_index;
      mapped_new_vertex_index = mapped_new_vertex_index + 1L;
    } # no 'else' needed, the rest stays at NA.
  }

  # Use the subset of the old vertices (simply grab coords).
  new_vertices = surface_mesh$vertices[old_vertex_indices_to_use, ];

  # Now for the faces.
  nf_old = nrow(surface_mesh$faces);
  new_faces = matrix(rep(NA, (nf_old*3L)), ncol=3L, nrow=nf_old); #over-allocate and remove invalid ones later.

  new_face_idx = 0L;
  for(old_face_idx in seq(nf_old)) {
    new_face_idx = new_face_idx + 1L;
    old_face_indices = surface_mesh$faces[old_face_idx, ];
    new_face_indices = vert_mapping[old_face_indices];
    new_faces[new_face_idx, ] = new_face_indices;
  }

  df = data.frame(new_faces);
  new_faces = data.matrix(df[stats::complete.cases(df),]); # remove all faces containing an NA vertex

  new_mesh = list('vertices'=new_vertices, 'faces'=new_faces); #, 'vert_mapping'=vert_mapping); # the sub mesh
  class(new_mesh) = c(class(new_mesh), 'fs.surface');

  if(ret_mappings) {
    nnv = nrow(new_vertices); # nnv = number of new vertices.
    rev_mapping = rep(-1L, nnv);
    for(map_idx in seq.int(nv_old)) {
      if(! is.na(vert_mapping[map_idx])) {
        rev_mapping[vert_mapping[map_idx]] = map_idx;
      }
    }
    if(any(rev_mapping < 0L)) {
      stop("Invalid mapping detected. That's a bug, please report.");
    }
    res = list('submesh'=new_mesh, 'vmap_full_to_submesh'=vert_mapping , 'vmap_submesh_to_full'=rev_mapping);
    return(res);
  }

  return(new_mesh);
}

