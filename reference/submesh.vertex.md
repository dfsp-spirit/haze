# Create a submesh including only the given vertices.

Create a submesh including only the given vertices.

## Usage

``` r
submesh.vertex(surface_mesh, old_vertex_indices_to_use, ret_mappings = FALSE)
```

## Arguments

- surface_mesh:

  an `fs.surface` instance, the original mesh. See
  [`read.fs.surface`](https://rdrr.io/pkg/freesurferformats/man/read.fs.surface.html)
  or
  [`subject.surface`](https://rdrr.io/pkg/fsbrain/man/subject.surface.html)
  to get one. Can also be an rgl tmesh, see
  [`tmesh3d`](https://dmurdoch.github.io/rgl/dev/reference/mesh3d.html).

- old_vertex_indices_to_use:

  integer vector, the vertex indices of the 'surface_mesh' that should
  be used to construct the new sub mesh.

- ret_mappings:

  whether to return the vertex mappings. If `TRUE`, the return value
  becomes a list with entries 'submesh', 'vmap_full_to_submesh', and
  'vmap_submesh_to_full'.

## Value

the new mesh, made up of the given 'old_vertex_indices_to_use' and all
(complete) faces that exist between the query vertices in the source
mesh. But see 'ret_mapping' parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
if(requireNamespace("fsbrain, quietly=T")) {
sjd = fsbrain::fsaverage.path(T);
sj = "fsaverage";
mesh = fsbrain::subject.surface(sjd, sj, hemi="lh");
lab = fsbrain::subject.label(sjd, sj, "cortex", hemi = "lh");
sm = submesh.vertex(mesh, lab);
fsbrain::vis.fs.surface(mesh); # show the full mesh.
fsbrain::vis.fs.surface(sm);   # show only the cortex.
}} # }
```
