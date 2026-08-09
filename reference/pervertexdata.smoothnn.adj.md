# Smooth per-vertex data using nearest-neighbor smoothing based on mesh adjacency information.

Smooth per-vertex data using nearest-neighbor smoothing based on mesh
adjacency information.

## Usage

``` r
pervertexdata.smoothnn.adj(
  mesh_adj,
  pvdata,
  num_iter,
  method = "C++",
  silent = getOption("haze.silent", default = TRUE)
)
```

## Arguments

- mesh_adj:

  list of vectors of integers, the adjacency list representation of the
  mesh. One can use the pre-computed adjacency for some special meshes,
  see
  [`mesh.neigh.pre`](https://dfsp-spirit.github.io/haze/reference/mesh.neigh.pre.md).
  Data for vertices should include the vertex itself.

- pvdata:

  numerical vector of per-vertex-data for the mesh, one value per
  vertex. Data values of `NA` will be ignored, allowing you to mask
  parts of the data. If you pass an `n x m` matrix or data.frame, the
  `n` rows will be treated as (independent) overlays that should be
  smoothed in parallel. To set the number of cores to use for parallel
  processing, set the 'mc_cores' option like this:
  `options("mc.cores"=22L);` before calling this function. Data.frames
  and matrices with a single row will be converted to vectors, and the
  return value will be a vector in that case.

- num_iter:

  positive integer, number of smoothing iterations.

- method:

  character string, one of 'C++' or 'R'. The C++ version is much faster
  (about 30 times faster on our test machine), and there is little
  reason to ever use the R version in production code, so just ignore
  this.

- silent:

  logical, whether to suppress output messages.

## Value

numerical vector, the smoothed data (for vector input). If `pvdata` is a
matrix or a data.frame (with more than a single column), the result is
also a matrix or data.frame.

## See also

[`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md)
if you have a mesh and still need the connectivity to be computed.

## Examples

``` r
if (FALSE) { # \dontrun{
mesh = rgl::tetrahedron3d();
mesh_adj = mesh.adj(mesh, k = 1L);
pvd = rnorm(nrow(mesh$vb), mean = 5.0, sd = 1.0);
pvd_smoothed = pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 30L);
} # }
```
