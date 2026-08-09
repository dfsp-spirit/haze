# Smooth matrix in parallel

Smooth matrix in parallel

## Usage

``` r
pervertexdata.smoothnn.adj.mat(
  mesh_adj,
  pvdata,
  num_iter,
  method = "C++",
  silent = getOption("haze.silent", default = FALSE)
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

  numerical matrix or data.frame of per-vertex-data for the mesh. Each
  row must be an overlay with one value per vertex. Rows are independent
  and smoothed in parallel. Data values of `NA` will be ignored,
  allowing you to mask parts of the data.

- num_iter:

  positive integer, number of smoothing iterations.

- method:

  character string, one of 'C++' or 'R'. The C++ version is much faster
  (about 30 times faster on our test machine), and there is little
  reason to ever use the R version in production code, so just ignore
  this.

- silent:

  logical, whether to suppress output messages.

## Note

There is no need to call this manually, it will be called by
`pervertexdata.smooth.adj` when its input is a matrix.

To set the number of cores to use, set the 'mc_cores' option like this:
`options("mc.cores"=22L);` before calling this function.

## Examples

``` r
if (FALSE) { # \dontrun{
mesh = rgl::tetrahedron3d();
mesh_adj = mesh.adj(mesh, k = 1L);
pvd1 = rnorm(nrow(mesh$vb), mean = 1.0, sd = 1.0);
pvd2 = rnorm(nrow(mesh$vb), mean = 4.0, sd = 1.0);
pvd3 = rnorm(nrow(mesh$vb), mean = 7.0, sd = 1.0);
pvd_matrix = rbind(pvd1, pvd2, pvd3);
# or: pvd_matrix = matrix(data=rnorm(nrow(mesh$vb)*5, mean=5.0, sd=1.0), nrow=5);
pvd_smoothed = haze:::pervertexdata.smoothnn.adj.mat(mesh_adj, pvd_matrix, num_iter = 30L);
} # }
```
