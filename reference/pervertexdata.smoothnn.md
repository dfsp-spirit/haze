# Smooth per-vertex data based on mesh.

Smooth per-vertex data based on mesh.

## Usage

``` r
pervertexdata.smoothnn(surface, pvdata, num_iter, k = 1L, method = "C++")
```

## Arguments

- surface:

  a mesh, represented as an `fs.surface` instance from the
  `freesurferformats` package or a `tmesh3d` instance from `rgl`, or a
  character string representing the path of a mesh to load with
  [`freesurferformats::read.fs.surface`](https://rdrr.io/pkg/freesurferformats/man/read.fs.surface.html).

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

- k:

  scalar positive integer, the k value for the k-ring neighborhood. For
  k=1, this function computes the adjacency list representation of the
  graph (where the neighbors include the vertex itself).

- method:

  character string, one of 'C++' or 'R'. The C++ version is much faster
  (about 30 times faster on our test machine), and there is little
  reason to ever use the R version in production code, so just ignore
  this.

## Value

numerical vector, the smoothed data.

## See also

[`pervertexdata.smoothnn.adj`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.adj.md)
if you already have pre-computed adjacency data for the mesh. Using that
data can increase performance considerably, especially if you need to
smooth many data sets.

## Examples

``` r
if (FALSE) { # \dontrun{
mesh = rgl::tetrahedron3d();
pvd = rnorm(nrow(mes2$vb), mean = 5.0, sd = 1.0);
pvd_smoothed = pervertexdata.smoothnn(mesh, pvd, num_iter = 30L);
} # }
```
