# Smooth data, C++ version.

Smooth data, C++ version.

## Usage

``` r
pervertexdata.smoothnn.adj.cpp(mesh_adj, pvdata, num_iter)
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

## Value

numerical vector, the smoothed data.
