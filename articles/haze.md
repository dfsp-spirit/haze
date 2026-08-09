# Getting Started with haze

## Introduction

`haze` provides fast smoothing of per-vertex data on triangular meshes,
along with utility functions for mesh manipulation and data
interpolation. The smoothing is implemented in C++ and supports parallel
processing for multiple overlays.

## Smoothing Per-Vertex Data

The main workhorse is
[`pervertexdata.smoothnn()`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md),
which performs nearest-neighbor smoothing using the $`k`$-ring
neighborhood of each vertex.

``` r

library("haze")

# Load a FreeSurfer mesh and per-vertex data
mesh <- freesurferformats::read.fs.surface(
  system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE))
thickness <- freesurferformats::read.fs.morph(
  system.file("extdata", "fsaverage_lh_thickness", package = "haze", mustWork = TRUE))

# Smooth: k=2 neighborhood, 300 iterations
smoothed <- pervertexdata.smoothnn(mesh, thickness, num_iter = 300L, k = 2L)
#> Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
#> Warning: 'rgl.init' failed, will use the null device.
#> See '?rgl.useNULL' for ways to avoid this warning.
```

## Re-Using the Mesh Neighborhood

When smoothing many datasets on the same mesh, pre-compute the adjacency
list once with
[`mesh.adj()`](https://dfsp-spirit.github.io/haze/reference/mesh.adj.md)
and then use
[`pervertexdata.smoothnn.adj()`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.adj.md):

``` r

mesh_adj <- mesh.adj(mesh, k = 1L)
```

This adjacency can be reused:

``` r

data1 <- rnorm(length(mesh_adj), mean = 1.0, sd = 0.1)
data2 <- rnorm(length(mesh_adj), mean = 5.0, sd = 0.1)

smoothed1 <- pervertexdata.smoothnn.adj(mesh_adj, data1, num_iter = 15L)
smoothed2 <- pervertexdata.smoothnn.adj(mesh_adj, data2, num_iter = 15L)
```

Multiple overlays (as a matrix, one per row) are automatically smoothed
in parallel:

``` r

pvd <- rbind(data1, data2)
options("mc.cores" = 2L)
smoothed_pvd <- pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L)
#> Handling 2 vertex overlays with 163842 values each in parallel on 2 CPU cores. (Use 'options("mc.cores"=N)' to request N cores).
```

## Creating Sub-Meshes

Use
[`submesh.vertex()`](https://dfsp-spirit.github.io/haze/reference/submesh.vertex.md)
to extract a patch of a mesh defined by vertex indices:

``` r

# Extract a sub-mesh of the first 1000 vertices
small_mesh <- submesh.vertex(mesh, old_vertex_indices_to_use = 1:1000)
```

## Working with k-d Trees

[`find_nv_kdtree()`](https://dfsp-spirit.github.io/haze/reference/find_nv_kdtree.md)
finds the nearest mesh vertex for given query coordinates:

``` r

query_points <- matrix(c(0, 0, 0, 10, 10, 10), ncol = 3, byrow = TRUE)
nearest_vertex_idx <- find_nv_kdtree(query_points, mesh)
```

[`nn_interpolate_kdtree()`](https://dfsp-spirit.github.io/haze/reference/nn_interpolate_kdtree.md)
assigns per-vertex data to query points via nearest-neighbor lookup,
while
[`linear_interpolate_kdtree()`](https://dfsp-spirit.github.io/haze/reference/linear_interpolate_kdtree.md)
performs linear interpolation — useful for mapping morphometry data
between subjects with aligned spherical meshes.

## Mesh Conversion Utilities

Convert between FreeSurfer surfaces and `rgl` `tmesh3d` objects using
the internal utilities `haze:::fs.surface.to.tmesh3d()` and
`haze:::tmesh3d.to.fs.surface()`:

``` r

tmesh <- haze:::fs.surface.to.tmesh3d(mesh)
fs_surface <- haze:::tmesh3d.to.fs.surface(tmesh)
```

## Further Reading

- Full API reference:
  [dfsp-spirit.github.io/haze/reference/](https://dfsp-spirit.github.io/haze/reference/)
- Run `example(<function>)` in R for interactive demos of any function.
