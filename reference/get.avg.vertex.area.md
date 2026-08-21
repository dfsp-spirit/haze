# Compute the mean vertex area of a triangular mesh.

Compute the mean vertex area of a triangular mesh.

## Usage

``` r
get.avg.vertex.area(surface = NULL, avg_vertex_area = NULL)
```

## Arguments

- surface:

  a mesh, represented as an `fs.surface` instance from the
  `freesurferformats` package, a `tmesh3d` instance from `rgl`, or a
  character string representing the path of a surface file to load.

- avg_vertex_area:

  numeric scalar, the mean vertex area in squared mesh units (e.g.,
  mm^2). If `NULL` (the default), it is computed from `surface`.

## Value

numeric scalar, the mean (average) vertex area of the mesh in squared
mesh units.
