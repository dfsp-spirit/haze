# Compute vertex neighborhoods for a mesh.

Compute vertex neighborhoods for a mesh.

## Usage

``` r
mesh.adj(surface, k = 1L)
```

## Arguments

- surface:

  a mesh, represented as an `fs.surface` instance from the
  `freesurferformats` package or a `tmesh3d` instance from `rgl`, or a
  character string representing the path of a mesh to load with
  [`freesurferformats::read.fs.surface`](https://rdrr.io/pkg/freesurferformats/man/read.fs.surface.html).

- k:

  scalar positive integer, the k value for the k-ring neighborhood. For
  k=1, this function computes the adjacency list representation of the
  graph (where the neighbors include the vertex itself).

## Value

list of integer vectors, the neighborhood data

## Examples

``` r
if (FALSE) { # \dontrun{
mesh = rgl::tetrahedron3d();
mesh_adj = mesh.adj(mesh, k = 1L);
} # }
```
