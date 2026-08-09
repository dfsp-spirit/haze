# Interpolate per-vertex data at the query points. Or map per-vertex data between subjects.

This method uses inverse distance weight interpolation within a
triangle. First, the face of the `mesh` that the `query_coordinate`
falls into is determined. Then results in 3 vertices with respective
per-vertex data, and a query coordinate. We then compute the distance to
all 3 vertices, and perform inverse distance weight interpolation with a
beta setting defined by parameter `iwd_beta`.

## Usage

``` r
linear_interpolate_kdtree(
  query_coordinates,
  mesh,
  pervertex_data,
  iwd_beta = 2,
  ...
)
```

## Arguments

- query_coordinates:

  nx3 numerical matrix of x,y,z coordinates. These are typically the
  vertex positions of a second (spherical!) mesh for that you need
  per-vertex data (e.g., the `fsaverage6` mesh).

- mesh:

  fs.surface instance, see
  [`read.fs.surface`](https://rdrr.io/pkg/freesurferformats/man/read.fs.surface.html)
  or
  [`subject.surface`](https://rdrr.io/pkg/fsbrain/man/subject.surface.html)
  to get one, or turn an `rgl` `tmesh` into one with
  [`tmesh3d.to.fs.surface`](https://rdrr.io/pkg/fsbrain/man/tmesh3d.to.fs.surface.html).

- pervertex_data:

  numerical vector, the continuous per-vertex data for the vertices of
  the mesh.

- iwd_beta:

  scalar double, the `beta` parameter for the inverse distance weight
  interpolation with the triangle. See details.

- ...:

  ignore, passed on to internal function `interp_tris`.

## Value

named list with entries: 'interp_values', the numerical vector of
interpolated data at the query_coordinates. 'nearest_vertex_in_face' the
nearest vertex in the face that the respective query coordinate falls
into, 'nearest_face' the index of the nearest face that the respective
query coordinate falls into.

## Note

The mesh must be spherical, and the `query_coordinates` must also be
located on the mesh sphere.
