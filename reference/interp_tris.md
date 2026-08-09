# Interpolate values within mesh triangles in R.

Interpolate values within mesh triangles in R.

## Usage

``` r
interp_tris(
  query_coordinates,
  mesh_vertices,
  nearest_face_vertices,
  pervertex_data,
  iwd_beta = 2,
  cpp = TRUE
)
```

## Arguments

- query_coordinates:

  nx3 numerical matrix of x,y,z coordinates. These are typically the
  vertex positions of a second (spherical!) mesh for that you need
  per-vertex data (e.g., the `fsaverage6` mesh).

- mesh_vertices:

  `nx3` matrix of x,y,z cartesian coordinates for the n mesh vertices
  (of the mesh that has the 'pervertex_data').

- pervertex_data:

  numerical vector, the continuous per-vertex data for the vertices of
  the mesh.

- iwd_beta:

  scalar double, the `beta` parameter for the inverse distance weight
  interpolation with the triangle. See details.

- cpp:

  logical, whether to use the much faster `C++` version. Leave alone
  unless you know what you are doing. Changing this to `FALSE` will make
  this a lot slower.

## Value

vector with one entry per query coordinate (i.e.,
`nrow(query_coordinates)`).
