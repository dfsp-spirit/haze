# Interpolate values within mesh triangles in C++.

Interpolate values within mesh triangles in C++.

## Usage

``` r
interp_tris_cpp_wrapper(
  query_coordinates,
  mesh_vertices,
  nearest_face_vertices,
  pervertex_data,
  iwd_beta = 2
)
```

## Arguments

- query_coordinates:

  nx3 numerical matrix of x,y,z coordinates. These are typically the
  vertex positions of a second (spherical!) mesh for that you need
  per-vertex data (e.g., the `fsaverage6` mesh).

- mesh_vertices:

  `nx3` matrix of x,y,z Cartesian coordinates for the n mesh vertices
  (of the mesh that has the 'pervertex_data').

- pervertex_data:

  numerical vector, the continuous per-vertex data for the vertices of
  the mesh.

- iwd_beta:

  scalar double, the `beta` parameter for the inverse distance weight
  interpolation with the triangle. See details.

## Value

vector with one entry per query coordinate (i.e.,
`nrow(query_coordinates)`).

## Note

No validation of input data is done, do this yourself before passing
values to this function to avoid hard crashes (which happen when junk is
passed to `C++`).
