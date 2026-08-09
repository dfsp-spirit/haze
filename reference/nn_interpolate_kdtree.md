# Get per-vertex data at vertices closest to the given query coordinates on the mesh.

Return per-vertex data at the vertices closest to the given query
points.

## Usage

``` r
nn_interpolate_kdtree(query_coordinates, mesh, pervertex_data)
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

## Value

the per-vertex data for the vertices closest to the query coordinates.
