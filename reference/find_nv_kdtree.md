# Find nearest mesh vertex for query coordinates using kdtree.

Find nearest mesh vertex for query coordinates using kdtree.

## Usage

``` r
find_nv_kdtree(query_coordinates, mesh, threads = parallel::detectCores())
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

- threads:

  integer, the number of threads to run in parallel.

## Value

named list with keys 'index' and 'distance'. 'index': integer vector,
the `n` vertex indices which are closest to the `nx3` matrix of
query_coordinates. 1-based indices for R are returned. 'distance':
double vector, the distances to the respective vertices in the 'index'
key.

## Note

@note The mesh must be spherical, and the query_coordinates must be
located on the mesh sphere.

## See also

`https://github.com/ThomasYeoLab/CBIG/blob/master/external_packages/SD/SDv1.5.1-svn593/BasicTools/MARS_findNV_kdTree.m`
