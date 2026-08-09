# Return pre-computed neighborhood data for specific meshes.

Return pre-computed neighborhood data for specific meshes.

## Usage

``` r
mesh.neigh.pre(meshname)
```

## Arguments

- meshname:

  a text identifier specifying the mesh you want connectivity data for.
  Currently supported meshes are listed here. 'lh_fsaverage': the left
  hemisphere of the FreeSurfer 6 fsaverage template. 'rh_fsaverage': the
  right hemisphere of the FreeSurfer 6 fsaverage template.
  'lh_fsaverage6': the left hemisphere of the FreeSurfer 6 fsaverage6
  template. 'rh_fsaverage6': the right hemisphere of the FreeSurfer 6
  fsaverage6 template.

## Value

list of vectors, the connectivity data as an adjacency list. The outer
list has length n, where n is the number of vertices in the graph. The
inner lists represent, for each vertex, all of its neighbors.
