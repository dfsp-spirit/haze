# smoothr
Fast smoothing of per-vertex data on triangular meshes for R.


## About

This package performs smoothing of per-vertex data on triangular meshes, as illustrated here:

![Vis](./web/smoothr.jpg?raw=true "Per-vertex data on a brain mesh before (left) and after (right) smoothing.")

**Fig.1**: *Per-vertex data on a brain mesh before (left) and after (right) smoothing. White represents NA values.*

Such smoothing is typically used to reduce high-frequency noise and improve SNR. To avoid any confusion: `smoothr` does **not** smooth the mesh itself. One can use `Rvcg::vcgSmooth` from the `Rvcg` R package to do that.


## Features

* nearest neighbor smoothing based on edge distance (e.g., `k`-ring neighborhood of each vertex, with arbitrary `k`)
* Gaussian smoothing based on geodesic distances on the mesh (geodesic computation is slow for large meshes), WIP

### Properties

* fast: the smoothing is done in C++
* works with various standard mesh file formats like PLY, OBJ, OFF as well as FreeSurfer brain meshes
* the internal mesh representation is `tmesh3d` from the `rgl` package, which allows for easy visualization
* supports re-use of neighborhood data for faster smoothing of several datasets on the same mesh
* ignores values set to `NA` during smoothing, which can be used to mask certain mesh areas (like the medial wall in neuroimaging)
* comes with pre-computed neighborhood data for meshes commonly used in surface-based neuroimaging (FreeSurfer's fsaverage and fsaverage6)


## Usage

```
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "smoothr", mustWork = TRUE));
per_vertex_data = freesurferformats::read.fs.morph(system.file("extdata", "fsaverage_lh_thickness", package = "smoothr", mustWork = TRUE));

smoothed_data = pervertexdata.smoothnn(mesh, per_vertex_data, num_iter = 300L);
```

This is WIP, come back another day.

