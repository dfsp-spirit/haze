# haze
Fast smoothing of per-vertex data on triangular meshes for R.


## About

This package performs smoothing of per-vertex data on triangular meshes, as illustrated here:

![Vis](./web/smoothr.jpg?raw=true "Per-vertex data on a brain mesh before (left) and after (right) smoothing.")

**Fig.1**: *Per-vertex data on a brain mesh before (left) and after (right) smoothing. White represents NA values.*

Such smoothing is typically used to reduce high-frequency noise and improve SNR. To avoid any confusion: `haze` does **not** smooth the mesh itself, use `Rvcg::vcgSmooth()` for that.


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

## Installation

This is WIP, come back another day.

## Usage

* The most important function in the package is `pervertexdata.smoothnn()`. 
* If you would like to pre-compute the mesh neighborhood once and re-use it for smoothing many datasets on the same mesh, use `mesh.adj` in combination with `pervertexdata.smoothnn.adj` instead.
* Help for a specific function can be accessed in the usual R manner: `?<function>`, where you replace `<function>` with a function name. Like this: `?pervertexdata.smoothnn`.
* The [unit tests](./tests/testthat/) that come with this package are essentially a list of examples that illustrate how to use the functions.
* Here are some examples to get you started:

```r
library("haze");

# Example 1: Smooth neuroimaging data on a human brain mesh in FreeSurfer format (see Fig.1 above):
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "smoothr", mustWork = TRUE));
per_vertex_data = freesurferformats::read.fs.morph(system.file("extdata", "fsaverage_lh_thickness", package = "smoothr", mustWork = TRUE));
smoothed_data = pervertexdata.smoothnn(mesh, per_vertex_data, num_iter = 300L, k = 2);

Example 2: Smooth random data on an rgl tetrahedron:
mesh2 = rgl::tetrahedron3d();
pvd = rnorm(nrow(mesh2$vb), mean = 5.0, sd = 1.0);
pvd_smoothed = pervertexdata.smoothnn(mesh2, pvd, num_iter = 30L);
```
