# haze
Fast smoothing of per-vertex data on triangular meshes for R.


## About

This package package performs smoothing of per-vertex data on triangular meshes, as illustrated here:

![Vis](./web/haze.jpg?raw=true "Per-vertex data on a brain mesh before (left) and after (right) smoothing.")

**Fig.1**: *Per-vertex data on a brain mesh before (left) and after (right) smoothing. White represents NA values.*

Such smoothing is typically used to reduce high-frequency noise and improve SNR.

### What haze is not

To avoid any confusion: haze does not smooth the mesh itself, use [Rvcg](https://github.com/zarquon42b/Rvcg)`::vcgSmooth()` for that.

## Features

* nearest neighbor smoothing based on edge distance (e.g., `k`-ring neighborhood of each vertex, with arbitrary `k`)
* Gaussian smoothing based on geodesic distances on the mesh (recommended for smaller meshes only)


### Properties

* fast: the smoothing is done in C++
* supports re-use of neighborhood data for faster smoothing of several datasets on the same mesh
* even faster for several overlays: several overlays can be smoothed in parallel on multi-core CPUs.
* works with various standard mesh file formats like PLY, OBJ, OFF as well as FreeSurfer brain meshes
* the internal mesh representation is `tmesh3d` from the `rgl` package, which allows for easy visualization
* ignores values set to `NA` during smoothing, which can be used to mask certain mesh areas (like the medial wall in neuroimaging)
* comes with pre-computed neighborhood data for meshes commonly used in surface-based neuroimaging (FreeSurfer's fsaverage and fsaverage6)


## Installation

```r
install.packages("devtools");
devtools::install_github("dfsp-spirit/Rvcg", ref="smooth_pervertex_data");
devtools::install_github("dfsp-spirit/haze");
```

## Documentation and Usage

* The most important function in the package is `pervertexdata.smoothnn()`. 
  - If you would like to pre-compute the mesh neighborhood once and re-use it for smoothing many datasets on the same mesh, use `mesh.adj` in combination with `pervertexdata.smoothnn.adj` instead.
* Help for a specific function can be accessed in the usual R manner: `?<function>`, where you replace `<function>` with a function name. Like this: `?pervertexdata.smoothnn`.
* Run `example(<function>)` to see a live demo that uses the function `<function>`. Like this: `example(pervertexdata.smoothnn)`.
* The [unit tests](./tests/testthat/) that come with this package are essentially a list of examples that illustrate how to use the functions.
* Here are some examples to get you started:

```r
library("haze");

# Example 1: Smooth neuroimaging data on a human brain mesh in FreeSurfer format (see Fig.1 above):
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
per_vertex_data = freesurferformats::read.fs.morph(system.file("extdata", "fsaverage_lh_thickness", package = "haze", mustWork = TRUE));
smoothed_data = pervertexdata.smoothnn(mesh, per_vertex_data, num_iter = 300L, k = 2);

# Example 2: Smooth random data on an rgl tetrahedron:
mesh2 = rgl::tetrahedron3d();
pvd = rnorm(nrow(mesh2$vb), mean = 5.0, sd = 1.0);
pvd_smoothed = pervertexdata.smoothnn(mesh2, pvd, num_iter = 30L);

# Example 3: Like 1, but re-use the adjacency list representation of the mesh to smooth several per-vertex data overlays on the same mesh:
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood
data1 = rnorm(length(mesh_adj), mean=1.0, sd=0.1); # generate random data
data2 = rnorm(length(mesh_adj), mean=5.0, sd=0.1); # generate more random data
smoothed_data1 = pervertexdata.smoothnn.adj(mesh_adj, data1, num_iter = 15L);
smoothed_data2 = pervertexdata.smoothnn.adj(mesh_adj, data2, num_iter = 15L);

# Example 4: Like 3, but smooth the overlays in parallel on several CPU cores:
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood
data1 = rnorm(length(mesh_adj), mean=1.0, sd=0.1); # generate random data
data2 = rnorm(length(mesh_adj), mean=5.0, sd=0.1); # generate more random data
pvd = rbind(data1, data2); # Turn your data into a matrix, one overlay per row.
options("mc.cores" = 2L);   # Request to run on 2 cores in parallel.
smoothed_pvd = pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L); # Compute the smoothed matrix.
```

## Credits

The fast mesh operations used in this package are implemented in the [Rvcg package](https://github.com/zarquon42b/Rvcg) by Stefan Schlager, which uses [VCGLIB](http://vcg.isti.cnr.it/vcglib/).


## Author and Getting help

The `haze` R package was written by [Tim Schäfer](http://rcmd.org/ts).

Please [open an issue](https://github.com/dfsp-spirit/haze/issues) here on GitHub if you have found a bug or need help.

