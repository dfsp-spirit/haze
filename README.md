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
smoothed_pvd = pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L); # Compute the smoothed matrix. When a matrix is passed, the rows are automatically handled in parallel, there is nothing more to do.
```

### Benchmarks

This requires the `microbenchmark` package.

#### Measuring C++ versus R performance

The smoothing is done in C++ by default, but it is possible to manually select a pure R version instead. The following benchmark compares the performance between the C++ and the R versions of the same algorithm.

```R
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood

pvd = rnorm(length(mesh_adj), mean=1.0, sd=0.1);

microbenchmark::microbenchmark(pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L, method="R"), pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L, method="C++"), times=5L);
```

On my machine, this shows that the C++ version is about 29 times faster for the given data.


#### Measuring multi-core performance

To get a rough idea on how much faster the C++ smoothing runs with several cores, benchmark it for your data. E.g.:

```R
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood

num_overlays = 25;
pvd = matrix(data=rnorm(length(mesh_adj)*num_overlays, mean=1.0, sd=0.1), nrow=num_overlays); # generate random data

# Compare with 1 versus 5 cores:

microbenchmark::microbenchmark({options("mc.cores" = 1L); pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L);}, {options("mc.cores" = 5L); pervertexdata.smoothnn.adj(mesh_adj, pvd, num_iter = 15L);}, times=5L);
```

On my machine, this shows that code is about 4 times faster when using 5 CPU cores instead of a single one.

So with 5 cores and the C++ version (Example 4 above), the smoothing runs about 120 times faster compared to the single-core, pure R version for the data above on my machine.

##### Number of cores versus execution time

One could also measure this using from 1 up to 10 cores and plot the resulting execution time.

```R
# Continued from the 'Measuring multi-core performance' example.
# The next 2 lines just save us some typing in the microbenchmark line.
rf <- function(nc, details) { options("mc.cores" = nc); pervertexdata.smoothnn.adj(details$mesh_adj, details$pvd, details$num_iter); }
d = list("mesh_adj"=mesh_adj, "pvd"=pvd, "num_iter"=15L);
mb = microbenchmark::microbenchmark(rf(1,d), rf(2,d), rf(3,d), rf(4,d), rf(5,d), rf(6,d), rf(7,d), rf(8,d), rf(9,d), rf(10,d), times=3L);
plot(mb, xaxt='n', xlab="Number of cores", ylab = "Execution time [ns]");
axis(1, at = seq(10));
```

![Vis](./web/haze_multicore.png?raw=true "Haze multi-core performance.")

For this data set and my machine, using more than 5 cores does not seem to help much.

## Credits

The fast mesh operations used in this package are implemented in the [Rvcg package](https://github.com/zarquon42b/Rvcg) by Stefan Schlager, which uses [VCGLIB](http://vcg.isti.cnr.it/vcglib/).


## Author and Getting help

The `haze` R package was written by [Tim Schäfer](http://rcmd.org/ts).

Please [open an issue](https://github.com/dfsp-spirit/haze/issues) here on GitHub if you have found a bug or need help.

