#!/usr/bin/which Rscript
#
# This script illustrates how to quickly smooth many vectors of per-vertex data
# on the same mesh.
#
#


library("haze");


mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage6_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = haze::mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood

nv = nrow(mesh$vertices); # number of mesh vertices

data1 = rnorm(nv, mean=1.0, sd=0.01); # Generate some random data.
data2 = rnorm(nv, mean=3.0, sd=0.01); # More!
data3 = rnorm(nv, mean=5.0, sd=0.01); # ...
data4 = rnorm(nv, mean=7.0, sd=0.01); # ...
data5 = rnorm(nv, mean=9.0, sd=0.01); # ...
data_matrix = rbind(data1, data2, data3, data4, data5); # your data as a matrix.

options("mc.cores" = 5L);   # Request 5 CPU cores.

# Get the smoothed matrix:
smoothed = haze::pervertexdata.smoothnn.adj(mesh_adj, data_matrix, num_iter = 15L);


