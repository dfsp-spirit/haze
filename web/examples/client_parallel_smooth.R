#!/usr/bin/which Rscript
#
# This script illustrates how to quickly smooth many vectors of per-vertex data
# on the same mesh.
#
#


library("haze");
library("foreach");
library("parallel");
library("doParallel");

mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage6_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = haze::mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood

nv = nrow(mesh$vertices); # number of mesh vertices

data1 = rnorm(nv, mean=1.0, sd=0.01); # Generate some random data.
data2 = rnorm(nv, mean=3.0, sd=0.01); # More!
data3 = rnorm(nv, mean=5.0, sd=0.01); # ...
data4 = rnorm(nv, mean=7.0, sd=0.01); # ...
data5 = rnorm(nv, mean=9.0, sd=0.01); # ...

data_matrix = rbind(data1, data2, data3, data4, data5); # your data as a matrix.
num_rows = nrow(data_matrix);

num_cores_to_use = parallel::detectCores() - 1L;
cat(sprintf("Smoothing %d overlays for %d vertices using %d cores in parallel.\n", nrow(data_matrix), ncol(data_matrix), num_cores_to_use));
cluster = parallel::makeCluster(num_cores_to_use);
doParallel::registerDoParallel(cluster);

smoothed_data_matrix <- foreach::foreach(mr=iter(data_matrix, by='row'), .combine=rbind, .packages="haze") %do% {
  smoothed_row = haze::pervertexdata.smoothnn.adj(mesh_adj, mr, num_iter = 15L);
  smoothed_row
}

cat(sprintf("Smoothed data matrix with dimensions (%d subjects, %d mesh vertices).\n", dim(smoothed_data_matrix)[1], dim(smoothed_data_matrix)[2]));

parallel::stopCluster(cluster);

