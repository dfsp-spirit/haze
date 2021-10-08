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

mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
mesh_adj = haze::mesh.adj(mesh, k = 1L); # compute 1-ring neighborhood

nv = nrow(mesh$vertices); # number of mesh vertices

data1 = rnorm(nv, 5.0, 1.0); # Generate some random data.
data2 = rnorm(nv, 5.5, 2.0); # More!
data3 = rnorm(nv, 6.0, 2.0); # ...
data4 = rnorm(nv, 6.5, 2.0); # ...
data5 = rnorm(nv, 7.0, 2.0); # ...

data_matrix = rbind(data1, data2, data3, data4, data5); # your data as a matrix.
num_rows = nrow(data_matrix);

num_cores_to_use = parallel::detectCores() - 1L;
cat(sprintf("Smoothing %d overlays using %d cores in parallel.\n", nrow(data_matrix), num_cores_to_use));
cluster = parallel::makeCluster(num_cores_to_use);
registerDoParallel(cluster);

smoothed_data_matrix = foreach::foreach(vec_idx=1:num_rows, .combine=cbind) %dopar% {
  library("haze");
  smoothed_row = haze::pervertexdata.smoothnn.adj(mesh_adj, data_matrix[vec_idx,], num_iter = 15L, k=1L);
  smoothed_row #Equivalent to smoothed_data_matrix = cbind(smoothed_data_matrix, smoothed_row)
}

cat(sprintf("Smoothed data matrix with dimensions (%d subjects, %d mesh vertices).\n", dim(smoothed_data_matrix)[1], dim(smoothed_data_matrix)[2]));

parallel::stopCluster(cluster);
