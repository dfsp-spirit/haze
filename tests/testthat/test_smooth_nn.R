
test_that("One can smooth data on fsaverage6 using pre-computed adjacency.", {
  vvfile = system.file("extdata", "fsaverage6_lh_white_meshdist_edge_1.vv", package = "smoothr", mustWork = TRUE);
  vv = read.vv(vvfile);

  num_verts = 40962L;
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothnn.adj(vv, data, 15L);

  testthat::expect_equal(length(smoothed_data), num_verts);
})


test_that("One can smooth data on fsaverage using the raw mesh.", {
  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "smoothr", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  num_verts = nrow(mesh$vertices);
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothnn(mesh, data, 15L);

  # This shows that the C++ version is a lot faster than the R version (55 times faster on our Linux server).
  #microbenchmark::microbenchmark(pervertexdata.smoothnn(mesh, data, 100L, method="C++"), pervertexdata.smoothnn(mesh, data, 100L, method="R"), times=10L)

  testthat::expect_equal(length(smoothed_data), num_verts);
})


