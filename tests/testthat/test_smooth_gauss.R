

test_that("One can Gaussian smooth data on fsaverage3 using the raw mesh.", {
  fsmesh_file = system.file("extdata", "fsaverage3_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  num_verts = nrow(mesh$vertices);
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothgauss(mesh, data, 15L);

  testthat::expect_equal(length(smoothed_data), num_verts);
})
