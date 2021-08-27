
test_that("One can smooth data on fsaverage6 using pre-computed adjacency.", {
  vvfile = system.file("extdata", "fsaverage6_lh_white_meshdist_edge_1.vv", package = "haze", mustWork = TRUE);
  vv = read.vv(vvfile);

  num_verts = 40962L;
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothnn.adj(vv, data, 15L);

  testthat::expect_equal(length(smoothed_data), num_verts);
})


test_that("One can smooth data on fsaverage using the raw mesh.", {
  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  num_verts = nrow(mesh$vertices);
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothnn(mesh, data, 15L);

  ## This shows that the C++ version is a lot faster than the R version (55 times faster on our Linux server).
  #microbenchmark::microbenchmark(pervertexdata.smoothnn(mesh, data, 100L, method="C++"), pervertexdata.smoothnn(mesh, data, 100L, method="R"), times=10L);

  testthat::expect_equal(length(smoothed_data), num_verts);
})


test_that("One can compute mesh adjacency and re-use it for smoothing.", {
  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);

  mesh_adj = mesh.adj(fsmesh_file, k = 1L);

  num_verts = length(mesh_adj);
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothnn.adj(mesh_adj, data, 15L);

  testthat::expect_equal(length(smoothed_data), num_verts);

  ## Compare performance with pre-computed mesh vs having to compute mesh (only makes sense if you smooth several data sets).
  ## Obviously, the pervertexdata.smoothnn.adj version should be faster, but the difference is most noticeable for small iteration counts:
  #mesh = freesurferformats::read.fs.surface(fsmesh_file);
  #microbenchmark::microbenchmark(pervertexdata.smoothnn.adj(mesh_adj, data, 100L), pervertexdata.smoothnn(mesh, data, 100L), times=5L);
  #microbenchmark::microbenchmark(pervertexdata.smoothnn.adj(mesh_adj, data, 10L), pervertexdata.smoothnn(mesh, data, 10L), times=5L);
})


test_that("Smoothing of thickness data looks plausible.", {
  # This test produces 2 figures, look at them.
  testthat::skip_on_cran();

  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);
  pervertexdata_file = system.file("extdata", "fsaverage_lh_thickness", package = "haze", mustWork = TRUE);

  thickness = freesurferformats::read.fs.morph(pervertexdata_file);
  smooth_thickness = pervertexdata.smoothnn(fsmesh_file, thickness, 300L);

  if(requireNamespace("fsbrain", quietly = TRUE)) {

    #cm1 = fsbrain::vis.data.on.fsaverage(morph_data_lh = thickness, morph_data_rh = NA);
    #cm2 = fsbrain::vis.data.on.fsaverage(morph_data_lh = smooth_thickness, morph_data_rh = NA);

    #fsbrain::vis.export.from.coloredmeshes(cm1, output_img = "~/haze_thickness_before.png");
    #fsbrain::vis.export.from.coloredmeshes(cm2, output_img = "~/haze_thickness_after.png");

  } else {
    testthat::skip("This test requires the optional 'fsbrain' package to be installed.");
  }
  testthat::expect_equal(1L, 1L); # Tests without checks would be skipped by testthat.
})


test_that("Ignoring NA values in the data works as expected.", {
  testthat::skip_on_cran();
  # This test produces 2 figures, look at them.

  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);
  pervertexdata_file = system.file("extdata", "fsaverage_lh_thickness", package = "haze", mustWork = TRUE);
  label_file = system.file("extdata", "fsaverage_lh_cortex_label", package = "haze", mustWork = TRUE);

  thickness = freesurferformats::read.fs.morph(pervertexdata_file);

  # Set medial wall thickness values to NA using the cortex label file.
  num_verts = length(thickness);
  label = freesurferformats::read.fs.label(label_file);
  mask = rep(TRUE, num_verts);
  mask[label] = FALSE;
  thickness[mask] = NA;
  testthat::expect_equal(length(which(is.na(thickness))), 13887L);


  smooth_thickness = pervertexdata.smoothnn(fsmesh_file, thickness, 300L);
  testthat::expect_equal(length(which(is.na(smooth_thickness))), 13887L);

  if(requireNamespace("fsbrain", quietly = TRUE)) {
    #cm1 = fsbrain::vis.data.on.fsaverage(morph_data_lh = thickness, morph_data_rh = NA);
    #cm2 = fsbrain::vis.data.on.fsaverage(morph_data_lh = smooth_thickness, morph_data_rh = NA);

    #fsbrain::vis.export.from.coloredmeshes(cm1, output_img = "~/haze_thickness_masked_before.png");
    #fsbrain::vis.export.from.coloredmeshes(cm2, output_img = "~/haze_thickness_masked_after.png");
  }

})


