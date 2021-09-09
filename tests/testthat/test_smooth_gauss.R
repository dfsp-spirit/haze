

test_that("One can Gaussian smooth data on fsaverage3 using the raw mesh.", {
  fsmesh_file = system.file("extdata", "fsaverage3_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  num_verts = nrow(mesh$vertices);
  data = rnorm(num_verts, 5.0, 1.0);

  smoothed_data = pervertexdata.smoothgauss(mesh, data, 15L);

  testthat::expect_equal(length(smoothed_data), num_verts);
})


test_that("One can Gaussian smooth data that includes NA values.", {
  fsmesh_file = system.file("extdata", "fsaverage3_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  num_verts = nrow(mesh$vertices);
  data = rnorm(num_verts, 5.0, 1.0);
  data[100:199] = NA;

  smoothed_data = pervertexdata.smoothgauss(mesh, data, 15L);

  testthat::expect_equal(length(smoothed_data), num_verts);

  testthat::expect_equal(length(which(is.na(smoothed_data))), 100);
  testthat::expect_true(all(is.na(smoothed_data[100:199])));
  testthat::expect_true(! (any(is.na(smoothed_data[0:99]))));
  testthat::expect_true(! (any(is.na(smoothed_data[200:length(smoothed_data)]))));
})


test_that("Gaussian smoothing of thickness data looks plausible with NA values.", {
  # This test produces 2 figures, look at them.
  testthat::skip_on_cran();


  if(requireNamespace("fsbrain", quietly = TRUE)) {

    fsbrain::download_fsaverage3(TRUE);
    sjd = fsbrain::fsaverage.path(TRUE);
    sj = "fsaverage3";
    mesh = fsbrain::subject.surface(sjd, sj, "white", hemi="lh");
    pvd = fsbrain::subject.morph.native(sjd, sj, "thickness", hemi="lh", cortex_only = TRUE); # mask non-cortex to get NAs
    smoothed_pvd = pervertexdata.smoothgauss(mesh, pvd, fwhm = 5.0);

    #cm1 =  fsbrain::vis.fs.surface(mesh, per_vertex_data = pvd);
    #cm2 =  fsbrain::vis.fs.surface(mesh, per_vertex_data = smoothed_pvd);

    #fsbrain::vis.export.from.coloredmeshes(cm1, output_img = "~/haze_thickness_before_gauss.png");
    #fsbrain::vis.export.from.coloredmeshes(cm2, output_img = "~/haze_thickness_after_gauss.png");

  } else {
    testthat::skip("This test requires the optional 'fsbrain' package to be installed.");
  }
  testthat::expect_equal(1L, 1L); # Tests without checks would be skipped by testthat.
})

