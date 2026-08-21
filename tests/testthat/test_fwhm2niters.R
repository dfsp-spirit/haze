# Tests for the conversion between FreeSurfer Gaussian smoothing (FWHM) and
# nearest-neighbor smoothing iteration counts (fwhm2niters / niters2fwhm).

testthat::test_that("fwhm2niters reproduces the FreeSurfer v6 fsaverage FWHM-to-iterations mapping.", {
  # These are the values FreeSurfer reports for the fsaverage templates when
  # smoothing per-vertex data with k=1 (1-ring) NN smoothing. They are based on
  # the FreeSurfer group-average mean vertex area (see README, section on
  # Gaussian vs NN smoothing). The value ~0.5015 mm^2 is what FreeSurfer
  # effectively uses for the fsaverage templates.
  freeSurfer_avg_vertex_area = 0.5015;

  testthat::expect_equal(fwhm2niters(2,  avg_vertex_area = freeSurfer_avg_vertex_area), 3L);
  testthat::expect_equal(fwhm2niters(5,  avg_vertex_area = freeSurfer_avg_vertex_area), 18L);
  testthat::expect_equal(fwhm2niters(10, avg_vertex_area = freeSurfer_avg_vertex_area), 74L);
  testthat::expect_equal(fwhm2niters(15, avg_vertex_area = freeSurfer_avg_vertex_area), 166L);
  testthat::expect_equal(fwhm2niters(20, avg_vertex_area = freeSurfer_avg_vertex_area), 294L);
  testthat::expect_equal(fwhm2niters(25, avg_vertex_area = freeSurfer_avg_vertex_area), 460L);
})


testthat::test_that("fwhm2niters computes the iterations from the actual mesh area when a surface is given.", {
  # When the mean vertex area is computed from the raw fsaverage mesh geometry
  # (and not from the FreeSurfer group-average normalization), the iteration
  # counts differ from the FreeSurfer values (e.g. 92 instead of 74 for
  # FWHM=10mm). This test pins the raw-mesh based values.
  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  # The mean vertex area computed from the raw mesh geometry:
  A = haze:::get.avg.vertex.area(mesh);
  testthat::expect_true(A > 0.3 && A < 0.5);   # ~0.4 mm^2 for fsaverage

  testthat::expect_equal(fwhm2niters(10, surface = mesh), 92L);            # surface is converted internally
  testthat::expect_equal(fwhm2niters(10, avg_vertex_area = A), 92L);       # same result when passing the area directly
  testthat::expect_equal(fwhm2niters(25, avg_vertex_area = A), 578L);

  # The inverse should give back roughly the requested FWHM.
  testthat::expect_equal(niters2fwhm(92, avg_vertex_area = A), 10, tolerance = 0.05);
})


testthat::test_that("niters2fwhm is the inverse of fwhm2niters (round trip).", {
  freeSurfer_avg_vertex_area = 0.5015;
  for(fwhm in c(2, 5, 10, 15, 20, 25)) {
    niters = fwhm2niters(fwhm, avg_vertex_area = freeSurfer_avg_vertex_area);
    fwhm_back = niters2fwhm(niters, avg_vertex_area = freeSurfer_avg_vertex_area);
    testthat::expect_equal(fwhm_back, fwhm, tolerance = 0.05);
  }
})


testthat::test_that("fwhm2niters returns a positive integer and is monotonically increasing in FWHM.", {
  A = 0.5;
  n = fwhm2niters(10, avg_vertex_area = A);
  testthat::expect_true(is.integer(n));
  testthat::expect_true(n >= 1L);
  testthat::expect_true(fwhm2niters(20, avg_vertex_area = A) > fwhm2niters(10, avg_vertex_area = A));
})


testthat::test_that("fwhm2niters warns when k is not 1, and validates its inputs.", {
  # Warning for k != 1 (the FreeSurfer calibration is only valid for k=1).
  testthat::expect_warning(fwhm2niters(10, avg_vertex_area = 0.5, k = 2L));
  testthat::expect_warning(niters2fwhm(10, avg_vertex_area = 0.5, k = 3L));

  # Need either a surface or an explicit mean vertex area.
  testthat::expect_error(fwhm2niters(10), "surface.*avg_vertex_area");
  testthat::expect_error(niters2fwhm(10));

  # FWHM must be a positive number.
  testthat::expect_error(fwhm2niters(-5, avg_vertex_area = 0.5), "fwhm");
  testthat::expect_error(fwhm2niters("abc", avg_vertex_area = 0.5), "fwhm");

  # The mean vertex area must be positive.
  testthat::expect_error(fwhm2niters(10, avg_vertex_area = -1), "avg_vertex_area");

  # The number of iterations must be a positive integer.
  testthat::expect_error(niters2fwhm(0, avg_vertex_area = 0.5), "niters");
  testthat::expect_error(niters2fwhm(-3, avg_vertex_area = 0.5), "niters");
})
