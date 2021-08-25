
test_that("The demo integer vv file can be read.", {
  vvfile = system.file("extdata", "fsaverage6_lh_white_meshdist_edge_1.vv", package = "smoothr", mustWork = TRUE);
  vv = read.vv(vvfile);

  testthat::expect_true(is.list(vv));
  testthat::expect_equal(length(vv), 40962L);
})
