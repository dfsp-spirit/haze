
testthat::test_that("The demo integer vv file can be read.", {
  vvfile = system.file("extdata", "fsaverage6_lh_white_meshdist_edge_1.vv", package = "haze", mustWork = TRUE);
  vv = read.vv(vvfile);

  testthat::expect_true(is.list(vv));
  testthat::expect_equal(length(vv), 40962L);

  neighbors_per_vertex = unlist(lapply(vv, length));
  testthat::expect_equal(min(neighbors_per_vertex), 6L);
  testthat::expect_equal(max(neighbors_per_vertex), 7L);
})
