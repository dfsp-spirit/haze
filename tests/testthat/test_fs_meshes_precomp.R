
test_that("The 1-ring graph neighborhood for fsaverage lh can be retrieved.", {
  vv = mesh.neigh.pre("lh_fsaverage");
  testthat::expect_true(is.list(vv));
  testthat::expect_equal(length(vv), 163842L);
})


test_that("The 1-ring graph neighborhood for fsaverage rh can be retrieved.", {
  vv = mesh.neigh.pre("rh_fsaverage");
  testthat::expect_true(is.list(vv));
  testthat::expect_equal(length(vv), 163842L);
})


test_that("The 1-ring graph neighborhood for fsaverage6 lh can be retrieved.", {
  vv = mesh.neigh.pre("lh_fsaverage6");
  testthat::expect_true(is.list(vv));
  testthat::expect_equal(length(vv), 40962L);
})


test_that("The 1-ring graph neighborhood for fsaverage6 rh can be retrieved.", {
  vv = mesh.neigh.pre("rh_fsaverage6");
  testthat::expect_true(is.list(vv));
  testthat::expect_equal(length(vv), 40962L);
})
