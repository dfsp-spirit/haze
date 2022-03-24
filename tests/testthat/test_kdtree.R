
testthat::test_that("One can find the vertices closest to given query coordinates on a mesh.", {
  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  query_coordinates = matrix(c(-65, -102,-42, 1, 64,75), ncol=3L, byrow = TRUE);
  res = haze:::find_nv_kdtree(query_coordinates, mesh);

  testthat::expect_true(is.list(res));
  testthat::expect_true("index" %in% names((res)));
  testthat::expect_true("distance" %in% names((res)));

  testthat::expect_true(is.matrix(res$index));
  testthat::expect_equal(length(res$index), 2L);

  testthat::expect_true(is.matrix(res$distance));
  testthat::expect_true(is.double(res$distance));
  testthat::expect_equal(length(res$distance), 2L);
  # See the next test for a prove that the reported vertices are really the closest ones.
})


testthat::test_that("One can retrieve per-vertex data for mesh vertices closest to query coordinates.", {
  fsmesh_file = system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE);
  mesh = freesurferformats::read.fs.surface(fsmesh_file);

  num_verts = nrow(mesh$vertices);
  pervertex_data = rnorm(num_verts, 5.0, 1.0);

  source_vertices = c(50L, 100L, 111111L);
  query_coordinates = mesh$vertices[source_vertices,] + 0.001; # Change coords slightly, the kdtree should find the closest vertices then.

  idata = nn_interpolate_kdtree(query_coordinates, mesh, pervertex_data);

  testthat::expect_true(is.vector(idata));
  testthat::expect_true(is.double(idata));
  testthat::expect_equal(length(idata), length(source_vertices));
  testthat::expect_equal(idata[1], pervertex_data[source_vertices[1]]);
  testthat::expect_equal(idata[2], pervertex_data[source_vertices[2]]);
  testthat::expect_equal(idata[3], pervertex_data[source_vertices[3]]);
})


testthat::test_that("One can map per-vertex data between spherical meshes.", {
  source_mesh_file = system.file("extdata", "fsaverage_mesh_lh_sphere", package = "haze", mustWork = TRUE);
  dest_mesh_file = system.file("extdata", "fsaverage6_mesh_lh_sphere", package = "haze", mustWork = TRUE);
  source_data_file = system.file("extdata", "fsaverage_lh_thickness", package = "haze", mustWork = TRUE);

  source_mesh = freesurferformats::read.fs.surface(source_mesh_file);
  source_pervertex_data = freesurferformats::read.fs.morph(source_data_file); # for source mesh
  dest_mesh = freesurferformats::read.fs.surface(dest_mesh_file);


  interp_res = linear_interpolate_kdtree(dest_mesh$vertices, source_mesh, source_pervertex_data);
  dest_pervertex_data = interp_res$interp_values;

  testthat::expect_true(is.vector(dest_pervertex_data));
  testthat::expect_true(is.double(dest_pervertex_data));
  testthat::expect_equal(length(dest_pervertex_data), length(source_pervertex_data));
})



