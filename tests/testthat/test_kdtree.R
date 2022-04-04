
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


# This test shows how one could map data from one subject to another. The prerequisite is that
# you have spherical, aligned meshes for both subjects. This is the case if you have run
# FreeSurfer's reconall on the subjects: use the files <subject>/surf/lh.sphere and rh.sphere.
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
  testthat::expect_equal(length(dest_pervertex_data), nrow(dest_mesh$vertices));

  # The tests above do not really show that the mapped data makes sense. It's hard to test that,
  # but it is very easy to see it. Here is how you can check:
  do_plot = FALSE;
  if(do_plot) {
    if(requireNamespace("fsbrain", quietly = TRUE)) {
      fsbrain::vis.fs.surface(source_mesh, per_vertex_data = source_pervertex_data);
      fsbrain::vis.fs.surface(dest_mesh, per_vertex_data = dest_pervertex_data);
    }
  }
})

# This test shows how one could map data from one subject to another. The prerequisite is that
# you have spherical, aligned meshes for both subjects. This is the case if you have run
# FreeSurfer's reconall on the subjects: use the files <subject>/surf/lh.sphere and rh.sphere.
testthat::test_that("The R and CPP versions of the linear_interpolate_kdtree function give (almost) identical results.", {
  source_mesh_file = system.file("extdata", "fsaverage_mesh_lh_sphere", package = "haze", mustWork = TRUE);
  dest_mesh_file = system.file("extdata", "fsaverage6_mesh_lh_sphere", package = "haze", mustWork = TRUE);
  source_data_file = system.file("extdata", "fsaverage_lh_thickness", package = "haze", mustWork = TRUE);

  source_mesh = freesurferformats::read.fs.surface(source_mesh_file);
  source_pervertex_data = freesurferformats::read.fs.morph(source_data_file); # for source mesh
  dest_mesh = freesurferformats::read.fs.surface(dest_mesh_file);


  interp_res_R = linear_interpolate_kdtree(dest_mesh$vertices, source_mesh, source_pervertex_data, cpp = FALSE);
  interp_res_CPP = linear_interpolate_kdtree(dest_mesh$vertices, source_mesh, source_pervertex_data, cpp = TRUE);
  dest_pervertex_data_R = interp_res_R$interp_values;
  dest_pervertex_data_CPP = interp_res_CPP$interp_values;

  testthat::expect_true(is.vector(dest_pervertex_data_R));
  testthat::expect_true(is.double(dest_pervertex_data_R));
  testthat::expect_equal(length(dest_pervertex_data_R), nrow(dest_mesh$vertices));

  testthat::expect_true(is.vector(dest_pervertex_data_CPP));
  testthat::expect_true(is.double(dest_pervertex_data_CPP));
  testthat::expect_equal(length(dest_pervertex_data_CPP), nrow(dest_mesh$vertices));

  # The min values should be zero for both.
  testthat::expect_equal(min(dest_pervertex_data_R), 0.0);
  testthat::expect_equal(min(dest_pervertex_data_CPP), 0.0);

  # The max values are not identical for some reason (precision?), but should be in same range.
  testthat::expect_true(max(dest_pervertex_data_R) > 4.6 & max(dest_pervertex_data_R) < 4.8);
  testthat::expect_true(max(dest_pervertex_data_CPP) > 4.6 & max(dest_pervertex_data_CPP) < 4.8);

  # The real test for similar values:
  testthat::expect_equal(dest_pervertex_data_R, dest_pervertex_data_CPP, tolerance=1e-1);

  # We can also visually show that the data look similar:
  do_plot = FALSE;
  if(do_plot) {
    if(requireNamespace("fsbrain", quietly = TRUE)) {
      fsbrain::vis.fs.surface(dest_mesh, per_vertex_data = dest_pervertex_data_R);
      fsbrain::vis.fs.surface(dest_mesh, per_vertex_data = dest_pervertex_data_CPP);
    }
  }
})




