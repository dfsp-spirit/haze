# Convert FreeSurfer Gaussian smoothing FWHM to number of NN smoothing iterations.

Convert a Gaussian smoothing kernel width, given as the full width at
half maximum (FWHM) in millimeters as used by FreeSurfer, into the
number of nearest-neighbor (NN) smoothing iterations to use with
[`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md)
(or
[`pervertexdata.smoothnn.adj`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.adj.md))
to achieve a similar amount of smoothing.

## Usage

``` r
fwhm2niters(fwhm, surface = NULL, avg_vertex_area = NULL, k = 1L)
```

## Arguments

- fwhm:

  numeric scalar, the desired full width at half maximum (FWHM) of the
  Gaussian smoothing kernel in millimeters, as used by FreeSurfer's
  `mris_fwhm --fwhm <fwhm>` or `mri_surf2surf --fwhm <fwhm>`.

- surface:

  a mesh, represented as an `fs.surface` instance from the
  `freesurferformats` package, a `tmesh3d` instance from `rgl`, or a
  character string representing the path of a surface file to load. Used
  to compute the mean vertex area. Ignored if `avg_vertex_area` is
  given.

- avg_vertex_area:

  numeric scalar, the mean vertex area of the mesh in square millimeters
  (squared mesh units). If `NULL` (the default), it is computed from
  `surface`. See details for how FreeSurfer's own numbers differ.

- k:

  positive integer, the k-ring neighborhood size to use for NN
  smoothing. The FreeSurfer-based calibration below is only valid for
  the 1-ring neighborhood (`k=1`); a warning is issued otherwise.

## Value

integer scalar, the number of NN smoothing iterations (the `num_iter`
parameter of
[`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md))
that approximates Gaussian smoothing with the requested FWHM.

## Details

True Gaussian smoothing on a triangular mesh would require computing
geodesic distances between (many) vertex pairs, which is very slow.
FreeSurfer avoids this by emulating the Gaussian with many iterations of
nearest-neighbor (1-ring) averaging: the whole trick is that *repeated
local averaging converges to a Gaussian kernel* (central limit theorem
on a graph), so no geodesic distances are needed at all. Each iteration
is a cheap local operation (average a vertex with its direct edge
neighbors), which is exactly what
[`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md)
implements.

The relationship between the effective Gaussian standard deviation
`gstd` after `n` iterations and the mean vertex area `avg_vertex_area`
was calibrated empirically by FreeSurfer in `MRISfwhm2niters()` (see
utils/mrisutils.cpp in the FreeSurfer source):


      gstd   = fwhm / sqrt(log(256))                     # FWHM -> Gaussian sigma (log(256) = 8*ln(2))
      niters = round(1.14 * (4 * pi * gstd^2) / (7 * avg_vertex_area))

where `7` approximates the average number of neighbors in the 1-ring of
a triangulated cortical surface, and `1.14` is an empirical fudge factor
that FreeSurfer fitted so that the measured FWHM of the smoothed output
matches the requested one. This function is the inverse of
[`niters2fwhm`](https://dfsp-spirit.github.io/haze/reference/niters2fwhm.md).

**Note on matching FreeSurfer's own outputs:** FreeSurfer normalizes the
mean vertex area by the *group-average* surface area stored in its
template surface files (`group_avg_surface_area`, about 82,000 mm^2 per
hemisphere for fsaverage), which is larger than the raw geometric area
of the fsaverage meshes (about 65,000 mm^2 for the packaged fsaverage lh
white surface). FreeSurfer therefore reports somewhat lower iteration
counts than you get when `avg_vertex_area` is computed directly from
these meshes (e.g., for FWHM = 10 mm on fsaverage: 74 iterations in
FreeSurfer vs. 92 from the raw mesh area). To reproduce FreeSurfer's
exact numbers, pass
`avg_vertex_area = group_avg_surface_area / nvertices` (about 0.5015
mm^2 for the fsaverage templates).

## See also

[`niters2fwhm`](https://dfsp-spirit.github.io/haze/reference/niters2fwhm.md)
for the inverse conversion,
[`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md)
and
[`pervertexdata.smoothnn.adj`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.adj.md)
for the actual NN smoothing.

## Examples

``` r
if (FALSE) { # \dontrun{
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
niters = fwhm2niters(10, surface = mesh);      # 92 iterations, approximating ~10 mm FWHM
smoothed = pervertexdata.smoothnn(mesh, data, num_iter = niters);
} # }
```
