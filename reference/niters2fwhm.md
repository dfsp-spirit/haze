# Convert number of NN smoothing iterations to a FreeSurfer Gaussian smoothing FWHM.

Inverse of
[`fwhm2niters`](https://dfsp-spirit.github.io/haze/reference/fwhm2niters.md):
convert a number of nearest-neighbor (NN) smoothing iterations into the
equivalent Gaussian smoothing kernel width, given as the full width at
half maximum (FWHM) in millimeters.

## Usage

``` r
niters2fwhm(niters, surface = NULL, avg_vertex_area = NULL, k = 1L)
```

## Arguments

- niters:

  positive integer, the number of NN smoothing iterations (the
  `num_iter` parameter of
  [`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md)).

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

numeric scalar, the FWHM in millimeters of the Gaussian smoothing kernel
that the given number of NN smoothing iterations approximates.

## Details

This is the inverse of
[`fwhm2niters`](https://dfsp-spirit.github.io/haze/reference/fwhm2niters.md),
using the same empirically calibrated FreeSurfer formula
(utils/mrisutils.cpp in the FreeSurfer source):


      gstd = sqrt((7 * avg_vertex_area * niters) / (1.14 * 4 * pi))
      fwhm = gstd * sqrt(log(256))

It is useful to check what amount of Gaussian smoothing (FWHM) a given
number of NN smoothing iterations corresponds to, e.g., to report the
equivalent FWHM in a paper or to compare settings across methods. The
calibration is only valid for the 1-ring neighborhood (`k=1`).

## See also

[`fwhm2niters`](https://dfsp-spirit.github.io/haze/reference/fwhm2niters.md)
for the forward conversion,
[`pervertexdata.smoothnn`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.md)
and
[`pervertexdata.smoothnn.adj`](https://dfsp-spirit.github.io/haze/reference/pervertexdata.smoothnn.adj.md)
for the actual NN smoothing.

## Examples

``` r
if (FALSE) { # \dontrun{
mesh = freesurferformats::read.fs.surface(system.file("extdata", "fsaverage_mesh_lh_white", package = "haze", mustWork = TRUE));
fwhm = niters2fwhm(92, surface = mesh);       # ~10 mm FWHM
} # }
```
