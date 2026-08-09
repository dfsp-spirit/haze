# Convert homogeneous coordinates to Cartesian coordinates.

Convert homogeneous coordinates to Cartesian coordinates.

## Usage

``` r
homogeneous_to_cartesian(homog)
```

## Arguments

- homog:

  nx4 numeric matrix of input coordinates

## Value

nx3 matrix of Cartesian coordinates

## Examples

``` r
if (FALSE) { # \dontrun{
homog = matrix(c(1,2,3,1,1,2,3,2), ncol=4, byrow=TRUE);
haze:::homogeneous_to_cartesian(homog);
} # }
```
