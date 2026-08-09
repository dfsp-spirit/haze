# Read vv binary file.

Read matrix-like data from vv files. This is a custom format from the
cpp_geodesic repo, designed to allow fast reading of vector-of-vectors
data. The format does NOT require that all inner vectors have the same
length, so it is NOT limited to matrices. The format is designed for
storing graphs as adjacency lists.

## Usage

``` r
read.vv(filepath)
```

## Arguments

- filepath:

  string. Full path to the input vv file.

## Value

list of vectors, the data. The vv files may can store double or int,
which is encoded in the file header and used accordingly.
