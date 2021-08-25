# smoothr
Fast smoothing of per-vertex data on triangular meshes for R.


## About

This package performs smoothing of per-vertex data:

![Vis](./web/smoothr.jpg?raw=true "Per-vertex data on a brain mesh before (left) and after (right) smoothing.")
**Fig.1**: *Per-vertex data on a brain mesh before (left) and after (right) smoothing.*

To avoid any confusion: `smoothr` does **not** smooth the mesh itself (use the `Rvcg` R package for that).

## Features

* nearest neighbor smoothing based on edge distance (e.g., `k`-ring neighborhood of each vertex, with arbitrary `k`)
* Gaussian smoothing based on geodesic distances on the mesh (geodesic computation is slow for large meshes), WIP
* re-use of neighborhood data for smoothing several datasets on the same mesh
* comes with pre-computed neighborhood data for meshes commonly used in surface-based neuroimaging (FreeSurfer's fsaverage and fsaverage6)

## Usage

This is WIP, come back another day.

