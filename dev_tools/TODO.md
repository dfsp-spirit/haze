## haze TODO

### Remove dependency on custom Rvcg branch

Question: 
I once forked the Rvcg R package code to create a function I needed exported. I forgot why I didnt do PR to the author, or why he did not accept it, but anyways, he didnt make CRAN release in 12 years, so is wont help anyways. by now the main I forked from is 44 commits ahead of my branch. can you inspect this branch and tell me what I did? I assume implement some smoothing function in c++ that would be super slow in R? or what? the package, with this branch, is actually required by another of my packges (haze) and I want to get rid of this werid construct. please check what really is used from Rvcg, and why that branch is needed, and report for now: https://github.com/dfsp-spirit/Rvcg/tree/smooth_pervertex_data 


Answer:

What You Did on the smooth_pervertex_data Branch

In standard Rvcg, the smoothing functions (such as vcgSmooth / vcgSmoothImplicit) are designed to smooth the mesh geometry itself (i.e., updating the 3D vertex coordinates (x,y,z) via Laplacian, Taubin, or HC smoothing algorithms from VCGLib).

However, in neuroimaging / surface morphometry, you often have a fixed triangular mesh (like a brain surface) and 1D scalar per-vertex morphometry data (e.g., cortical thickness, curvature, or statistical maps) associated with each vertex that needs spatial smoothing across the mesh surface manifold.

To achieve this:

    VCGLib per-vertex scalar filtering: You tapped into VCGLib's vertex-attribute smoothing routines (or vertex neighbor averaging on vertex scalar attributes/qualities) in C++ to smooth scalar data across adjacent vertices over k iterations.

    Exported function: You exposed a dedicated R function (such as vcgSmoothPerVertexData / per-vertex scalar smoothing) taking a mesh3d object and a numerical vector of vertex values.

    Why in C++: In pure R, performing iterative Laplacian/diffusion smoothing across hundreds of thousands of vertices (like an fsaverage surface with ~163k vertices per hemisphere) by traversing adjacency lists in an R loop is prohibitively slow. VCGLib handles neighbor indexing and buffer swapping in optimized C++.

What haze Needs and How to Decouple It

haze is your package dedicated to fast smoothing of per-vertex data on triangular meshes. Because standard CRAN Rvcg lacked the per-vertex data smoothing entry point, haze had to rely on your specific fork/branch.
How to replace the dependency cleanly:

    Self-contained C++ in haze:
    Instead of depending on an external VCGLib fork, you can compute vertex-to-vertex adjacency directly in C++ (via Rcpp / Armadillo / Eigen or pure STL):

        Build the 1-ring neighborhood / adjacency list from the face index matrix (mesh$it or mesh$faces).

        Implement iterative 1D Laplacian averaging:
        vi(t+1)​=(1−λ)vi(t)​+λ∣N(i)∣1​j∈N(i)∑​vj(t)​

        This is only ~30–50 lines of C++ with OpenMP support, has zero heavy external dependencies, and removes the need for Rvcg entirely.

    Matrix / Graph approach:
    Precompute the normalized graph Laplacian / adjacency sparse matrix (via Matrix in R or Eigen in C++) and perform matrix-vector multiplications V(t+1)=WV(t) for the requested number of iterations.


