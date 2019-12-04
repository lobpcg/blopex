BLOPEX tools avialable in MATLAB/Octave
===============================================================

laplacian
-------
[![View Laplacian in 1D, 2D, or 3D on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d)

The code computes the exact eigenpairs of (1-3)D negative Laplacian on a rectangular finite-difference grid for combinations of Dirichlet, Neumann, and Periodic boundary conditions using explicit formulas from
http://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_of_the_second_derivative
The code can also compute the sparse matrix itself, using Kronecker sums of 1D Laplacians. For more information on tensor sums, see
http://en.wikipedia.org/wiki/Kronecker_sum_of_discrete_Laplacians

Example, compute everything for 3D negative Laplacian with mixed boundary conditions:
[lambda,V,A] = laplacian([100,45,55],{'DD' 'NN' 'P'}, 20);
Compute only the eigenvalues:
lambda = laplacian([100,45,55],{'DD' 'NN' 'P'}, 20);
Compute the matrix only:
[~,~,A] = laplacian([100,45,55],{'DD' 'NN' 'P'});

GNU OCTAVE compatible.

This code is a part of the BLOPEX eigensolver package, see
http://en.wikipedia.org/wiki/BLOPEX
or go directly to
http://code.google.com/p/blopex/

Copyright owners: Bryan C. Smith and Andrew V. Knyazev
Cite As
Andrew Knyazev (2019). Laplacian in 1D, 2D, or 3D (https://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d), MATLAB Central File Exchange. Retrieved December 4, 2019.


lobpcg
-------
[![View Locally Optimal Block Preconditioned Conjugate Gradient on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/48-locally-optimal-block-preconditioned-conjugate-gradient)

This main function LOBPCG is a version of the preconditioned conjugate gradient method (Algorithm 5.1) described in A. V. Knyazev, Toward the Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method, SIAM Journal on Scientific Computing 23 (2001), no. 2, pp. 517-541. http://dx.doi.org/10.1137/S1064827500366124
A C-version of this code is a part of the https://github.com/lobpcg/blopex
package and is available, e.g., in SLEPc and HYPRE. A scipy version is https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html
Tested in MATLAB 6.5-7.13-R2019a and available Octave 3.2.3-3.4.2.
Cite As
Andrew Knyazev (2019). Locally Optimal Block Preconditioned Conjugate Gradient (https://www.github.com/lobpcg/blopex), GitHub. Retrieved December 4, 2019.

matlab2hypre
-------
[![View matlab2hypre and hypre2matlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/27437-matlab2hypre-and-hypre2matlab)



MATLAB sparse matrix converter to Hypre
MATLAB vector converter to Hypre
Hypre sparse matrix converter to MATLAB
Hypre vector converter to MATLAB

The Parallel High Performance Preconditioners (hypre) is a library of routines for scalable (parallel) solution of linear systems. The built-in BLOPEX package in addition allows solving eigenvalue problems. See http://en.wikipedia.org/wiki/Hypre and http://en.wikipedia.org/wiki/BLOPEX

Hypre currently does not have a direct MATLAB interface. These codes allow one to transfer data between hypre and MATLAB/OCTAVE.

Our previous code,
http://www.mathworks.com/matlabcentral/fileexchange/7421-matlab-to-hypre-and-hypre-to-matlab-matrix-and-vector-converter
is not Hypre 2.6.0b compatible.

GNU OCTAVE compatible. Tested with Octave 3.2.3 and Hypre 2.6.0b.

This code is a part of the BLOPEX package:
http://en.wikipedia.org/wiki/BLOPEX or directly http://code.google.com/p/blopex/

Copyright 2010 Bryan C. Smith, Diana Zakaryan, Andrew V. Knyazev
Cite As

Andrew Knyazev (2019). matlab2hypre and hypre2matlab (https://www.mathworks.com/matlabcentral/fileexchange/27437-matlab2hypre-and-hypre2matlab), MATLAB Central File Exchange. Retrieved December 4, 2019.
