BLOPEX: Block Locally Optimal Preconditioned Eigenvalue Xolvers
===============================================================

Authors: A. V. Knyazev, M. E. Argentati, I. Lashuk, E. E. Ovtchinnikov


General
-------

Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) is a package, written in C and MATLAB/OCTAVE, that includes an eigensolver implemented with the Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG). Its main features are: a matrix-free iterative method for computing several extreme eigenpairs of symmetric positive generalized eigenproblems; a user-defined symmetric positive preconditioner; robustness with respect to random initial approximations, variable preconditioners, and ill-conditioning of the stiffness matrix; and apparently optimal convergence speed.

BLOPEX supports parallel MPI-based computations. BLOPEX is incorporated in the HYPRE package and is available as an external block to the SLEPc package. [PHAML](http://math.nist.gov/~WMitchell/phaml/phaml.html) has an interface to call BLOPEX eigensolvers, as well as DevTools by [Visual Kinematics](http://vki.com).

How to use it
-------------

- This code does NOT compute ALL eigenvalues. It is designed to compute no more than ~20%, i.e., if the matrix size is 100, do NOT attempt to compute more that ~20 eigenpairs, otherwise, the code will crash. 

- We do NOT provide the binaries, only the source, mostly in C, some in FORTRAN and MATLAB/OCTAVE. 

- To get the source, checkout the repository. Do NOT use the "Downloads" tab --- the downloads are for developers only. 

- In HYPRE, the BLOPEX code is already incorporated. Just download and compile HYPRE. 

- In [SLEPc][SLEPc] release 3.2 and later, the BLOPEX code is available as an external package. Download and install [PETSc][PETSc] first. Then download SLEPc, and use the option `--download-blopex=1` in the SLEPc configure, prior to SLEPc make. 

- To build a stand-alone version, one needs the core of BLOPEX, which is `blopex_abstract`, which contains the eigensolver codes and can be compiled into the BLOPEX library, using the sample makefiles provided. One also needs `blopex_serial_double` and/or `blopex_serial_complex`. The current stand-alone version is for serial computations only. For parallel computations, one must use SLEPc or HYPRE. 

- A native MATLAB/Octave code is available at `blopex_tool/matlab/lobpcg`, see also [Matlab Central][matlabcentral1]. 

- A native Java code is available from a [sister project][sparse-eigensolvers-java].

- An interface of our C `blopex_abstract` code to MATLAB (both 32 and 64 bit) is available at `blopex_matlab`. This interface is mainly designed for testing of our C `blopex_abstract` code. For actual computations in MATLAB/OCTAVE, is it instead recommended to use the native MATLAB/OCTAVE code above. 

- `blopex_petsc` and `blopex_hypre` are provided for reference only. 

We also provide some useful testing tools in `blopex_tools/matlab`:

- `laplacian` builds the matrix of a 1-,2-, or 3-D [Laplacian][wikipedia1] and computes their eigenpairs using the [explicit formulas][wikipedia2]. The same code is also available [here][matlabcentral2].

- `matlab2hypre` converts matrices between HYPRE and MATLAB/OCTAVE. The same code is also available [here][matlabcentral3].

[SLEPc]:     http://slepc.upv.es
[PETSc]:     http://www.mcs.anl.gov/petsc
[sparse-eigensolvers-java]:  http://code.google.com/p/sparse-eigensolvers-java/ 
[wikipedia1]: http://en.wikipedia.org/wiki/Kronecker_sum_of_discrete_Laplacians
[wikipedia2]: http://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_of_the_second_derivative
[matlabcentral1]:  http://www.mathworks.com/matlabcentral/fileexchange/48-lobpcg-m
[matlabcentral2]:  http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d
[matlabcentral3]:  http://www.mathworks.com/matlabcentral/fileexchange/27437-matlab2hypre-and-hypre2matlab

BLOPEX examples
---------------

There are a few codes, used to test the library, which can be viewed as examples, called "drivers", see, e.g., the code for single-thread double precision at `blopex_serial_double/driver/serial_driver.c`.

HYPRE-BLOPEX has its own example, search HYPRE documentation for LOBPCG or eigenvalue. Several HYPRE drivers are used for testing in the BLOPEX 2007 paper listed below.

SLEPc has lots of examples. Use the option `--download-blopex=1` in the SLEPc configure, prior to SLEPc make and read SLEPc docs how to run SLEPc examples using BLOPEX eigenvalue solver LOBPCG.

Please contact HYPRE and SLEPc developers directly if you have any questions about their examples of using LOBPCG-BLOPEX. 

Can BLOPEX be made to run faster?
---------------------------------

Yes, BLOPEX can be made to run much faster. None of current BLOPEX implementations, including those in HYPRE and SLEPc, uses a "true multivector". A "multivector" is the main structure in BLOPEX, representing blocks of vectors. LOBPCG code requires performing basic algebraic operations, e.g., sums, scalar products, and application of functions, with multivectors. However, to our knowledge, in all current (2019) implementations of BLOPEX, the multivector is faked, i.e., is simply substituted by a collection of individual vectors. Thus, e.g., a sum of fake multivectors is a loop of sums of vectors, which is only Level 1 BLAS. A "true multivector" would likely accelerate the code by orders of magnitude in parallel runs.

We currently (2019) do not have resources and thus plans to implement the "true multivector" in serial versions of BLOPEX, even though that would make the code run faster several times due to the use of Level 3 BLAS. If you volunteer to do it, we can help with advice, and would be glad to add you to the team of developers. It is not that difficult, but is surely time consuming and requires proper qualifications. That could be a fun project, e.g., for INTEL and AMD developers, to implement a true multivector library, make a BLOPEX interface to it, and test it on important applications.

HYPRE has an incomplete (as of 2019) implementation of the "true multivector," with some key functions still missing, and apparently no plans to complete it.

PETSc developers plan (Nov. 2019) to implement the "true multivector" in the next PETSc release in the format of `TAIJ` matrix that will allow `MatMult`, `MatSolve`, etc., using aggregated communication for distribution of the multivector between the processors and contiguous local structure suitable for high level BLAS. Please contact `petsc-maint@mcs.anl.gov` for questions/comments. 

Related Projects
----------------

- [HYPRE](http://www.llnl.gov/CASC/hypre).
- [SLEPc](http://slepc.upv.es) and [PETSc](http://www.mcs.anl.gov/petsc).
- [sparse-eigensolvers-java](http://code.google.com/p/sparse-eigensolvers-java).


References
----------

1. A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov, *Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in hypre and PETSc*, SIAM Journal on Scientific Computing 25(5): 2224-2239 (2007) [DOI](http://dx.doi.org/10.1137/060661624)

2. A. V. Knyazev, *Toward the Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method*, SIAM Journal on Scientific Computing 23(2): 517-541 (2001) [DOI](http://dx.doi.org/10.1137/S1064827500366124)


Code license
------------
MIT / Apache-2.0 License

Contributors
------------
Prashanth Nadukandi, Jose E. Roman

-----------------

This material was based upon work supported by the National Science Foundation under Grant No. 111 5734. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
