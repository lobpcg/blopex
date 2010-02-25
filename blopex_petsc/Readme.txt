
---- Setup Petsc ----  

To use blopex_solve_double or blopex_solve_complex using the 
Beta release of Petsc first download the development version from 
the Petsc website at ftp://info.mcs.anl.gov/pub/petsc-dev.tar.gz
 
To install files issue command:

>> gunzip -c petsc-dev.tar.gz | tar -xof -

This creates files in directory petsc-dev.  
Now two environment variables must be setup:

>> PETSC_DIR=/home/grads/dmccuan/petsc-dev
>> export PETSC_DIR
>> PETSC_ARCH=linux-gnu-c-debug
>> export PETSC_ARCH

We will also need these set when compiling the blopex_petsc interface.   
To avoid setting everytime you logon you can put them in your .bashrc 
file as follows:

export PETSC_DIR=/home/grads/dmccuan/petsc-dev
export PETSC_ARCH=linux-gnu-c-debug

The next step is to configure Petsc.  Shift to the Petsc directory and 
issue command

>> ./config/configure.py  --with-shared --with-scalar-type=complex

Since we are not using the source for blopex included in Petsc we don't need the 
options --download-blopex=1 or --download-hypre=1.  

Petsc runs with either double or complex objects but not both.  
To test the blopex routines for double configure Petsc without the 
option --with-scalar-type=complex. 


After configuration is finished issue the command

>> make all test

---- Setup blopex_petsc ----  

Source for blopex_abstract and blopex_petsc must be present on your system.
Either obtained via tarballs or svn repository.

In order to compile the blopex_petsc drivers for testing we must first have 
blopex_abstract installed and compiled.  

Then, some adjustments to the makefiles for blopex_petsc are needed. 

The ../blopex_petsc/Makefile.inc must specify the location of petsc_abstract
For example:

LOBPCG_ROOT_DIR = ../../blopex_abstract

when blopex_abstract and petsc_abstract directories are defined under the same directory. 

Then ../blopex_petsc/driver/Makefile, ../blopex_petsc/driver_fiedler/Makefile, and 
../blopex_petsc/petsc-interface must have an include for the Petsc variable definitions.   
For the Beta version use

include ${PETSC_DIR}/conf/variables

note:  ${CLINKER}, ${PETSC_COMPILE_SINGLE}, and ${PETSC_KSP_LIB} are set by this include. 

After running make, executables for driver and driver_fiedler are created in their 
respective directories.  

----- Execution of Tests -----

There are two test drivers located in subdirectories of blopex_petsc:
driver an driver_fiedler.

driver builds a 7pt laplacian for solution and calls either lobpcg_solve_complex 
if Petsc is configured for complex (this is controlled by PETSC_USE_COMPLEX preprocessor
 variable) or lobpcg_solve_doube if Petsc is configured for double. 

For example:

>> mpirun -np 2 driver -n_eigs 3 -itr 20 
or 
>> ./driver -n_eigs 3 -itr 20

driver_fiedler accepts as input matrices in Petsc format. 

For example:

>> ./driver_fiedler -matrix DL-matrix-complex.petsc -n_eigs 3 -itr 20

The option -matrix specifies the matrix to solve.  
The initial eigenvectors are generated randomly. 

The matrix file is in a Petsc format.  
These can be setup via some Matlab programs in the PETSc socket interface to Matlab; 
PetscBinaryRead.m and PetscBinaryWrite.m.  
These programs read and write Matlab matrices and vectors to files formated for Petsc. 
The version from Petsc only supports double.   
We have modified these programs to  also support complex.  
The complex versions are included in the .../blopex_petsc directory along with 
PetscWriteReadExample.m to illustrate how to use them. 

The double files are L-matrix-double.petsc (65536x65536 diagonal values 0-4)
and DL-matrix-double.petsc (65536x65536 tridiagonal values 0-4).  

The complex files are complex versions of the L and DL double matrices (with imag part zero) 
plus test_complex1.petsc (40x40 134 nz hpd random),
test_complex2.petsc (1000x1000 50786 nz hpd random) and 
test_complex3 (10876 nz hpd random). 
