#!/bin/bash

# Apply this script to develop.20050606 Hypre version for testing.
# It is not guaranteed to work with other versions. 

# Ellen, replace the following files in krylov:

rm -f $1/krylov/lobpcg.*
cp ../blopex_abstract/krylov/lobpcg.* $1/krylov

rm -f $1/krylov/HYPRE_lobpcg.*
cp krylov/HYPRE_lobpcg.* $1/krylov

# Ellen, the file below is new. Please add it to krylov: 

cp krylov/HYPRE_MatvecFunctions.h $1/krylov

# Ellen, the same file as above. We put it in hypre/include
# by hand just for our testing. It will need to be copied
# to hypre/include by Hypre's make in the actual hypre distibution. 

mkdir $1/hypre         #gives a warning if the directory already exists 
mkdir $1/hypre/include #gives a warning if the directory already exists
cp krylov/HYPRE_MatvecFunctions.h $1/hypre/include

# Ellen, please remove in multivector the file:  

# rm $1/multivector/HYPRE_interpreter.h

# It is our file and we no longer need it. 
 

# Ellen, replace the following files in multivector: 

rm -f $1/multivector/multivector.*
cp ../blopex_abstract/multivector/multivector.* $1/multivector

rm -f $1/multivector/temp_multivector.*
cp ../blopex_abstract/multivector/temp_multivector.* $1/multivector

cp ../blopex_abstract/multivector/interpreter.h $1/multivector

# Ellen, replace the following files in utilities and in *_ls: 

rm -f $1/utilities/fortran_matrix.*
cp ../blopex_abstract/utilities/fortran_matrix.* $1/utilities

rm -f $1/parcsr_ls/HYPRE_parcsr_int.*
cp parcsr_ls/HYPRE_parcsr_int.* $1/parcsr_ls

rm -f $1/struct_ls/HYPRE_struct_int.*
cp struct_ls/HYPRE_struct_int.* $1/struct_ls

rm -f $1/sstruct_ls/HYPRE_sstruct_int.*
cp sstruct_ls/HYPRE_sstruct_int.* $1/sstruct_ls

# Ellen, we also need to modify the drivers, below, 
# but please check that the current drivers were not updated since 
# develop.20050606.tar.gz that we used here, before the copying.  

rm -f $1/test/ij.c
cp test/ij.c $1/test

rm -f $1/test/struct.c
cp test/struct.c $1/test

rm -f $1/test/sstruct.c
cp test/sstruct.c $1/test

