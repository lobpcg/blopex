Starting with version 1.1 of BLOPEX, source code for BLOPEX and
documentation is maintained at Google code site http://code.google.com/p/blopex/ .

This site contains Wiki's which document testing on various systems and a tar
of test files.  

The Wiki test documents are for testing of BLOPEX changes prior to incorporation into
the distribution tar for Petsc.  

To configure BLOPEX as distributed use option "--download-blopex=1" rather than the 
tars as described in the Wiki's.

Version 1.1 of BLOPEX adds support for complex numbers and 64bit integers.

There are 3 drivers supplied with the Petsc distribution;  driver, driver_fiedler,
and driver_diag.  These provide examples of using BLOPEX eigenvalue solver with 
Petsc.

Refer to Wiki's for details of configuration, driver compilation, and testing.

The test files for driver_fiedler are not part of the Petsc distribution but can 
be obtained by by downloading  "Petsc test files" from the Google code site. 

To install the test files: 
(1) Download "blopex_petsc_test.tar.gz" to your home directory
(2) cd $PETSC_DIR
(3) tar -zxvf $HOME/blopex_petsc_test.tar.gz

This places the test files in $PETSC_DIR/src/contrib/blopex/driver_fiedler.
Again, refer to the Wiki's for details of executing these tests.
