This is a c program (serial_driver) which serves as an example and test of 
lobpcg_solve_double in the blopex_abstract package.

It calls lobpcg_solve_double with a 40x40 diagonal matrix (double precision) 
with values of 1 thru 40 on the diagonal.  

It solves for 5 eigenvalues, initial eigenvectors are random, tolerance 1e-6 and max iterations 100. 

This is the same serial_driver as before the addition of complex numbers.  
The only change is to call call lobpex_solve_double instead of logpcg_solve. 

To create the executable for serial_driver, the blopex_abstract objects must be compiled and 
Makefile.inc in ../blopex_serial_double updated with appropriate locations for blopex_abstract, 
compile flags, and lapack library.

For example: 

LOBPCG_ROOT_DIR =  ../../blopex_abstract
CFLAGS = -Wall -g
LAPACKLIB = -llapack

Then, issue command 

>> make  

This should create a test executable "driver/serial_driver" 

To run the executible (which is in ../blopex_serial_double/driver), 
shift to that directory and enter command  

>> ./serial_driver

