This is a c program  (serial_driver) which was created for testing of the 
complex additions to blopex_abstract.  It is similar to blopex_serial_double but 
(1) deals only with complex matrices, 
(2) matmultivec, pcg_multi, and multivector routines have been modified to 
    handle complex matrices and vectors, 
(3) calls lobpcg_solve_complex, and 
(4) the serial_driver.c has been enhanced for additional test options.

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

To run the executible (which is in ../blopex_serial_double/driver), shift to that directory and enter command  

>> ./serial_driver R 100 10 

or

>> ./serial_driver L test1a.txt test1x.txt

The first form creates a 100x100 random hermitian matrix and solves for 10 eigenvalues.

The second form loads matrix from test1a.txt, initial eigenvectors from test1x.txt and 
solves for the the same number of eigenvalues as initial eigenvectors.  

These files are formated as follows: ascii lines with cr.  
1st line is '%d %d\n'  number of rows, number of columns. 
All other lines, '%22.16e %22.16e\n'   real part of number, imaginary part of number. 
Numbers are in standard fortran row rank order.  

These files can be created in matlab by 1st setting up the matrix and saving using the 
following (matrix in variable x and x is complex).

fid = fopen('c:\text1a.txt','w');
[row,col]=size(x);
fprintf(fid,'%d %d\n',row,col);
for j=1:4 
   for i=1:6 
      fprintf(fid,'%22.16e   %22.16e\n',real(x(i,j)),imag(x(i,j)));
   end
end
fclose(fid);


The matrix to solve should be hermitian positive definite.  
The matrix of initial values should have the same number of rows as the matrix to solve and
 the number of columns determines the number of eivenvalues to solve for. 


