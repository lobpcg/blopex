
% n=40;
% e = ones(n,1);
% A = spdiags([i*e -2*e -i*e], -1:1, n, n);
name='C:/test_complex1';

% n=10;
% A=rand(n,n)+rand(n,n)*i;
% A=sparse(A+A');
% name='C:/test_complex2';

% A=[1 1+i 0;1-i 2 0;0 0 3];
% A=sparse(A);
% name='C:/test_complex3';

PetscBinaryWrite(name,A)
A1=PetscBinaryRead(name,'complex');
% norm(A2-A,1)