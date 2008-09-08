% name='C:/L-matrix-complex.petsc';
% Construct 1000x1000 Hermitian Pos Def Matrix
% Random values 0-100
[A,err]=hpd_matrix(1000,.01);
if err > 0;
   err
   break;
end
A=A*100;
[i,j,v]=find(A);
nz=nnz(A);
[m,n]=size(A);

% Now force it to Sparce Complex
B=sparse(i,j,complex(v,0),m,n,nz);

% Write it to file, then read back in and compare
name='C:/test_complex4'; 
PetscBinaryWrite(name,B);
C=PetscBinaryRead(name,'complex');
norm(B-C,1)
 