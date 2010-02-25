% script for testing "blopex in matlab", real input
%attempts to compute several first eigenpairs of the Poisson operator (5-point Laplacian) 
%in a 2x2 square with the mesh size 1/10 using incomplete Cholesky factorization as a preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_LAPLACE2D.M, real input ...\n');
fprintf('\nComputing several first eigenpairs of the Poisson operator (5-point Laplacian) in 2x2 square with the mesh size 1/10,\n');
fprintf(' random initial guess, incomplete Cholesky as preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')

global R_cholinc % Cholesky factor for use in precondIC.m
blockSize = 10;

operatorA = 100.*delsq(numgrid('S',21)); [n,n]=size(operatorA);
R_cholinc=cholinc(operatorA,1e-3);

operatorB = [];%not a generalized HEVP
Y=[];  % no constraints
X=rand(n,blockSize);

% measure time for solution with "blopex_matlab"
tic
% precondIC.m implements incomplete Cholesky as a preconditioner
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,'precondIC',Y,1e-5,50,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
% precondIC.m implements incomplete Cholesky as a preconditioner
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,'precondIC',Y,1e-5,50,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)

% compare lobpcg/blopex solution with matlab's eigs function
opts.issym=1;opts.isreal=1;opts.tol=1e-5;opts.maxit=50;opts.disp=0;
D = eigs(operatorA,blockSize,'SA',opts);
fprintf('\nEigenvalues computed by Matlab EIGS:\n');
fprintf('\n%17.16e',D(1:blockSize));
fprintf('\n');
clear('R_cholinc');
