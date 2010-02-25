% script for testing "blopex in matlab", complex arithmetic
% operatorA is sparse complex Hermitian, no preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_RANDOM_SPARSE_COMPLEX.M, complex input ...\n');
fprintf('\nSolving standard HEVP, operatorA is random sparse complex Hermitian, random complex initial guess, no preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')

n = 3000;
density = .1;
blockSize = 7;
% generate random sparse complex Hermitian matrix, nonzeroes ~ density*n^2
operatorA = sprandsym(n,density);
operatorA_imag = sprandsym(operatorA).*1i; % generate imaginary part for operatorA with the same sparsity pattern
operatorA = operatorA + tril(operatorA_imag,-1) + conj(triu(operatorA_imag,1));
clear('operatorA_imag');

operatorB = []; % not a generalized HEVP
Y=[];  % no constraints
X=rand(n,blockSize) + 1i*rand(n,blockSize);

% measure time for solution with "blopex_matlab"
tic
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,Y,1e-2,50,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,Y,1e-2,50,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)

% compare lobpcg/blopex solution with matlab's eigs function
opts.issym=1;opts.isreal=0;opts.tol=1e-2;opts.maxit=50;opts.disp=0;
D = eigs(operatorA,blockSize,'SR',opts);
fprintf('\nEigenvalues computed by Matlab EIGS:\n');
fprintf('\n%17.16e',D(1:blockSize));
fprintf('\n');
