% script for testing "blopex in matlab", real input
% operatorA is sparse real symmetric, no preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_RANDOM_SPARSE_REAL.M, real input ...\n');
fprintf('\nSolving standard HEVP, operatorA is random sparse real symmetric, random initial guess, no preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')

n = 10000;
density = .1;
blockSize = 10;
% generate random sparse real symmetric matrix, nonzeroes ~ density*n^2
operatorA = sprandsym(n,density);

operatorB = [];%not a generalized HEVP
Y=[];  % no constraints
X=rand(n,blockSize);

% measure time for solution with "blopex_matlab"
tic
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,Y,1e-2,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,Y,1e-2,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)

% compare lobpcg/blopex solution with matlab's eigs function
opts.issym=1;opts.isreal=1;opts.tol=1e-2;opts.maxit=70;opts.disp=0;
D = eigs(operatorA,blockSize,'SA',opts);
fprintf('\nEigenvalues computed by Matlab EIGS:\n');
fprintf('\n%17.16e',D(1:blockSize));
fprintf('\n');
