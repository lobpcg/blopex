% script for testing "blopex in matlab", complex arithmetic
% operatorA is dense complex Hermitian, diagonal preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_RANDOM_DENSE_COMPLEX.M, complex input ...\n');
fprintf('\nSolving standard HEVP, operatorA is random dense complex Hermitian, random complex initial guess, diagonal preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')

global operatorA
n = 500;
blockSize = 5;
% generate random complex Hermitian matrix with spectrum
% distributed uniformly over the interval (1e-1,10)  
a = 1e-1; 
b = 10;
D = a + (b-a).*rand(n,1);
V = orth(rand(n)+1i*rand(n)); 
operatorA = V*diag(D)*V';
clear('V');

operatorB = []; %not a generalized HEVP
Y=[];  % no constraints
X=rand(n,blockSize) + 1i*rand(n,blockSize);

% measure time for solution with "blopex_matlab"
tic
% 'precondD.m' implements "inv(diag(operatorA))"
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,'precondD',Y,1e-3,50,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
% 'precondD.m' implements "inv(diag(operatorA))"
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,'precondD',Y,1e-3,50,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)

% compare lobpcg/blopex solution with matlab's eig function
D = sort(eig(operatorA));
fprintf('\nEigenvalues computed by Matlab EIG:\n');
fprintf('\n%17.16e',D(1:blockSize));
fprintf('\n');
