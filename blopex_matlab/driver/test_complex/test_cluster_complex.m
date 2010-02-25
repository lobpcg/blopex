% script for testing "blopex in matlab", complex arithmetic
% operatorA is complex Hermitian with multiple and clustered eigenvalues: 0,0,1e-20,1e-10,1e-5,1e-5 + 1e-10, 1, 1 + 1e-10, 1 + 1e-9, 1 + 1e-8; 
% diagonal preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_CLUSTER_COMPLEX.M, real input ...\n');
fprintf('\nSolving standard HEVP, operatorA is dense complex with clustered and multiple eigenvalues,\n');
fprintf(' random complex initial guess, diagonal preconditioner ...\n');
fprintf('\n==================================================================================================================================\n')

global operatorA
n = 100;
blockSize = 11;
num_clustered = 10;

% generate diagonal matrix with prescribed clustered and multiple num_clustered smallest eigenvalues
% the rest (n-num_clustered) diagonal elements are distributed randomly (uniformly) on [2 ... 100]
D = zeros(n,1);
D(3) = 1e-20;
D(4) = 1e-10;
D(5) = 1e-5;
D(6) = 1e-5 + 1e-10;
D(7) = 1.;
D(8) = 1 + 1e-10;
D(9) = 1 + 1e-9;
D(10)= 1 + 1e-8; 
D(num_clustered+1:n) = 2 + 98.*rand(n-num_clustered,1);
V = orth(rand(n)+1i*rand(n)); 
operatorA = V*diag(D)*V';
clear('V');

operatorB = [];%not a generalized HEVP
Y=[];  % no constraints
X=rand(n,blockSize);

% measure time for solution with "blopex_matlab"
tic
% 'precondD.m' implements "inv(diag(operatorA))"
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,'precondD',Y,1e-10,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
% 'precondD.m' implements "inv(diag(operatorA))"
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,'precondD',Y,1e-10,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)

% compare lobpcg/blopex solution with matlab's eig function
D = sort(eig(operatorA));
fprintf('\nEigenvalues computed by Matlab EIG:\n');
fprintf('\n%17.16e',D(1:blockSize));
fprintf('\n');
