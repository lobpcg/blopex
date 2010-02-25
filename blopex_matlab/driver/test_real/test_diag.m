% script for testing "blopex in matlab", real input
% operatorA is diag(1,2,...,n), inv(operatorA) is implemented as a preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_DIAG.M, real input ...\n');
fprintf('\nSolving standard HEVP, operatorA is real diagonal, deflation w.r.t. first 3 eigenvectors, inv(operatorA) as preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')


clear
global operatorA

n = 10000;
operatorA = spdiags((1:n)',0,n,n);
operatorB = [];
Y=eye(n,3);  %constraints
X=rand(n,3);

% measure time for solution with "blopex_matlab"
tic
% 'precondD.m' implements "inv(operatorA)"
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,'precondD',Y,1e-4,40,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)
lambda

% measure time for solution with "lobpcg.m"
tic
% 'precondD.m' implements "inv(operatorA)"
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,'precondD',Y,1e-4,40,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)
lambda2
