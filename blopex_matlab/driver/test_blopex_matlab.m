% script for testing "blopex in matlab"

clear
n = 10000;
operatorA = spdiags((1:n)',0,n,n);
operatorB = speye(n,n);
Y=eye(n,3);  %constraints
X=rand(n,3);

% measure time for solution with "blopex_matlab"
tic
% 'precond.m' implements "inv(operatorA)"
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,'precond',Y,1e-4,40,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
% 'precond.m' implements "inv(operatorA)"
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,'precond',Y,1e-4,40,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)