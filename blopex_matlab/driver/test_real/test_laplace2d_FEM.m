% script for testing "blopex in matlab", real input
%attempts to compute several first eigenpairs of 2D Laplacian discretized using bilinear quadrilateral FE's  
%in a 1x1 square with the mesh size 1/100; diagonal preconditioner
%for testing purposes operatorA (stiffness matrix) and operatorB (mass matrix) are implemented as functions 'multA.m' and 'multB.m' correspondingly

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_LAPLACE2D_FEM.M, real input ...\n');
fprintf('\nComputing several first eigenpairs of the 2D Laplacian discretized using bilinear quadrilateral FEs in 1x1 square\n');
fprintf('with the mesh size 1/100; random initial guess, diagonal  preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')

global operatorA
global operatorB
blockSize = 10;

[operatorA,operatorB] = FEMDiscretizeSq(100,1);

Y=[];  % no constraints
X=rand(size(operatorA,1),blockSize);

% measure time for solution with "blopex_matlab"
tic
%'precondD.m' implements "inv(diag(operatorA))"
% 'multA.m' and 'multB.m' implement multiplication by operatorA and operatorB respectively
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,'multA','multB','precondD',Y,1e-2,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
%'precondD.m' implements "inv(diag(operatorA))"
% 'multA.m' and 'multB.m' implement multiplication by operatorA and operatorB respectively
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,'multA','multB','precondD',Y,1e-2,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)

