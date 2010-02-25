% script for testing "blopex in matlab", complex input
%attempts to compute several first eigenpairs of generalized HEVP where operatorA and operatorB come from
%2D Laplacian discretization using bilinear quadrilateral FE's  
%in a 1x1 square with the mesh size 1/100, under imaginary  Hermitian perturbations; incomplete Cholesky as preconditioner

clear
fprintf('\n==================================================================================================================================\n')
fprintf('\nRunning TEST_LAPLACE2D_PERTURBED_COMPLEX.M, real input ...\n');
fprintf('\nComputing several first eigenpairs of the 2D Laplacian discretized using bilinear quadrilateral FEs in 1x1 square\n');
fprintf('with the mesh size 1/100 under small imaginary Hermitian perturbations; random complex initial guess, incomplete Cholesky as preconditioner ...\n')
fprintf('\n==================================================================================================================================\n')

global R_cholinc % Cholesky factor for use in precondIC.m

blockSize = 10;
[operatorA,operatorB] = FEMDiscretizeSq(100,1);

%add imaginary Hermitian perturbations to operatorA and operatorB 
operatorA_imag = sprandsym(operatorA) .* 1e-10i; 
operatorA = operatorA + tril(operatorA_imag,-1) + conj(triu(operatorA_imag,1));
clear('operatorA_imag');
operatorB_imag = sprandsym(operatorB) .* 1e-10i;
operatorB = operatorB + tril(operatorB_imag,-1) + conj(triu(operatorB_imag,1));
clear('operatorB_imag');
R_cholinc=cholinc(operatorA,1e-3);

Y=[];  % no constraints
X=rand(size(operatorA,1),blockSize) + 1i*rand(size(operatorA,1),blockSize);

% measure time for solution with "blopex_matlab"
tic
% precondIC.m implements incomplete Cholesky as a preconditioner
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
  blopex_matlab(X,operatorA,operatorB,'precondIC',Y,1e-2,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, blopex_matlab: %g\n',time);
disp(s)


% measure time for solution with "lobpcg.m"
tic
% precondIC.m implements incomplete Cholesky as a preconditioner
[blockVectorX2,lambda2,failureFlag2,lambdaHistory2,residualNormsHistory2]=...
  lobpcg(X,operatorA,operatorB,'precondIC',Y,1e-2,70,1);
time = toc;
s = sprintf('\nSOLUTION TIME, lobpcg.m: %g\n',time);
disp(s)
clear('R_cholinc');
