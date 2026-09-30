%LOBPCG_EXAMPLE Compute the four smallest eigenpairs of a 2-D Laplacian.
%
% Run this script from this folder, or add this folder to the MATLAB path.

gridSize = 20;
blockSize = 4;
oneDimensionalLaplacian = spdiags(...
    [-ones(gridSize, 1), 2 * ones(gridSize, 1), -ones(gridSize, 1)], ...
    -1:1, gridSize, gridSize);
operatorA = kron(speye(gridSize), oneDimensionalLaplacian) + ...
    kron(oneDimensionalLaplacian, speye(gridSize));
preconditioner = ichol(operatorA, struct('michol', 'on'));
operatorT = @(vectors) preconditioner \ (preconditioner' \ vectors);

rng(0);
initialVectors = randn(size(operatorA, 1), blockSize);
[eigenvectors, eigenvalues, failureFlag, ~, residualNormsHistory] = ...
    lobpcg(initialVectors, operatorA, [], operatorT, 1e-8, 100, 0);

fprintf('Convergence flag: %d (0 means all eigenpairs converged)\n', failureFlag);
disp('Smallest eigenvalues:');
disp(eigenvalues.');
fprintf('Largest final residual norm: %.3e\n', ...
    max(residualNormsHistory(:, end)));
