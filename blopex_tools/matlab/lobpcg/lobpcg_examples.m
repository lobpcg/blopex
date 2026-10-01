%LOBPCG_EXAMPLE Demonstrate the examples documented in lobpcg.m.
%
% Call lobpcg_examples from this folder, or add this folder to the MATLAB path.
% The codistributed example runs only when Parallel Computing Toolbox is
% installed and licensed. Examples that intentionally use difficult
% preconditioners or clustered spectra may not fully converge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   License:  MIT / Apache-2.0
%   Copyright (c) 2026 A.V. Knyazev, Andrew.Knyazev@ucdenver.edu
%   $Revision: 1.0 $  $Date: 1-November-2026
%   Tested in GNU Octave Version: 11.3.0 and
%   MATLAB with Paralel Computing Toolbx 26.2 (R2026b), but is
%   expected to work on any post 2008 MATLAB and Octave.
function lobpcg_examples

previousRngState = rng;
rngCleanup = onCleanup(@() rng(previousRngState));
rng(0);
exampleSettings = struct(...
    'gridSize', 20, 'laplacianBlockSize', 8, ...
    'standardTolerance', 1e-5, 'verbosity', 2, 'quiet', 0, ...
    'unpreconditionedMaxIterations', 50, ...
    'sequentialBlockSize', 2, 'sequentialBlockCount', 4, ...
    'sequentialMaxIterations', 200, ...
    'preconditionedMaxIterations', 60, 'identityBMaxIterations', 50, ...
    'diagonalSize', 1000, 'diagonalDensity', 0.1, ...
    'diagonalBlockSize', 5, 'diagonalMaxIterations', 15, ...
    'parallelSize', 100, 'parallelBlockSize', 2, ...
    'parallelWorkerCount', 2, 'parallelMaxIterations', 5, ...
    'singleSize', 100, 'singleBlockSize', 2, 'singleMaxIterations', 15, ...
    'clusteredMatrixSize', 100, 'clusteredValues', [1, 1.99, 2:99], ...
    'clusteredBlockSizes', [1, 2, 3], ...
    'clusteredTolerance', 1e-10, 'clusteredMaxIterations', 80, ...
    'nonsymmetricEntry', 2, 'selectiveShift', 3.5);
gridSize = exampleSettings.gridSize;
blockSize = exampleSettings.laplacianBlockSize;
oneDimensionalLaplacian = spdiags(...
    [-ones(gridSize, 1), 2 * ones(gridSize, 1), -ones(gridSize, 1)], ...
    -1:1, gridSize, gridSize);
operatorA = kron(speye(gridSize), oneDimensionalLaplacian) + ...
    kron(oneDimensionalLaplacian, speye(gridSize));
n = size(operatorA, 1);
fprintf('Using a %d-by-%d interior-grid Laplacian (%d unknowns).\n', ...
    gridSize, gridSize, n);
% Equivalent PDE-grid construction from the docstring, when available:
% operatorA = delsq(numgrid('S', gridSize + 1));
% Or, with the laplacian.m File Exchange function on the MATLAB path:
% [~, ~, operatorA] = laplacian([gridSize - 1, gridSize - 1]);

%% Example 1: Eight eigenpairs without preconditioning
[X_noPreconditioner, lambda_noPreconditioner, failure_noPreconditioner] = ...
    lobpcg(randn(n, blockSize), operatorA, exampleSettings.standardTolerance, ...
    exampleSettings.unpreconditionedMaxIterations, exampleSettings.verbosity);
fprintf('Unpreconditioned solve convergence flag: %d\n', ...
    failure_noPreconditioner);
fprintf('Returned %d eigenvectors.\n', size(X_noPreconditioner, 2));
disp('Unpreconditioned eigenvalues:');
disp(lambda_noPreconditioner.');

%% Example 2: Compute the same eigenpairs in four constrained blocks
blockSizeSequential = exampleSettings.sequentialBlockSize;
sequentialBlockCount = exampleSettings.sequentialBlockCount;
X_sequential = zeros(n, blockSizeSequential * sequentialBlockCount);
lambdaAll = zeros(blockSizeSequential * sequentialBlockCount, 1);
for blockIndex = 1:sequentialBlockCount
    [X_block, lambda_block] = lobpcg(randn(n, blockSizeSequential), ...
        operatorA, X_sequential(:, 1:(blockIndex - 1) * blockSizeSequential), ...
        exampleSettings.standardTolerance, ...
        exampleSettings.sequentialMaxIterations, exampleSettings.verbosity);
    blockRange = (blockIndex - 1) * blockSizeSequential + (1:blockSizeSequential);
    X_sequential(:, blockRange) = X_block;
    lambdaAll(blockRange) = lambda_block;
end
fprintf('Sequential solve returned %d eigenvectors.\n', size(X_sequential, 2));
disp('Eigenvalues from sequential constrained blocks:');
disp(lambdaAll.');

%% Example 3: Incomplete-Cholesky preconditioning
preconditionerL = ichol(operatorA, struct('michol', 'on'));
preconditionerFunction = @(vectors) ...
    preconditionerL' \ (preconditionerL \ vectors);
[X_preconditioned, lambda_preconditioned, failure_preconditioned] = ...
    lobpcg(randn(n, blockSize), operatorA, [], ...
    preconditionerFunction, exampleSettings.standardTolerance, ...
    exampleSettings.preconditionedMaxIterations, exampleSettings.verbosity);
fprintf('Preconditioned solve convergence flag: %d\n', failure_preconditioned);
disp('Preconditioned eigenvalues:');
disp(lambda_preconditioned.');
fprintf('Preconditioned output class: %s\n', class(X_preconditioned));

%% Example 4: Generalized problem with B = I
[X_identityB, lambda_identityB, failure_identityB] = ...
    lobpcg(randn(n, blockSize), operatorA, speye(n), ...
    preconditionerFunction, exampleSettings.standardTolerance, ...
    exampleSettings.identityBMaxIterations, exampleSettings.verbosity);
fprintf('Generalized B=I convergence flag: %d\n', failure_identityB);
fprintf('Generalized B=I returned %d eigenvectors.\n', size(X_identityB, 2));
disp('Generalized B=I eigenvalues:');
disp(lambda_identityB.');

%% Example 5: Diagonally dominant sparse matrix and preconditioners
diagonalSize = exampleSettings.diagonalSize;
diagonalMatrix = spdiags((1:diagonalSize)', 0, diagonalSize, diagonalSize);
diagonalOperatorA = diagonalMatrix + ...
    sprandsym(diagonalSize, exampleSettings.diagonalDensity);
initialDiagonalVectors = randn(diagonalSize, exampleSettings.diagonalBlockSize);
maxIterations = exampleSettings.diagonalMaxIterations;
residualTolerance = exampleSettings.standardTolerance;
[~, ~, ~, ~, residualsNoPreconditioner] = lobpcg(...
    initialDiagonalVectors, diagonalOperatorA, residualTolerance, maxIterations, 1);
[~, ~, ~, ~, residualsDiagonal] = lobpcg(initialDiagonalVectors, ...
    diagonalOperatorA, [], @(vectors) diagonalMatrix \ vectors, ...
    residualTolerance, maxIterations, 1);

figure;
subplot(2, 2, 1);
plotNoPreconditioner = semilogy(max(residualsNoPreconditioner, [], 1)); hold on;
plotVariant = semilogy(max(residualsDiagonal, [], 1), ':');
legend([plotNoPreconditioner, plotVariant], ...
    {'No preconditioning', 'Diagonal preconditioner'}, 'Location', 'best');
title('Diagonal preconditioner'); axis tight;

% Nonsymmetric preconditioning, as in the documented example.
preconditionerMatrix = diagonalMatrix;
preconditionerMatrix(1, 2) = exampleSettings.nonsymmetricEntry;
[~, ~, ~, ~, residualsNonsymmetric] = lobpcg(initialDiagonalVectors, ...
    diagonalOperatorA, [], @(vectors) preconditionerMatrix \ vectors, ...
    residualTolerance, maxIterations, 1);
subplot(2, 2, 2);
plotNoPreconditioner = semilogy(max(residualsNoPreconditioner, [], 1)); hold on;
plotVariant = semilogy(max(residualsNonsymmetric, [], 1), '--s');
legend([plotNoPreconditioner, plotVariant], ...
    {'No preconditioning', 'Nonsymmetric preconditioner'}, 'Location', 'best');
title('Nonsymmetric preconditioner'); axis tight;

% Nonlinear preconditioning.
preconditionerMatrix = diagonalMatrix;
[~, ~, ~, ~, residualsNonlinear] = lobpcg(initialDiagonalVectors, ...
    diagonalOperatorA, [], ...
    @(vectors) preconditionerMatrix \ (vectors + 10 * sin(vectors)), ...
    residualTolerance, maxIterations, 1);
subplot(2, 2, 3);
plotNoPreconditioner = semilogy(max(residualsNoPreconditioner, [], 1)); hold on;
plotVariant = semilogy(max(residualsNonlinear, [], 1), '-.*');
legend([plotNoPreconditioner, plotVariant], ...
    {'No preconditioning', 'Nonlinear preconditioner'}, 'Location', 'best');
title('Nonlinear preconditioner'); axis tight;

% Selective preconditioning.
preconditionerMatrix = abs(diagonalMatrix - ...
    exampleSettings.selectiveShift * speye(diagonalSize));
[~, ~, ~, ~, residualsSelective] = lobpcg(initialDiagonalVectors, ...
    diagonalOperatorA, [], @(vectors) preconditionerMatrix \ vectors, ...
    residualTolerance, maxIterations, 1);
subplot(2, 2, 4);
plotNoPreconditioner = semilogy(max(residualsNoPreconditioner, [], 1)); hold on;
plotVariant = semilogy(max(residualsSelective, [], 1), '-d');
legend([plotNoPreconditioner, plotVariant], ...
    {'No preconditioning', 'Selective preconditioner'}, 'Location', 'best');
title('Selective preconditioner'); axis tight;

%% Example 6: Codistributed operators
if license('test', 'Distrib_Computing_Toolbox') && ...
        exist('codistributed', 'class') == 8
    existingPool = gcp('nocreate');
    createdPoolForExample = isempty(existingPool);
    if createdPoolForExample
        parallelPool = parpool('local', exampleSettings.parallelWorkerCount);
    else
        parallelPool = existingPool;
    end
    try
        operatorA_distributed = codistributed(diag(1:exampleSettings.parallelSize));
        operatorB_distributed = codistributed(diag(...
            exampleSettings.parallelSize + (1:exampleSettings.parallelSize)));
        [X_distributed, lambda_distributed] = lobpcg(...
            randn(exampleSettings.parallelSize, exampleSettings.parallelBlockSize), ...
            operatorA_distributed, operatorB_distributed, ...
            exampleSettings.standardTolerance, ...
            exampleSettings.parallelMaxIterations, exampleSettings.verbosity);
        disp('Codistributed eigenvalues:');
        disp(gather(lambda_distributed).');
        fprintf('Codistributed output has %d eigenvectors.\n', ...
            size(gather(X_distributed), 2));
    catch exampleError
        if createdPoolForExample && isvalid(parallelPool)
            delete(parallelPool);
        end
        rethrow(exampleError);
    end
    if createdPoolForExample && isvalid(parallelPool)
        delete(parallelPool);
    end
else
    fprintf('Skipping codistributed example: Parallel Computing Toolbox unavailable.\n');
end

%% Example 7: Single-precision inputs
singleOperatorA = single(diag(1:exampleSettings.singleSize));
singleOperatorB = single(diag(...
    exampleSettings.singleSize + (1:exampleSettings.singleSize)));
singleInitialVectors = single(randn(...
    exampleSettings.singleSize, exampleSettings.singleBlockSize));
[X_single, lambda_single] = lobpcg(singleInitialVectors, ...
    singleOperatorA, singleOperatorB, exampleSettings.standardTolerance, ...
    exampleSettings.singleMaxIterations, exampleSettings.verbosity);
fprintf('All-single input/output classes: X=%s, A=%s, B=%s, lambda=%s\n', ...
    class(X_single), class(singleOperatorA), class(singleOperatorB), ...
    class(lambda_single));
fprintf('All-single solve returned %d eigenvectors.\n', size(X_single, 2));

%% Example 8: Sensitivity to clustered eigenvalues and block size
clusteredOperatorA = diag(exampleSettings.clusteredValues);
clusteredBlockSizes = exampleSettings.clusteredBlockSizes;
clusteredResults = cell(size(clusteredBlockSizes));
for blockIndex = 1:numel(clusteredBlockSizes)
    clusteredBlockSize = clusteredBlockSizes(blockIndex);
    [clusteredVectors, clusteredValues, clusteredFailureFlag] = lobpcg(...
        randn(exampleSettings.clusteredMatrixSize, clusteredBlockSize), ...
        clusteredOperatorA, exampleSettings.clusteredTolerance, ...
        exampleSettings.clusteredMaxIterations, exampleSettings.verbosity);
    clusteredOperatorVectors = clusteredOperatorA * clusteredVectors;
    clusteredResiduals = clusteredOperatorVectors - ...
        bsxfun(@times, clusteredVectors, clusteredValues');
    clusteredResidualScale = max(norm(clusteredOperatorVectors, 'fro'), eps);
    clusteredRelativeResidual = norm(clusteredResiduals, 'fro') / ...
        clusteredResidualScale;
    clusteredResults{blockIndex} = struct(...
        'vectors', clusteredVectors, 'values', clusteredValues, ...
        'failureFlag', clusteredFailureFlag, ...
        'relativeResidual', clusteredRelativeResidual);
    assert(size(clusteredVectors, 2) == clusteredBlockSize);
    assert(all(isfinite(clusteredValues)) && isfinite(clusteredRelativeResidual));
    fprintf(['Clustered block size %d: flag %d, relative residual %.3e, ', ...
        '%d eigenpairs returned.\n'], clusteredBlockSize, ...
        clusteredFailureFlag, clusteredRelativeResidual, size(clusteredVectors, 2));
end
    clear rngCleanup
    end
