function accuracyResults = benchmark_lobpcg_attainable_accuracy( ...
    unpreconditionedMaxIterations, preconditionedMaxIterations, ...
    referenceJsonPath, diffJsonPath)
%BENCHMARK_LOBPCG_ATTAINABLE_ACCURACY Measure LOBPCG accuracy near its precision limit.
%
%   RESULTS = BENCHMARK_LOBPCG_ATTAINABLE_ACCURACY() runs 16 solves over
%   two 3-D Dirichlet grids, block sizes 4 and 19, double/single precision,
%   and preconditioned/unpreconditioned modes. Residual tolerance is zero.
%   Iteration caps default to 1000 unpreconditioned and 400 preconditioned.
%   Zero tolerance can cause rank deficiency and ill-conditioning; results
%   may stagnate or lose accuracy near the precision limit.
%
%   RESULTS = BENCHMARK_LOBPCG_ATTAINABLE_ACCURACY(UNPRECONDITIONEDMAXITERATIONS,
%   PRECONDITIONEDMAXITERATIONS, REFERENCEJSONPATH, DIFFJSONPATH) sets the
%   two iteration caps and JSON paths. Existing reference JSON is preserved;
%   a missing baseline is created and the diff file is overwritten each run.
%
%   RESULTS is a 1-by-16 struct array with configuration, iteration/failure
%   metadata, and eigenvalue, subspace, residual, and orthogonality errors.
%   The function prints accuracy metrics only, not eigenvalues, vectors, or
%   timings. Full vectors and histories are not retained.
%
%   Accuracy metrics are computed after sorting the returned eigenvalues and
%   matching eigenvectors in ascending eigenvalue order. Exact eigenvalues
%   come from LAPLACIAN_ND. The exact reference vectors include the complete
%   eigenspace at the block's cutoff eigenvalue, so repeated eigenvalues are
%   compared as a subspace rather than by individual eigenvector directions.
%
%   MAXABSOLUTEEIGENVALUEERROR is max_i(abs(lambda_i - lambdaExact_i)).
%   MAXRELATIVEEIGENVALUEERROR is max_i(abs(lambda_i - lambdaExact_i) /
%   max(abs(lambdaExact_i), eps('double'))).
%   EIGENVECTORSUBSPACEERROR is norm(V - Q*(Q'*V), 'fro') / sqrt(k), where
%   V contains the k computed vectors and Q is an orthonormal basis for the
%   exact reference eigenspace. It measures computed-vector leakage outside
%   that exact subspace; it is not a one-to-one vector difference.
%   RELATIVERESIDUAL is norm(A*V - V*diag(lambda), 'fro') /
%   max(norm(A*V, 'fro'), eps('double')). For single-precision cases, this
%   diagnostic is evaluated in double precision using the original double A.
%   ORTHOGONALITYERROR is norm(V'*V - eye(k), 'fro'). These matrix norms are
%   Frobenius norms. ITERATIONCOUNT is the number of residual-history steps
%   plus one; FAILUREFLAG is the solver's convergence status (0 success,
%   nonzero failure or incomplete convergence).
%   Octave does not support sparse single matrices, so its single-precision
%   cases use double sparse storage and cast operator/preconditioner results
%   to single. The incomplete-Cholesky factorization itself remains double.
%
%   Printed and JSON accuracy deltas are current metric minus the matching
%   reference metric; positive values therefore indicate a larger error.
%
%   MATLAB Command Window: benchmark_lobpcg_attainable_accuracy()
%
%   Windows PowerShell with MATLAB:
%       & 'C:\Program Files\MATLAB\R2026b\bin\matlab.exe' `
%         -sd 'C:\Users\uurun\Documents\GitHub\blopex\blopex_tools\matlab\lobpcg' `
%         -batch 'benchmark_lobpcg_attainable_accuracy'
%
%   Windows PowerShell with Octave:
%       & 'C:\Program Files\GNU Octave\Octave-10.2.0\mingw64\bin\octave-cli.exe' `
%         --quiet `
%         --path 'C:/Users/uurun/Documents/GitHub/blopex/blopex_tools/matlab/lobpcg' `
%         --eval 'benchmark_lobpcg_attainable_accuracy'

if nargin < 1 || isempty(unpreconditionedMaxIterations)
    unpreconditionedMaxIterations = 1000;
end
if nargin < 2 || isempty(preconditionedMaxIterations)
    preconditionedMaxIterations = 400;
end
if nargin < 3
    referenceJsonPath = '';
end
if nargin < 4
    diffJsonPath = '';
end
referenceJsonPath = normalize_json_path(referenceJsonPath);
diffJsonPath = normalize_json_path(diffJsonPath);

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(referenceJsonPath)
    referenceJsonPath = fullfile(scriptFolder, ...
        'benchmark_lobpcg_attainable_accuracy_reference.json');
end
if isempty(diffJsonPath)
    diffJsonPath = make_diff_path(referenceJsonPath);
end
laplacianFolder = fullfile(fileparts(scriptFolder), 'laplacian');
folders = {scriptFolder, laplacianFolder};
oldPathEntries = strsplit(path, pathsep);
addedFolders = {};
for folderIndex = 1:numel(folders)
    if ~any(strcmp(oldPathEntries, folders{folderIndex}))
        addpath(folders{folderIndex}, '-begin');
        addedFolders{end + 1} = folders{folderIndex}; %#ok<AGROW>
    end
end
pathCleanup = onCleanup(@() remove_added_paths(addedFolders));
rehash;
assert(strcmpi(which('lobpcg'), fullfile(scriptFolder, 'lobpcg.m')), ...
    'MATLAB/Octave did not resolve the sibling lobpcg.m.');

oldRng = rng;
rngCleanup = onCleanup(@() rng(oldRng));

gridShapes = [49, 50, 51; 50, 50, 50];
blockSizes = [4, 19];
precisionNames = {'double', 'single'};
solverModes = {'unpreconditioned', 'modified incomplete-Cholesky'};
referenceClusterMargin = 16;
caseCount = size(gridShapes, 1) * numel(blockSizes) * ...
    numel(precisionNames) * numel(solverModes);
emptyCase = struct('precision', '', 'gridShape', [], 'blockSize', 0, ...
    'preconditioner', '', 'maxIterations', 0, 'iterationCount', 0, ...
    'failureFlag', 1, 'residualTolerance', 0, ...
    'maxAbsoluteEigenvalueError', NaN, 'maxRelativeEigenvalueError', NaN, ...
    'eigenvectorSubspaceError', NaN, 'relativeResidual', NaN, ...
    'orthogonalityError', NaN);
accuracyResults = repmat(emptyCase, 1, caseCount);
resultIndex = 0;

for gridIndex = 1:size(gridShapes, 1)
    gridShape = gridShapes(gridIndex, :);
    boundaryConditions = repmat({'DD'}, 1, 3);
    matrixOrder = prod(gridShape);
    referenceModeCount = min(matrixOrder, max(blockSizes) + referenceClusterMargin);
    [operatorADouble, expandedEigenvalues, expandedEigenvectors] = ...
        laplacian_nd(gridShape, boundaryConditions, referenceModeCount);

    for blockIndex = 1:numel(blockSizes)
        blockSize = blockSizes(blockIndex);
        exactEigenvalues = expandedEigenvalues(1:blockSize);
        cutoffTolerance = 100 * eps('double') * max(1, abs(exactEigenvalues(end)));
        referenceSubspaceCount = find(...
            expandedEigenvalues <= exactEigenvalues(end) + cutoffTolerance, 1, 'last');
        if isempty(referenceSubspaceCount)
            error('BLOPEX:benchmark_lobpcg_accuracy:ReferenceSubspaceNotFound', ...
                'No exact reference subspace for grid %s, block %d.', ...
                mat2str(gridShape), blockSize);
        end
        referenceVectors = expandedEigenvectors(:, 1:referenceSubspaceCount);
        [referenceVectors, ~] = qr(referenceVectors, 0);

        rng(1729 + gridIndex * numel(blockSizes) + blockIndex);
        initialVectorsDouble = randn(matrixOrder, blockSize);
        for precisionIndex = 1:numel(precisionNames)
            precisionName = precisionNames{precisionIndex};
            useOctaveSingle = strcmp(precisionName, 'single') && ...
                exist('OCTAVE_VERSION', 'builtin') ~= 0;
            if strcmp(precisionName, 'single')
                initialVectors = single(initialVectorsDouble);
                if useOctaveSingle
                    operatorA = @(vectors) single(operatorADouble * double(vectors));
                    preconditionerMatrix = operatorADouble;
                else
                    operatorA = single(operatorADouble);
                    preconditionerMatrix = operatorA;
                end
            else
                operatorA = operatorADouble;
                initialVectors = initialVectorsDouble;
                preconditionerMatrix = operatorA;
            end

            preconditionerL = ichol(preconditionerMatrix, struct('michol', 'on'));
            if useOctaveSingle
                preconditioner = @(vectors) single(preconditionerL' \ ...
                    (preconditionerL \ double(vectors)));
            else
                preconditioner = @(vectors) preconditionerL' \ ...
                    (preconditionerL \ vectors);
            end
            for modeIndex = 1:numel(solverModes)
                usePreconditioner = modeIndex == 2;
                if usePreconditioner
                    currentPreconditioner = preconditioner;
                    iterationLimit = preconditionedMaxIterations;
                else
                    currentPreconditioner = [];
                    iterationLimit = unpreconditionedMaxIterations;
                end

                metrics = measure_accuracy(operatorA, operatorADouble, ...
                    initialVectors, exactEigenvalues, referenceVectors, ...
                    iterationLimit, currentPreconditioner);
                resultIndex = resultIndex + 1;
                accuracyResults(resultIndex) = struct(...
                    'precision', precisionName, 'gridShape', gridShape, ...
                    'blockSize', blockSize, 'preconditioner', solverModes{modeIndex}, ...
                    'maxIterations', iterationLimit, ...
                    'iterationCount', metrics.iterationCount, ...
                    'failureFlag', metrics.failureFlag, 'residualTolerance', 0, ...
                    'maxAbsoluteEigenvalueError', metrics.maxAbsoluteEigenvalueError, ...
                    'maxRelativeEigenvalueError', metrics.maxRelativeEigenvalueError, ...
                    'eigenvectorSubspaceError', metrics.eigenvectorSubspaceError, ...
                    'relativeResidual', metrics.relativeResidual, ...
                    'orthogonalityError', metrics.orthogonalityError);
                fprintf(['%s grid %s block %d %s: abs-lambda %.3e, ', ...
                    'rel-lambda %.3e, subspace %.3e, residual %.3e, ', ...
                    'orthogonality %.3e\n'], precisionName, mat2str(gridShape), ...
                    blockSize, solverModes{modeIndex}, ...
                    metrics.maxAbsoluteEigenvalueError, ...
                    metrics.maxRelativeEigenvalueError, ...
                    metrics.eigenvectorSubspaceError, metrics.relativeResidual, ...
                    metrics.orthogonalityError);
            end
            clear preconditionerL preconditioner
        end
    end
end

if exist(referenceJsonPath, 'file') == 2
    reference = jsondecode(fileread(referenceJsonPath));
else
    write_accuracy_reference(accuracyResults, referenceJsonPath);
    reference = jsondecode(fileread(referenceJsonPath));
end
accuracyDiff = compare_accuracy_results(accuracyResults, reference.cases);
print_accuracy_diff(accuracyDiff, referenceJsonPath);
write_json_output(accuracyDiff, diffJsonPath);
fprintf('Saved attainable-accuracy diff: %s\n', diffJsonPath);
end

function metrics = measure_accuracy(operatorA, operatorADouble, initialVectors, ...
        exactEigenvalues, exactVectors, maxIterations, preconditioner)
warningState = warning;
warningCleanup = onCleanup(@() warning(warningState));
warning('off', 'all');
if isempty(preconditioner)
    [vectors, values, failureFlag, ~, residualHistory] = ...
        lobpcg(initialVectors, operatorA, 0, maxIterations, 0);
else
    [vectors, values, failureFlag, ~, residualHistory] = ...
        lobpcg(initialVectors, operatorA, [], preconditioner, 0, maxIterations, 0);
end
iterationCount = size(residualHistory, 2) + 1;
if exist('gather', 'file') ~= 0 || exist('gather', 'builtin') ~= 0
    vectors = gather(vectors);
    values = gather(values);
end
vectors = double(vectors);
values = double(values(:));
imaginaryTolerance = 100 * eps('double') * max(1, norm(values, 'fro'));
if max(abs(imag(values))) > imaginaryTolerance
    error('BLOPEX:benchmark_lobpcg_accuracy:ComplexEigenvalues', ...
        'LOBPCG returned non-real eigenvalues for a real symmetric operator.');
end
values = real(values);
[values, order] = sort(values);
vectors = vectors(:, order);
absoluteErrors = abs(values - exactEigenvalues);
relativeErrors = absoluteErrors ./ max(abs(exactEigenvalues), eps('double'));
projection = exactVectors' * vectors;
operatorAVectors = operatorADouble * vectors;
residuals = operatorAVectors - bsxfun(@times, vectors, values');
residualScale = max(norm(operatorAVectors, 'fro'), eps('double'));
metrics = struct('iterationCount', iterationCount, 'failureFlag', failureFlag, ...
    'maxAbsoluteEigenvalueError', max(absoluteErrors), ...
    'maxRelativeEigenvalueError', max(relativeErrors), ...
    'eigenvectorSubspaceError', ...
        norm(vectors - exactVectors * projection, 'fro') / sqrt(size(vectors, 2)), ...
    'relativeResidual', norm(residuals, 'fro') / residualScale, ...
    'orthogonalityError', norm(vectors' * vectors - eye(size(vectors, 2)), 'fro'));
end

function write_accuracy_reference(results, filePath)
reference = struct('schemaVersion', 1, ...
    'benchmarkName', 'LOBPCG attainable-accuracy comparison', ...
    'generatedAt', accuracy_timestamp(), 'matlabVersion', version, ...
    'residualTolerance', 0, 'cases', {results});
write_json_output(reference, filePath);
fprintf('Saved attainable-accuracy reference: %s\n', filePath);
end

function diffReport = compare_accuracy_results(currentResults, referenceCases)
comparisonTemplate = struct('precision', '', 'gridShape', [], 'blockSize', 0, ...
    'preconditioner', '', 'referenceFound', false, ...
    'maxAbsoluteEigenvalueErrorDelta', NaN, 'maxRelativeEigenvalueErrorDelta', NaN, ...
    'eigenvectorSubspaceErrorDelta', NaN, 'relativeResidualDelta', NaN, ...
    'orthogonalityErrorDelta', NaN);
comparisons = repmat(comparisonTemplate, 1, numel(currentResults));
referenceMatched = false(1, numel(referenceCases));
for currentIndex = 1:numel(currentResults)
    current = currentResults(currentIndex);
    matchIndex = [];
    for referenceIndex = 1:numel(referenceCases)
        reference = referenceCases(referenceIndex);
        if strcmp(current.precision, reference.precision) && ...
                isequal(current.gridShape, reshape(double(reference.gridShape), 1, [])) && ...
                current.blockSize == reference.blockSize && ...
                strcmp(current.preconditioner, reference.preconditioner)
            matchIndex = referenceIndex;
            break
        end
    end
    comparisons(currentIndex).precision = current.precision;
    comparisons(currentIndex).gridShape = current.gridShape;
    comparisons(currentIndex).blockSize = current.blockSize;
    comparisons(currentIndex).preconditioner = current.preconditioner;
    comparisons(currentIndex).referenceFound = ~isempty(matchIndex);
    if ~isempty(matchIndex)
        referenceMatched(matchIndex) = true;
        reference = referenceCases(matchIndex);
        comparisons(currentIndex).maxAbsoluteEigenvalueErrorDelta = ...
            current.maxAbsoluteEigenvalueError - reference.maxAbsoluteEigenvalueError;
        comparisons(currentIndex).maxRelativeEigenvalueErrorDelta = ...
            current.maxRelativeEigenvalueError - reference.maxRelativeEigenvalueError;
        comparisons(currentIndex).eigenvectorSubspaceErrorDelta = ...
            current.eigenvectorSubspaceError - reference.eigenvectorSubspaceError;
        comparisons(currentIndex).relativeResidualDelta = ...
            current.relativeResidual - reference.relativeResidual;
        comparisons(currentIndex).orthogonalityErrorDelta = ...
            current.orthogonalityError - reference.orthogonalityError;
    end
end
unmatched = find(~referenceMatched);
referenceOnly = repmat(struct('precision', '', 'gridShape', [], ...
    'blockSize', 0, 'preconditioner', ''), 1, numel(unmatched));
for index = 1:numel(unmatched)
    item = referenceCases(unmatched(index));
    referenceOnly(index) = struct('precision', item.precision, ...
        'gridShape', reshape(double(item.gridShape), 1, []), ...
        'blockSize', item.blockSize, 'preconditioner', item.preconditioner);
end
diffReport = struct('schemaVersion', 1, 'generatedAt', accuracy_timestamp(), ...
    'matlabVersion', version, 'comparisons', comparisons, ...
    'referenceOnly', referenceOnly);
end

function print_accuracy_diff(report, referencePath)
fprintf('\nAccuracy deltas versus %s\n', referencePath);
for index = 1:numel(report.comparisons)
    item = report.comparisons(index);
    fprintf('%s grid %s block %d %s: ', item.precision, ...
        mat2str(item.gridShape), item.blockSize, item.preconditioner);
    if item.referenceFound
        fprintf('abs-lambda %+.2e, rel-lambda %+.2e, subspace %+.2e, ', ...
            item.maxAbsoluteEigenvalueErrorDelta, ...
            item.maxRelativeEigenvalueErrorDelta, ...
            item.eigenvectorSubspaceErrorDelta);
        fprintf('residual %+.2e, orthogonality %+.2e\n', ...
            item.relativeResidualDelta, item.orthogonalityErrorDelta);
    else
        fprintf('no matching reference\n');
    end
end
end

function write_json_output(value, filePath)
try
    jsonText = jsonencode(value, 'PrettyPrint', true);
catch
    jsonText = jsonencode(value);
end
fileId = fopen(filePath, 'w');
assert(fileId ~= -1, 'Could not open JSON file %s.', filePath);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '%s\n', jsonText);
end

function pathValue = normalize_json_path(pathValue)
if isempty(pathValue)
    pathValue = '';
elseif isstring(pathValue) && isscalar(pathValue)
    pathValue = char(pathValue);
elseif ~ischar(pathValue)
    error('BLOPEX:benchmark_lobpcg_accuracy:InvalidJsonPath', ...
        'JSON paths must be character vectors or scalar strings.');
end
end

function filePath = make_diff_path(referencePath)
[folder, name, extension] = fileparts(referencePath);
if isempty(extension)
    extension = '.json';
end
filePath = fullfile(folder, [name, '_diff', extension]);
end

function timestamp = accuracy_timestamp()
if exist('OCTAVE_VERSION', 'builtin')
    timestamp = datestr(now, 30);
else
    timestamp = char(datetime("now", 'Format', "yyyyMMdd'T'HHmmss"));
end
end

function remove_added_paths(folders)
for folderIndex = 1:numel(folders)
    currentEntries = strsplit(path, pathsep);
    if any(strcmp(currentEntries, folders{folderIndex}))
        rmpath(folders{folderIndex});
    end
end
end
