%TEST_LOBPCG_INPUT_FORMATS Check LOBPCG numerical consistency by input format.
%
% Call test_lobpcg_input_formats to solve the constrained generalized problem
%
%     A*x = lambda*B*x,  A = diag(1:n),  B = 2*I,  Y = e_1,
%
% for the two smallest eigenpairs. The constraint excludes the first
% coordinate, so the expected eigenvalues are 1 and 1.5. A diagonal
% preconditioner, T = A^-1, is supplied in every solve.
%
% SERIAL FORMAT COMBINATIONS
% For each of double, single, and complex input values, the script runs:
%   1. Full numeric A, B, initial vectors X, and constraints Y.
%   2. Sparse A only.
%   3. Sparse B only.
%   4. Sparse X only.
%   5. Sparse Y only.
%   6. Sparse A, B, X, and Y together.
%   7. Function-handle A and B with a function-handle preconditioner.
%   8. Character-name A, B, and preconditioner functions.
%
% These eight combinations across three value types produce 24 serial test
% cases. Function-handle and character-name helpers support the same
% numeric type as the tested inputs.
%
% CODISTRIBUTED COMBINATIONS
% When Parallel Computing Toolbox is installed and licensed, the script
% starts a two-worker local pool if needed and tests double-precision
% codistributed A only, B only, X only, Y only, and all four inputs
% together. It closes only a pool that it created.
%
% VALIDATION
% Every result is sorted by eigenvalue and its eigenvector phases are
% normalized before comparison with a full double-precision reference.
% Serial cases require failureFlag == 0. Codistributed cases require
% matching eigenvalues and eigenvectors plus a relative generalized
% residual norm no greater than 1e-5, because the convergence flag can be
% conservative for mixed local/codistributed execution.
%
% The function temporarily prioritizes the sibling lobpcg.m on the MATLAB
% path and reports the resolved implementation path. It removes only the
% path entries it adds and leaves caller variables and figures untouched.
function test_lobpcg_input_formats
parameters = struct('matrixSize', 30, 'blockSize', 2, ...
    'constraintCount', 1, 'generalizedScale', 2, ...
    'maxIterations', 100, 'referenceTolerance', 1e-7, ...
    'doubleResidualTolerance', 1e-8, 'singleResidualTolerance', 1e-4, ...
    'singleComparisonTolerance', 5e-3, 'comparisonTolerance', 1e-6, ...
    'relativeResidualToleranceFactor', 10, 'parallelWorkerCount', 2, ...
    'quietVerbosity', 0, 'expectedGeneralizedValues', [1; 1.5]);

previousRngState = rng;
rngCleanup = onCleanup(@() rng(previousRngState));
rng(1729);

scriptFolder = fileparts(mfilename('fullpath'));
pathEntries = strsplit(path, pathsep);
addedScriptFolder = ~any(strcmp(pathEntries, scriptFolder));
if addedScriptFolder
    addpath(scriptFolder, '-begin');
end
pathCleanup = onCleanup(@() remove_path_if_added(scriptFolder, addedScriptFolder));
clear lobpcg
rehash;
expectedLobpcgPath = fullfile(scriptFolder, 'lobpcg.m');
resolvedLobpcgPath = which('lobpcg');
assert(strcmpi(resolvedLobpcgPath, expectedLobpcgPath), ...
    'MATLAB resolved lobpcg to %s instead of %s.', ...
    resolvedLobpcgPath, expectedLobpcgPath);
fprintf('Using LOBPCG implementation: %s\n', resolvedLobpcgPath);

matrixSize = parameters.matrixSize;
blockSize = parameters.blockSize;
initialReal = randn(matrixSize, blockSize);
initialImaginary = randn(matrixSize, blockSize);
helperPath = create_character_operator_helpers(parameters.generalizedScale);
helperCleanup = onCleanup(@() remove_character_operator_helpers(helperPath));

referenceInputs = make_inputs('double', initialReal, initialImaginary, parameters);
reference = solve_lobpcg(referenceInputs, referenceInputs.operatorA, ...
    referenceInputs.operatorB, referenceInputs.preconditioner);
assert_expected_solution(reference, parameters);

precisionNames = {'double', 'single', 'complex'};
sparseFields = {'operatorA', 'operatorB', 'initialVectors', 'constraints'};
resultsPerPrecision = 1 + numel(sparseFields) + 3;
emptyResult = struct('name', '', 'passed', false, 'message', '');
results = repmat(emptyResult, 1, numel(precisionNames) * resultsPerPrecision + 1);
resultIndex = 0;
for precisionIndex = 1:numel(precisionNames)
    precisionName = precisionNames{precisionIndex};
    inputs = make_inputs(precisionName, initialReal, initialImaginary, parameters);
    comparisonTolerance = tolerance_for(precisionName, parameters);

    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' full numeric inputs'], ...
        @() assert_consistent(solve_lobpcg(inputs, inputs.operatorA, ...
        inputs.operatorB, inputs.preconditioner), ...
        reference, comparisonTolerance));

    for fieldIndex = 1:numel(sparseFields)
        fieldName = sparseFields{fieldIndex};
        candidate = sparse_input(inputs, fieldName);
        resultIndex = resultIndex + 1;
        results(resultIndex) = run_case(sprintf('%s sparse %s', ...
            precisionName, fieldName), @() assert_consistent( ...
            solve_lobpcg(candidate, candidate.operatorA, candidate.operatorB, ...
            candidate.preconditioner), reference, comparisonTolerance));
    end

    candidate = sparse_all_numeric_inputs(inputs);
    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' all sparse numeric inputs'], ...
        @() assert_consistent(solve_lobpcg(candidate, candidate.operatorA, ...
        candidate.operatorB, candidate.preconditioner), ...
        reference, comparisonTolerance));
    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' function-handle operators'], ...
        @() assert_consistent(solve_lobpcg(inputs, ...
        @(vectors) inputs.operatorA * vectors, ...
        @(vectors) inputs.operatorB * vectors, inputs.preconditioner), ...
        reference, comparisonTolerance));
    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' character-name operators'], ...
        @() assert_consistent(solve_lobpcg(inputs, ...
        'lobpcg_format_operator_a', 'lobpcg_format_operator_b', ...
        'lobpcg_format_preconditioner'), ...
        reference, comparisonTolerance));
end

resultIndex = resultIndex + 1;
results(resultIndex) = run_case('codistributed matrix and vector inputs', ...
    @() verify_codistributed_inputs(referenceInputs, reference));

passed = sum([results.passed]);
fprintf('\nLOBPCG input-format consistency suite: %d passed, %d failed.\n', ...
    passed, numel(results) - passed);
for result = results(~[results.passed])
    fprintf('  FAIL: %s\n%s\n', result.name, result.message);
end
clear helperCleanup pathCleanup rngCleanup
assert(passed == numel(results), 'LOBPCG input-format consistency suite failed.');
end

function result = run_case(name, testFunction)
    result = struct('name', name, 'passed', false, 'message', '');
    try
        testFunction();
        result.passed = true;
        fprintf('PASS: %s\n', name);
    catch exception
        result.message = getReport(exception, 'extended', 'hyperlinks', 'off');
        fprintf('FAIL: %s\n', name);
    end
end

function inputs = make_inputs(precisionName, initialReal, initialImaginary, parameters)
    matrixSize = size(initialReal, 1);
    diagonal = (1:matrixSize)';
    switch precisionName
        case 'double'
            inputs.operatorA = diag(diagonal);
            inputs.operatorB = parameters.generalizedScale * eye(matrixSize);
            inputs.initialVectors = initialReal;
            inputs.constraints = eye(matrixSize, parameters.constraintCount);
            inputs.diagonal = diagonal;
            inputs.residualTolerance = parameters.doubleResidualTolerance;
        case 'single'
            inputs.operatorA = single(diag(diagonal));
            inputs.operatorB = single(parameters.generalizedScale * eye(matrixSize));
            inputs.initialVectors = single(initialReal);
            inputs.constraints = single(eye(matrixSize, parameters.constraintCount));
            inputs.diagonal = single(diagonal);
            inputs.residualTolerance = single(parameters.singleResidualTolerance);
        case 'complex'
            inputs.operatorA = complex(diag(diagonal));
            inputs.operatorB = complex(parameters.generalizedScale * eye(matrixSize));
            inputs.initialVectors = complex(initialReal, initialImaginary);
            inputs.constraints = complex(eye(matrixSize, parameters.constraintCount));
            inputs.diagonal = complex(diagonal);
            inputs.residualTolerance = parameters.doubleResidualTolerance;
        otherwise
            error('LOBPCG:test:UnknownPrecision', ...
                'Unknown precision %s.', precisionName);
    end
    inputs.preconditioner = @format_preconditioner;
    inputs.maxIterations = parameters.maxIterations;
    inputs.expectedClass = class(inputs.initialVectors);
    inputs.expectsComplex = strcmp(precisionName, 'complex');
    inputs.comparisonTolerance = tolerance_for(precisionName, parameters);
    inputs.relativeResidualToleranceFactor = parameters.relativeResidualToleranceFactor;
    inputs.quietVerbosity = parameters.quietVerbosity;
    inputs.parallelWorkerCount = parameters.parallelWorkerCount;
end

function candidate = sparse_input(inputs, fieldName)
    candidate = inputs;
    candidate.(fieldName) = sparse(candidate.(fieldName));
end

function candidate = sparse_all_numeric_inputs(inputs)
    candidate = inputs;
    candidate.operatorA = sparse(candidate.operatorA);
    candidate.operatorB = sparse(candidate.operatorB);
    candidate.initialVectors = sparse(candidate.initialVectors);
    candidate.constraints = sparse(candidate.constraints);
end

function result = solve_lobpcg(inputs, operatorA, operatorB, operatorT)
    if isnumeric(operatorB) || isa(operatorB, 'codistributed')
        % Numeric B is followed by constraints, then the preconditioner handle.
        [vectors, values, failureFlag] = lobpcg(inputs.initialVectors, ...
            operatorA, operatorB, inputs.constraints, operatorT, ...
            inputs.residualTolerance, inputs.maxIterations, inputs.quietVerbosity);
    else
        % Function/name B and T are parsed before the constraint matrix.
        [vectors, values, failureFlag] = lobpcg(inputs.initialVectors, ...
            operatorA, operatorB, operatorT, inputs.constraints, ...
            inputs.residualTolerance, inputs.maxIterations, inputs.quietVerbosity);
    end
    result = normalize_result(vectors, values, failureFlag, inputs);
end

function result = normalize_result(vectors, values, failureFlag, inputs)
    values = gather(values);
    vectors = gather(vectors);
    expectedVectorSize = size(inputs.initialVectors);
    assert(isequal(size(vectors), expectedVectorSize), ...
        'Expected eigenvector size %s, got %s.', ...
        mat2str(expectedVectorSize), mat2str(size(vectors)));
    assert(numel(values) == expectedVectorSize(2), ...
        'Expected %d eigenvalues, got %d.', expectedVectorSize(2), numel(values));
    assert(isnumeric(vectors) && isnumeric(values));
    assert(strcmp(class(vectors), inputs.expectedClass));
    assert(strcmp(class(values), inputs.expectedClass));
    assert(isreal(vectors) ~= inputs.expectsComplex, ...
        'Eigenvector real/complex type did not match the input format.');
    imaginaryTolerance = 100 * eps(class(values)) * max(1, norm(values, 'fro'));
    assert(max(abs(imag(values))) <= imaginaryTolerance, ...
        'Eigenvalues have non-negligible imaginary parts (max %.3e).', ...
        max(abs(imag(values))));
    values = real(values);
    [values, order] = sort(values);
    vectors = vectors(:, order);
    for vectorIndex = 1:size(vectors, 2)
        [largestEntry, pivot] = max(abs(vectors(:, vectorIndex)));
        assert(largestEntry > 0, 'An eigenvector is unexpectedly zero.');
        if isreal(vectors)
            if vectors(pivot, vectorIndex) < 0
                vectors(:, vectorIndex) = -vectors(:, vectorIndex);
            end
        else
            vectors(:, vectorIndex) = vectors(:, vectorIndex) * ...
                exp(-1i * angle(vectors(pivot, vectorIndex)));
        end
    end
    operatorA = gather(inputs.operatorA);
    operatorB = gather(inputs.operatorB);
    operatorAVectors = operatorA * vectors;
    operatorBVectors = operatorB * vectors;
    % bsxfun preserves compatibility with MATLAB releases before implicit expansion.
    scaledOperatorBVectors = bsxfun(@times, operatorBVectors, values');
    residuals = operatorAVectors - scaledOperatorBVectors;
    residualScale = max([norm(operatorAVectors, 'fro'), ...
        norm(scaledOperatorBVectors, 'fro'), eps(class(values))]);
    relativeResidualNorm = norm(residuals, 'fro') / residualScale;
    result = struct('vectors', double(vectors), 'values', double(values), ...
        'failureFlag', failureFlag, ...
        'relativeResidualNorm', double(relativeResidualNorm), ...
        'relativeResidualToleranceFactor', ...
        inputs.relativeResidualToleranceFactor);
end

function assert_expected_solution(result, parameters)
    assert(result.failureFlag == 0, 'Reference problem did not converge.');
    assert_close(result.values, parameters.expectedGeneralizedValues, ...
        parameters.referenceTolerance);
end

function assert_consistent(result, reference, tolerance)
    assert(result.failureFlag == 0, 'The format-specific problem did not converge.');
    assert_numerically_consistent(result, reference, tolerance);
end

function assert_numerically_consistent(result, reference, tolerance)
    assert(isequal(size(result.values), size(reference.values)), ...
        'Eigenvalue vector sizes differ.');
    assert(isequal(size(result.vectors), size(reference.vectors)), ...
        'Eigenvector matrix sizes differ.');
    assert_relative_close(result.values, reference.values, tolerance);
    residualLimit = reference.relativeResidualToleranceFactor * tolerance;
    assert_relative_close(result.vectors, reference.vectors, residualLimit);
    assert(result.relativeResidualNorm <= residualLimit, ...
        'Relative generalized residual norm %.3e exceeds %.3e.', ...
        result.relativeResidualNorm, residualLimit);
end

function verify_codistributed_inputs(referenceInputs, reference)
    if ~has_parallel_support()
        fprintf('SKIP: Parallel Computing Toolbox is unavailable.\n');
        return
    end

    pool = gcp('nocreate');
    createdPool = isempty(pool);
    if createdPool
        pool = parpool('local', referenceInputs.parallelWorkerCount);
    end
    poolCleanup = onCleanup(@() stop_created_pool(pool, createdPool));

    distributedFields = {'operatorA', 'operatorB', 'initialVectors', 'constraints'};
    for fieldIndex = 1:numel(distributedFields)
        fieldName = distributedFields{fieldIndex};
        candidate = referenceInputs;
        candidate.(fieldName) = codistributed(candidate.(fieldName));
        result = solve_lobpcg(candidate, candidate.operatorA, ...
            candidate.operatorB, candidate.preconditioner);
        print_codistributed_result(fieldName, result);
        assert(result.failureFlag == 0, ...
            'Codistributed %s input failed to converge.', fieldName);
        assert_numerically_consistent(result, reference, ...
            referenceInputs.comparisonTolerance);
    end

    candidate = referenceInputs;
    candidate.operatorA = codistributed(candidate.operatorA);
    candidate.operatorB = codistributed(candidate.operatorB);
    candidate.initialVectors = codistributed(candidate.initialVectors);
    candidate.constraints = codistributed(candidate.constraints);
    result = solve_lobpcg(candidate, candidate.operatorA, ...
        candidate.operatorB, candidate.preconditioner);
    print_codistributed_result('all inputs', result);
    assert(result.failureFlag == 0, ...
        'All-codistributed inputs failed to converge.');
    assert_numerically_consistent(result, reference, ...
        referenceInputs.comparisonTolerance);
    clear poolCleanup
end

function print_codistributed_result(label, result)
    fprintf('Codistributed %s: flag %d, eigenvalues %s, relative residual %.3e\n', ...
        label, result.failureFlag, mat2str(result.values(:).', 8), ...
        result.relativeResidualNorm);
end

function output = format_preconditioner(input)
    diagonal = (1:size(input, 1))';
    if isa(input, 'codistributed')
        diagonal = codistributed(diagonal);
    else
        diagonal = cast(diagonal, 'like', input);
    end
    output = bsxfun(@rdivide, input, diagonal);
end

function available = has_parallel_support
    available = license('test', 'Distrib_Computing_Toolbox') && ...
        exist('codistributed', 'class') == 8;
end

function stop_created_pool(pool, createdPool)
    if createdPool && ~isempty(pool) && isvalid(pool)
        delete(pool);
    end
end

function helperPath = create_character_operator_helpers(generalizedScale)
    helperPath = tempname;
    assert(mkdir(helperPath), 'Unable to create temporary helper directory.');
    write_helper(fullfile(helperPath, 'lobpcg_format_operator_a.m'), sprintf([ ...
        'function output = lobpcg_format_operator_a(input)\n', ...
        'diagonal = cast((1:size(input, 1))'', ''like'', input);\n', ...
        'output = bsxfun(@times, diagonal, input);\n', ...
        'end\n']));
    write_helper(fullfile(helperPath, 'lobpcg_format_operator_b.m'), sprintf([ ...
        'function output = lobpcg_format_operator_b(input)\n', ...
        'output = %.17g * input;\n', ...
        'end\n'], generalizedScale));
    write_helper(fullfile(helperPath, 'lobpcg_format_preconditioner.m'), sprintf([ ...
        'function output = lobpcg_format_preconditioner(input)\n', ...
        'diagonal = cast((1:size(input, 1))'', ''like'', input);\n', ...
        'output = bsxfun(@rdivide, input, diagonal);\n', ...
        'end\n']));
    addpath(helperPath);
end

function write_helper(fileName, source)
    fileId = fopen(fileName, 'w');
    assert(fileId ~= -1, 'Unable to create %s.', fileName);
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '%s', source);
    clear cleanup
end

function remove_character_operator_helpers(helperPath)
    if exist(helperPath, 'dir')
        rmpath(helperPath);
        rmdir(helperPath, 's');
    end
end

function remove_path_if_added(folder, wasAdded)
    if wasAdded
        currentPathEntries = strsplit(path, pathsep);
        if any(strcmp(currentPathEntries, folder))
            rmpath(folder);
        end
    end
end

function tolerance = tolerance_for(precisionName, parameters)
    if strcmp(precisionName, 'single')
        tolerance = parameters.singleComparisonTolerance;
    else
        tolerance = parameters.comparisonTolerance;
    end
end

function assert_close(actual, expected, tolerance)
    difference = norm(actual - expected, 'fro');
    assert(difference <= tolerance, ...
        'Expected difference <= %.3e, got %.3e.', tolerance, difference);
end

function assert_relative_close(actual, expected, relativeTolerance)
    assert(isequal(size(actual), size(expected)), ...
        'Size mismatch: actual is %s; expected is %s.', ...
        mat2str(size(actual)), mat2str(size(expected)));
    referenceScale = max(norm(expected, 'fro'), eps(class(expected)));
    relativeDifference = norm(actual - expected, 'fro') / referenceScale;
    assert(relativeDifference <= relativeTolerance, ...
        'Expected relative difference <= %.3e, got %.3e.', ...
        relativeTolerance, relativeDifference);
end