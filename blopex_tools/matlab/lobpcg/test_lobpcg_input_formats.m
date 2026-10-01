%TEST_LOBPCG_INPUT_FORMATS Check LOBPCG numerical consistency by input format.
%
% This standalone script solves the constrained generalized problem
%
%     A*x = lambda*B*x,  A = diag(1:30),  B = 2*I,  Y = e_1,
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
% matching eigenvalues and eigenvectors plus a generalized residual norm
% no greater than 1e-5, because the convergence flag can be conservative
% for mixed local/codistributed execution.
%
% The script temporarily prioritizes the sibling lobpcg.m on the MATLAB
% path, reports the resolved implementation path, and restores the
% caller's original path when it finishes.

clear;
close all;
rng(1729);

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
scriptFolder = fileparts(mfilename('fullpath'));
addpath(scriptFolder, '-begin');
clear lobpcg
rehash;
expectedLobpcgPath = fullfile(scriptFolder, 'lobpcg.m');
resolvedLobpcgPath = which('lobpcg');
assert(strcmpi(resolvedLobpcgPath, expectedLobpcgPath), ...
    'MATLAB resolved lobpcg to %s instead of %s.', ...
    resolvedLobpcgPath, expectedLobpcgPath);
fprintf('Using LOBPCG implementation: %s\n', resolvedLobpcgPath);

matrixSize = 30;
blockSize = 2;
initialReal = randn(matrixSize, blockSize);
initialImaginary = randn(matrixSize, blockSize);
helperPath = create_character_operator_helpers();
helperCleanup = onCleanup(@() remove_character_operator_helpers(helperPath));

referenceInputs = make_inputs('double', initialReal, initialImaginary);
reference = solve_numeric_inputs(referenceInputs, true);
assert_expected_solution(reference, 1e-7);

precisionNames = {'double', 'single', 'complex'};
sparseFields = {'operatorA', 'operatorB', 'initialVectors', 'constraints'};
resultsPerPrecision = 1 + numel(sparseFields) + 3;
emptyResult = struct('name', '', 'passed', false, 'message', '');
results = repmat(emptyResult, 1, numel(precisionNames) * resultsPerPrecision + 1);
resultIndex = 0;
for precisionIndex = 1:numel(precisionNames)
    precisionName = precisionNames{precisionIndex};
    inputs = make_inputs(precisionName, initialReal, initialImaginary);
    comparisonTolerance = tolerance_for(precisionName);

    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' full numeric inputs'], ...
        @() assert_consistent(solve_numeric_inputs(inputs, true), ...
        reference, comparisonTolerance));

    for fieldIndex = 1:numel(sparseFields)
        fieldName = sparseFields{fieldIndex};
        candidate = sparse_input(inputs, fieldName);
        resultIndex = resultIndex + 1;
        results(resultIndex) = run_case(sprintf('%s sparse %s', ...
            precisionName, fieldName), @() assert_consistent( ...
            solve_numeric_inputs(candidate, true), reference, comparisonTolerance));
    end

    candidate = sparse_all_numeric_inputs(inputs);
    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' all sparse numeric inputs'], ...
        @() assert_consistent(solve_numeric_inputs(candidate, true), ...
        reference, comparisonTolerance));
    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' function-handle operators'], ...
        @() assert_consistent(solve_function_inputs(inputs), ...
        reference, comparisonTolerance));
    resultIndex = resultIndex + 1;
    results(resultIndex) = run_case([precisionName ' character-name operators'], ...
        @() assert_consistent(solve_character_inputs(inputs), ...
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
assert(passed == numel(results), 'LOBPCG input-format consistency suite failed.');

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

function inputs = make_inputs(precisionName, initialReal, initialImaginary)
    matrixSize = size(initialReal, 1);
    diagonal = (1:matrixSize)';
    switch precisionName
        case 'double'
            inputs.operatorA = diag(diagonal);
            inputs.operatorB = 2 * eye(matrixSize);
            inputs.initialVectors = initialReal;
            inputs.constraints = eye(matrixSize, 1);
            inputs.diagonal = diagonal;
            inputs.residualTolerance = 1e-8;
        case 'single'
            inputs.operatorA = single(diag(diagonal));
            inputs.operatorB = single(2 * eye(matrixSize));
            inputs.initialVectors = single(initialReal);
            inputs.constraints = single(eye(matrixSize, 1));
            inputs.diagonal = single(diagonal);
            inputs.residualTolerance = single(1e-4);
        case 'complex'
            inputs.operatorA = complex(diag(diagonal));
            inputs.operatorB = complex(2 * eye(matrixSize));
            inputs.initialVectors = complex(initialReal, initialImaginary);
            inputs.constraints = complex(eye(matrixSize, 1));
            inputs.diagonal = complex(diagonal);
            inputs.residualTolerance = 1e-8;
        otherwise
            error('LOBPCG:test:UnknownPrecision', ...
                'Unknown precision %s.', precisionName);
    end
    inputs.preconditioner = @format_preconditioner;
    inputs.maxIterations = 100;
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

function result = solve_numeric_inputs(inputs, usePreconditioner)
    if usePreconditioner
        [vectors, values, failureFlag] = lobpcg(inputs.initialVectors, ...
            inputs.operatorA, inputs.operatorB, inputs.constraints, ...
            inputs.preconditioner, inputs.residualTolerance, ...
            inputs.maxIterations, 0);
    else
        [vectors, values, failureFlag] = lobpcg(inputs.initialVectors, ...
            inputs.operatorA, inputs.operatorB, inputs.constraints, ...
            inputs.residualTolerance, inputs.maxIterations, 0);
    end
    result = normalize_result(vectors, values, failureFlag, inputs);
end

function result = solve_function_inputs(inputs)
    operatorA = @(vectors) inputs.operatorA * vectors;
    operatorB = @(vectors) inputs.operatorB * vectors;
    [vectors, values, failureFlag] = lobpcg(inputs.initialVectors, ...
        operatorA, operatorB, inputs.preconditioner, inputs.constraints, ...
        inputs.residualTolerance, inputs.maxIterations, 0);
    result = normalize_result(vectors, values, failureFlag, inputs);
end

function result = solve_character_inputs(inputs)
    [vectors, values, failureFlag] = lobpcg(inputs.initialVectors, ...
        'lobpcg_format_operator_a', 'lobpcg_format_operator_b', ...
        'lobpcg_format_preconditioner', inputs.constraints, ...
        inputs.residualTolerance, inputs.maxIterations, 0);
    result = normalize_result(vectors, values, failureFlag, inputs);
end

function result = normalize_result(vectors, values, failureFlag, inputs)
    values = gather(values);
    vectors = gather(vectors);
    [values, order] = sort(real(values));
    vectors = vectors(:, order);
    for vectorIndex = 1:size(vectors, 2)
        [largestEntry, pivot] = max(abs(vectors(:, vectorIndex)));
        assert(largestEntry > 0, 'An eigenvector is unexpectedly zero.');
        vectors(:, vectorIndex) = vectors(:, vectorIndex) * ...
            exp(-1i * angle(vectors(pivot, vectorIndex)));
    end
    operatorA = gather(inputs.operatorA);
    operatorB = gather(inputs.operatorB);
    residuals = operatorA * vectors - bsxfun(@times, operatorB * vectors, values');
    residualNorm = max(sqrt(sum(conj(residuals) .* residuals)));
    result = struct('vectors', double(vectors), 'values', double(values), ...
        'failureFlag', failureFlag, 'residualNorm', double(residualNorm));
end

function assert_expected_solution(result, tolerance)
    assert(result.failureFlag == 0, 'Reference problem did not converge.');
    assert_close(result.values, [1; 1.5], tolerance);
end

function assert_consistent(result, reference, tolerance)
    assert(result.failureFlag == 0, 'The format-specific problem did not converge.');
    assert_numerically_consistent(result, reference, tolerance);
end

function assert_numerically_consistent(result, reference, tolerance)
    assert_close(result.values, reference.values, tolerance);
    assert_close(result.vectors, reference.vectors, 10 * tolerance);
    assert(result.residualNorm <= 10 * tolerance, ...
        'The generalized residual norm %.3e exceeds %.3e.', ...
        result.residualNorm, 10 * tolerance);
end

function verify_codistributed_inputs(referenceInputs, reference)
    if ~has_parallel_support()
        fprintf('SKIP: Parallel Computing Toolbox is unavailable.\n');
        return
    end

    pool = gcp('nocreate');
    createdPool = isempty(pool);
    if createdPool
        pool = parpool('local', 2);
    end
    poolCleanup = onCleanup(@() stop_created_pool(pool, createdPool));

    distributedFields = {'operatorA', 'operatorB', 'initialVectors', 'constraints'};
    for fieldIndex = 1:numel(distributedFields)
        fieldName = distributedFields{fieldIndex};
        candidate = referenceInputs;
        candidate.(fieldName) = codistributed(candidate.(fieldName));
        result = solve_numeric_inputs(candidate, true);
        print_codistributed_result(fieldName, result);
        assert_numerically_consistent(result, reference, 1e-6);
    end

    candidate = referenceInputs;
    candidate.operatorA = codistributed(candidate.operatorA);
    candidate.operatorB = codistributed(candidate.operatorB);
    candidate.initialVectors = codistributed(candidate.initialVectors);
    candidate.constraints = codistributed(candidate.constraints);
    result = solve_numeric_inputs(candidate, true);
    print_codistributed_result('all inputs', result);
    assert_numerically_consistent(result, reference, 1e-6);
    clear poolCleanup
end

function print_codistributed_result(label, result)
    fprintf(['Codistributed %s: flag %d, eigenvalues [%.8g %.8g], ', ...
        'residual %.3e\n'], label, result.failureFlag, result.values(1), ...
        result.values(2), result.residualNorm);
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

function helperPath = create_character_operator_helpers
    helperPath = tempname;
    assert(mkdir(helperPath), 'Unable to create temporary helper directory.');
    write_helper(fullfile(helperPath, 'lobpcg_format_operator_a.m'), sprintf([ ...
        'function output = lobpcg_format_operator_a(input)\n', ...
        'diagonal = cast((1:size(input, 1))'', ''like'', input);\n', ...
        'output = bsxfun(@times, diagonal, input);\n', ...
        'end\n']));
    write_helper(fullfile(helperPath, 'lobpcg_format_operator_b.m'), sprintf([ ...
        'function output = lobpcg_format_operator_b(input)\n', ...
        'output = 2 * input;\n', ...
        'end\n']));
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

function tolerance = tolerance_for(precisionName)
    if strcmp(precisionName, 'single')
        tolerance = 5e-3;
    else
        tolerance = 1e-6;
    end
end

function assert_close(actual, expected, tolerance)
    difference = norm(actual - expected, 'fro');
    assert(difference <= tolerance, ...
        'Expected difference <= %.3e, got %.3e.', tolerance, difference);
end