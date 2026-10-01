%TEST_LOBPCG Exercise the reachable public branches of lobpcg.m.
%
% Run this script from the folder containing lobpcg.m. It uses only core
% MATLAB functionality and raises an error when any test fails.

clear;
close all;
rng(42);
helperPath = create_character_operator_helpers();
helperCleanup = onCleanup(@() remove_character_operator_helpers(helperPath));

results = [ ...
    run_case('requires two inputs', @test_missing_inputs), ...
    run_case('requires numeric initial vectors', @test_non_numeric_vectors), ...
    run_case('rejects fat initial vectors', @test_fat_vectors), ...
    run_case('rejects small matrices', @test_small_matrix), ...
    run_case('rejects undersized block problems', @test_small_block_problem), ...
    run_case('rejects duplicate matrix inputs', @test_duplicate_matrix_input), ...
    run_case('rejects duplicate constraints', @test_duplicate_constraints), ...
    run_case('rejects unrecognized inputs', @test_unrecognized_input), ...
    run_case('rejects ambiguous empty inputs', @test_unrecognized_empty_input), ...
    run_case('rejects rank-deficient constrained input', @test_constraints_too_tight), ...
    run_case('rejects non-positive-definite B', @test_non_positive_definite_b), ...
    run_case('warns for excess inputs', @test_too_many_inputs_warning), ...
    run_case('warns for excess scalars', @test_too_many_scalars_warning), ...
    run_case('warns for excess function handles', @test_too_many_handles_warning), ...
    run_case('uses default scalar parameters', @test_default_parameters), ...
    run_case('reports every supported operator form', @test_verbosity_reporting), ...
    run_case('solves dense standard problem', @test_dense_standard_problem), ...
    run_case('solves sparse standard problem', @test_sparse_standard_problem), ...
    run_case('solves function-handle A problem', @test_function_handle_a), ...
    run_case('solves character-name A problem', @test_character_name_a), ...
    run_case('solves generalized numeric B problem', @test_numeric_b), ...
    run_case('solves generalized function-handle B problem', @test_function_handle_b), ...
    run_case('solves generalized character-name B problem', @test_character_name_b), ...
    run_case('applies function-handle preconditioner', @test_handle_preconditioner), ...
    run_case('applies character-name preconditioner', @test_character_preconditioner), ...
    run_case('enforces standard constraints', @test_standard_constraints), ...
    run_case('enforces generalized constraints', @test_generalized_constraints), ...
    run_case('supports complex Hermitian inputs', @test_complex_problem), ...
    run_case('supports single-precision inputs', @test_single_precision_problem), ...
    run_case('freezes converged eigenpairs', @test_partial_convergence), ...
    run_case('detects immediate convergence', @test_immediate_convergence), ...
    run_case('reports incomplete convergence', @test_incomplete_convergence), ...
    run_case('warns for rank-deficient residuals', @test_rank_deficient_residual), ...
    run_case('uses explicit Gram matrices near convergence', @test_explicit_gram_path), ...
    run_case('returns histories and produces diagnostics', @test_histories_and_diagnostics), ...
    run_case('reports sparse-constraint mismatch', @test_sparsity_warning)];

passed = sum([results.passed]);
fprintf('\nLOBPCG test suite: %d passed, %d failed.\n', ...
    passed, numel(results) - passed);
for result = results(~[results.passed])
    fprintf('  FAIL: %s\n%s\n', result.name, result.message);
end

assert(passed == numel(results), 'LOBPCG test suite failed.');
fprintf(['Note: the duplicate matrix-size parser condition has been ', ...
    'removed. The hard-coded SYMMETRIC_CONSTRAINTS branch is not ', ...
    'externally reachable, and distributed-array branches ', ...
    'require the Parallel Computing Toolbox, and the gather fallback ', ...
    'requires a pre-R2016a MATLAB runtime.\n']);

function helperPath = create_character_operator_helpers
    helperPath = tempname;
    assert(mkdir(helperPath), 'Unable to create temporary helper directory.');
    write_helper(fullfile(helperPath, 'lobpcg_test_operator.m'), ...
        sprintf(['function output = lobpcg_test_operator(input)\n', ...
         'output = diag(1:size(input, 1)) * input;\n', ...
         'end\n']));
    write_helper(fullfile(helperPath, 'lobpcg_test_b_operator.m'), ...
        sprintf(['function output = lobpcg_test_b_operator(input)\n', ...
         'output = 2 * input;\n', ...
         'end\n']));
    write_helper(fullfile(helperPath, 'lobpcg_test_preconditioner.m'), ...
        sprintf(['function output = lobpcg_test_preconditioner(input)\n', ...
         'output = diag(1:size(input, 1)) \\ input;\n', ...
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

function result = run_case(name, test_function)
    result = struct('name', name, 'passed', false, 'message', '');
    try
        test_function();
        result.passed = true;
        fprintf('PASS: %s\n', name);
    catch exception
        result.message = getReport(exception, 'extended', 'hyperlinks', 'off');
        fprintf('FAIL: %s\n', name);
    end
end

function test_missing_inputs
    expect_error(@() lobpcg(), 'BLOPEX:lobpcg:NotEnoughInputs');
end

function test_non_numeric_vectors
    [~, operatorA] = problem_data();
    expect_error(@() lobpcg('invalid', operatorA), ...
        'BLOPEX:lobpcg:FirstInputNotNumeric');
end

function test_fat_vectors
    expect_error(@() lobpcg(randn(6, 7), eye(6)), ...
        'BLOPEX:lobpcg:FirstInputFat');
end

function test_small_matrix
    expect_error(@() lobpcg(randn(5, 1), eye(5)), ...
        'BLOPEX:lobpcg:MatrixTooSmall');
end

function test_small_block_problem
    expect_error(@() lobpcg(randn(20, 5), diag(1:20)), ...
        'BLOPEX:lobpcg:MatrixTooSmall');
end

function test_duplicate_matrix_input
    [initialVectors, operatorA] = problem_data();
    expect_error(@() lobpcg(initialVectors, operatorA, eye(30), eye(30)), ...
        'BLOPEX:lobpcg:TooManyMatrixInputs');
end

function test_duplicate_constraints
    [initialVectors, operatorA] = problem_data();
    constraints = eye(30, 1);
    expect_error(@() lobpcg(initialVectors, operatorA, constraints, constraints), ...
        'BLOPEX:lobpcg:WrongConstraintsFormat');
end

function test_unrecognized_input
    [initialVectors, operatorA] = problem_data();
    expect_error(@() lobpcg(initialVectors, operatorA, ones(2, 3)), ...
        'BLOPEX:lobpcg:UnrecognizedInput');
end

function test_unrecognized_empty_input
    [initialVectors, operatorA] = problem_data();
    expect_error(@() lobpcg(initialVectors, operatorA, eye(30), eye(30, 1), []), ...
        'BLOPEX:lobpcg:UnrecognizedEmptyInput');
end

function test_constraints_too_tight
    [initialVectors, operatorA] = problem_data();
    expect_error(@() lobpcg(initialVectors(:, 1), operatorA, initialVectors(:, 1)), ...
        'BLOPEX:lobpcg:ConstraintsTooTight');
end

function test_non_positive_definite_b
    [initialVectors, operatorA] = problem_data();
    expect_error(@() lobpcg(initialVectors, operatorA, -eye(30)), ...
        'BLOPEX:lobpcg:InitialNotFullRank');
end

function test_too_many_inputs_warning
    [initialVectors, operatorA] = problem_data();
    lastwarn('');
    lobpcg(initialVectors, operatorA, [], [], [], [], [], [], []);
    [~, warningId] = lastwarn;
    assert(strcmp(warningId, 'BLOPEX:lobpcg:TooManyInputs'));
end

function test_too_many_scalars_warning
    [initialVectors, operatorA] = problem_data();
    lastwarn('');
    lobpcg(initialVectors, operatorA, 1e-5, 1, 0, 42);
    [~, warningId] = lastwarn;
    assert(strcmp(warningId, 'BLOPEX:lobpcg:TooManyScalarInputs'));
end

function test_too_many_handles_warning
    [initialVectors, operatorA] = problem_data();
    lastwarn('');
    lobpcg(initialVectors, operatorA, [], @diagonal_preconditioner, [], ...
        @diagonal_preconditioner, 1e-5, 1, 0);
    [~, warningId] = lastwarn;
    assert(strcmp(warningId, 'BLOPEX:lobpcg:TooManyStringFunctionHandleInputs'));
end

function test_default_parameters
    [initialVectors, operatorA] = problem_data();
    [~, eigenvalues, ~, lambdaHistory, residualHistory] = ...
        lobpcg(initialVectors, operatorA);
    assert(all(isfinite(eigenvalues)));
    assert(size(lambdaHistory, 2) <= 20);
    assert(size(residualHistory, 2) <= 19);
end

function test_verbosity_reporting
    [initialVectors, operatorA] = problem_data();
    constraints = eye(30, 1);
    assert(~isempty(initialVectors) && ~isempty(operatorA) && ~isempty(constraints));

    handleOutput = evalc(['lobpcg(initialVectors, @test_operator, ', ...
        '@test_b_operator, @diagonal_preconditioner, constraints, 1e-5, 1, 1);']);
    assert(contains(handleOutput, 'main operator is detected as an M-function'));
    assert(contains(handleOutput, 'second operator of the generalized eigenproblem'));
    assert(contains(handleOutput, 'preconditioner is detected as an M-function'));
    assert(contains(handleOutput, 'full matrix of 1 constraints'));

    sparseOutput = evalc(['lobpcg(sparse(initialVectors), sparse(operatorA), ', ...
        'sparse(2 * eye(30)), sparse(constraints), 1e-5, 1, 1);']);
    assert(contains(sparseOutput, 'sparse initial guess'));
    assert(contains(sparseOutput, 'main operator is detected as a sparse matrix'));
    assert(contains(sparseOutput, 'sparse matrix of 1 constraints'));

    characterOutput = evalc(['lobpcg(initialVectors, ''lobpcg_test_operator'', ', ...
        '''lobpcg_test_b_operator'', ''lobpcg_test_preconditioner'', ', ...
        'constraints, 1e-5, 1, 1);']);
    assert(contains(characterOutput, 'lobpcg_test_operator'));
    assert(contains(characterOutput, 'lobpcg_test_b_operator'));
    assert(contains(characterOutput, 'lobpcg_test_preconditioner'));
end

function test_dense_standard_problem
    [initialVectors, operatorA] = problem_data();
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 2], 1e-6);
    assert_close(eigenvectors' * eigenvectors, eye(2), 1e-8);
end

function test_sparse_standard_problem
    [initialVectors, operatorA] = problem_data();
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(sparse(initialVectors), sparse(operatorA), 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert(issparse(eigenvectors));
    assert_close(eigenvalues, [1; 2], 1e-6);
end

function test_function_handle_a
    [initialVectors, ~] = problem_data();
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, @test_operator, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 2], 1e-6);
end

function test_character_name_a
    [initialVectors, ~] = problem_data();
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, 'lobpcg_test_operator', 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 2], 1e-6);
end

function test_numeric_b
    [initialVectors, operatorA] = problem_data();
    operatorB = 2 * eye(30);
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, operatorB, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [0.5; 1], 1e-6);
    assert_close(eigenvectors' * operatorB * eigenvectors, eye(2), 1e-8);
end

function test_function_handle_b
    [initialVectors, operatorA] = problem_data();
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, @test_b_operator, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [0.5; 1], 1e-6);
end

function test_character_name_b
    [initialVectors, operatorA] = problem_data();
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, 'lobpcg_test_b_operator', 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [0.5; 1], 1e-6);
end

function test_handle_preconditioner
    preconditionerCalls = 0;
    [initialVectors, operatorA] = problem_data();
    [~, eigenvalues, failureFlag] = lobpcg(initialVectors, operatorA, [], ...
        @counting_preconditioner, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert(preconditionerCalls > 0);
    assert_close(eigenvalues, [1; 2], 1e-6);

    function output = counting_preconditioner(input)
        preconditionerCalls = preconditionerCalls + 1;
        output = diagonal_preconditioner(input);
    end
end

function test_character_preconditioner
    [initialVectors, operatorA] = problem_data();
    [~, eigenvalues, failureFlag] = lobpcg(initialVectors, operatorA, [], ...
        'lobpcg_test_preconditioner', 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 2], 1e-6);
end

function test_standard_constraints
    [initialVectors, operatorA] = problem_data();
    constraints = eye(30, 1);
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, constraints, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [2; 3], 1e-6);
    assert_close(constraints' * eigenvectors, zeros(1, 2), 1e-8);
end

function test_generalized_constraints
    [initialVectors, operatorA] = problem_data();
    operatorB = 2 * eye(30);
    constraints = eye(30, 1);
    [eigenvectors, eigenvalues, failureFlag] = lobpcg(initialVectors, operatorA, ...
        operatorB, constraints, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 1.5], 1e-6);
    assert_close(constraints' * operatorB * eigenvectors, zeros(1, 2), 1e-8);
end

function test_complex_problem
    [initialVectors, operatorA] = problem_data();
    operatorA(1, 2) = 1i;
    operatorA(2, 1) = -1i;
    initialVectors = initialVectors + 1i * randn(size(initialVectors));
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, 1e-6, 120, 0);
    assert(failureFlag == 0);
    assert(max(abs(imag(eigenvalues))) < 1e-10);
    assert_close(eigenvectors' * eigenvectors, eye(2), 1e-8);
end

function test_single_precision_problem
    [initialVectors, operatorA] = problem_data();
    [~, eigenvalues] = lobpcg(single(initialVectors), single(operatorA), ...
        single(1e-4), 80, 0);
    assert(isa(eigenvalues, 'single'));
    assert_close(double(eigenvalues), [1; 2], 2e-3);
end

function test_partial_convergence
    [initialVectors, operatorA] = problem_data();
    initialVectors(:, 1) = eye(30, 1);
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, 1e-8, 80, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 2], 1e-6);
end

function test_immediate_convergence
    [~, operatorA] = problem_data();
    initialVectors = eye(30, 2);
    [~, eigenvalues, failureFlag, lambdaHistory, residualHistory] = ...
        lobpcg(initialVectors, operatorA, 1e-8, 5, 0);
    assert(failureFlag == 0);
    assert_close(eigenvalues, [1; 2], 1e-12);
    assert(size(lambdaHistory, 2) == 1);
    assert(isempty(residualHistory));
end

function test_incomplete_convergence
    [initialVectors, operatorA] = problem_data();
    [~, ~, failureFlag] = lobpcg(initialVectors, operatorA, 1e-14, 1, 0);
    assert(failureFlag == 1);
end

function test_rank_deficient_residual
    [initialVectors, operatorA] = problem_data();
    lastwarn('');
    lobpcg(initialVectors, operatorA, [], @rank_deficient_preconditioner, 1e-8, 5, 0);
    [~, warningId] = lastwarn;
    assert(strcmp(warningId, 'BLOPEX:lobpcg:ResidualNotFullRank'));
end

function test_explicit_gram_path
    [~, operatorA] = problem_data();
    initialVectors = eye(30, 2) + 1e-12 * randn(30, 2);
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, 1e-14, 1, 0);
    assert(failureFlag == 1);
    assert_close(eigenvalues, [1; 2], 1e-8);
end

function test_histories_and_diagnostics
    [initialVectors, operatorA] = problem_data();
    assert(~isempty(initialVectors) && ~isempty(operatorA));
    previousFigureVisibility = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cleanup = onCleanup(@() set(0, 'DefaultFigureVisible', previousFigureVisibility));
    evalc(['[~, ~, ~, lambdaHistory, residualHistory] = ', ...
        'lobpcg(initialVectors, operatorA, 1e-8, 80, 2);']);
    assert(size(lambdaHistory, 1) == 2);
    assert(size(residualHistory, 1) == 2);
    assert(isgraphics(491) && isgraphics(492));
    close([491, 492]);
    clear cleanup
end

function test_sparsity_warning
    [initialVectors, operatorA] = problem_data();
    lastwarn('');
    lobpcg(sparse(initialVectors), operatorA, eye(30, 1), 1e-8, 1, 1);
    [~, warningId] = lastwarn;
    assert(strcmp(warningId, 'BLOPEX:lobpcg:SparsityInconsistent'));
end

function [initialVectors, operatorA] = problem_data
    initialVectors = randn(30, 2);
    operatorA = diag(1:30);
end

function output = test_operator(input)
    output = diag(1:size(input, 1)) * input;
end

function output = test_b_operator(input)
    output = 2 * input;
end

function output = diagonal_preconditioner(input)
    output = diag(1:size(input, 1)) \ input;
end

function output = rank_deficient_preconditioner(input)
    output = repmat(input(:, 1), 1, size(input, 2));
end

function expect_error(action, expectedIdentifier)
    try
        action();
    catch exception
        assert(strcmp(exception.identifier, expectedIdentifier), ...
            'Expected %s, got %s.', expectedIdentifier, exception.identifier);
        return
    end
    error('LOBPCG:test:ExpectedError', 'Expected error %s was not raised.', ...
        expectedIdentifier);
end

function assert_close(actual, expected, tolerance)
    assert(norm(actual - expected, 'fro') <= tolerance, ...
        'Expected error <= %.3e, got %.3e.', tolerance, norm(actual - expected, 'fro'));
end