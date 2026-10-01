%TEST_LOBPCG Exercise the reachable public branches of lobpcg.m.
%
% Run test_lobpcg() from this folder. It uses core MATLAB functionality and
% raises an error when any test fails.
% Tested in MATLAB Version: 26.2.0.3386108 (R2026b)
function test_lobpcg
    C = struct('matrixSize', 30, 'blockSize', 2, ...
        'standardTolerance', 1e-8, 'standardMaxIterations', 80, ...
        'quiet', 0, 'smallMatrixSize', 5, 'fatVectorRows', 6, ...
        'fatVectorColumns', 7, ...
        'smallProblemSize', 20, 'smallProblemBlockSize', 5, ...
        'shortMaxIterations', 1, 'immediateMaxIterations', 5, ...
        'complexTolerance', 1e-6, 'complexMaxIterations', 120, ...
        'singleTolerance', 1e-4, 'singleValueTolerance', 2e-3, ...
        'strictTolerance', 1e-14, 'constraintTolerance', 1e-8, ...
        'warningTolerance', 1e-5, 'residualRankMaxIterations', 5, ...
        'defaultMaxIterations', 20, 'nearExactPerturbation', 1e-12, ...
        'realEigenvalueTolerance', 1e-10, 'valueComparisonTolerance', 1e-6, ...
        'historyTolerance', 1e-12, 'smallestEigenvalues', [1; 2], ...
        'smallestGeneralizedEigenvalues', [0.5; 1], ...
        'constrainedEigenvalues', [2; 3], ...
        'constrainedGeneralizedEigenvalues', [1; 1.5]);
    previousRngState = rng;
    rngCleanup = onCleanup(@() rng(previousRngState));
    rng(42);
    originalPath = path;
    helperPath = create_character_operator_helpers();
    helperCleanup = onCleanup(@() restore_test_environment(helperPath, originalPath));
    
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
    fprintf(['Note: the duplicate matrix-size parser condition has been removed. ', ...
        'The hard-coded SYMMETRIC_CONSTRAINTS branch is not externally ', ...
        'reachable. Distributed-array branches require Parallel Computing ', ...
        'Toolbox; the gather fallback requires pre-R2016a MATLAB.\n']);
    clear helperCleanup rngCleanup

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

function restore_test_environment(helperPath, originalPath)
    path(originalPath);
    if exist(helperPath, 'dir')
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
    [~, operatorA] = makeProblemData();
    expect_error(@() lobpcg('invalid', operatorA), ...
        'BLOPEX:lobpcg:FirstInputNotNumeric');
end

function test_fat_vectors
    expect_error(@() lobpcg(randn(C.fatVectorRows, C.fatVectorColumns), ...
        eye(C.fatVectorRows)), ...
        'BLOPEX:lobpcg:FirstInputFat');
end

function test_small_matrix
    expect_error(@() lobpcg(randn(C.smallMatrixSize, 1), ...
        eye(C.smallMatrixSize)), ...
        'BLOPEX:lobpcg:MatrixTooSmall');
end

function test_small_block_problem
    expect_error(@() lobpcg(randn(C.smallProblemSize, C.smallProblemBlockSize), ...
        diag(1:C.smallProblemSize)), ...
        'BLOPEX:lobpcg:MatrixTooSmall');
end

function test_duplicate_matrix_input
    [initialVectors, operatorA] = makeProblemData();
    expect_error(@() lobpcg(initialVectors, operatorA, ...
        eye(C.matrixSize), eye(C.matrixSize)), ...
        'BLOPEX:lobpcg:TooManyMatrixInputs');
end

function test_duplicate_constraints
    [initialVectors, operatorA] = makeProblemData();
    constraints = eye(C.matrixSize, 1);
    expect_error(@() lobpcg(initialVectors, operatorA, constraints, constraints), ...
        'BLOPEX:lobpcg:WrongConstraintsFormat');
end

function test_unrecognized_input
    [initialVectors, operatorA] = makeProblemData();
    expect_error(@() lobpcg(initialVectors, operatorA, ...
        ones(C.blockSize, C.blockSize + 1)), ...
        'BLOPEX:lobpcg:UnrecognizedInput');
end

function test_unrecognized_empty_input
    [initialVectors, operatorA] = makeProblemData();
    expect_error(@() lobpcg(initialVectors, operatorA, eye(C.matrixSize), ...
        eye(C.matrixSize, 1), []), ...
        'BLOPEX:lobpcg:UnrecognizedEmptyInput');
end

function test_constraints_too_tight
    [initialVectors, operatorA] = makeProblemData();
    expect_error(@() lobpcg(initialVectors(:, 1), operatorA, initialVectors(:, 1)), ...
        'BLOPEX:lobpcg:ConstraintsTooTight');
end

function test_non_positive_definite_b
    [initialVectors, operatorA] = makeProblemData();
    expect_error(@() lobpcg(initialVectors, operatorA, -eye(C.matrixSize)), ...
        'BLOPEX:lobpcg:InitialNotFullRank');
end

function test_too_many_inputs_warning
    [initialVectors, operatorA] = makeProblemData();
    assertWarningId(@() lobpcg(initialVectors, operatorA, [], [], [], [], [], [], []), ...
        'BLOPEX:lobpcg:TooManyInputs');
end

function test_too_many_scalars_warning
    [initialVectors, operatorA] = makeProblemData();
    assertWarningId(@() lobpcg(initialVectors, operatorA, C.warningTolerance, ...
        C.shortMaxIterations, C.quiet, 42), 'BLOPEX:lobpcg:TooManyScalarInputs');
end

function test_too_many_handles_warning
    [initialVectors, operatorA] = makeProblemData();
    assertWarningId(@() lobpcg(initialVectors, operatorA, [], ...
        @diagonal_preconditioner, [], @diagonal_preconditioner, ...
        C.warningTolerance, C.shortMaxIterations, C.quiet), ...
        'BLOPEX:lobpcg:TooManyStringFunctionHandleInputs');
end

function test_default_parameters
    [initialVectors, operatorA] = makeProblemData();
    [eigenvectors, eigenvalues, ~, lambdaHistory, residualHistory] = ...
        lobpcg(initialVectors, operatorA);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert(all(isfinite(eigenvalues)));
    assert(size(lambdaHistory, 2) <= C.defaultMaxIterations);
    assert(size(residualHistory, 2) <= C.defaultMaxIterations - 1);
end

function test_verbosity_reporting
    [initialVectors, operatorA] = makeProblemData();
    constraints = eye(C.matrixSize, 1);
    assert(~isempty(initialVectors) && ~isempty(operatorA) && ~isempty(constraints));

    handleOutput = captureSolverOutput(@() lobpcg(initialVectors, ...
        @test_operator, @test_b_operator, @diagonal_preconditioner, ...
        constraints, C.warningTolerance, C.shortMaxIterations, 1));
    assert(contains(handleOutput, 'main operator is detected as an M-function'));
    assert(contains(handleOutput, 'second operator of the generalized eigenproblem'));
    assert(contains(handleOutput, 'preconditioner is detected as an M-function'));
    assert(contains(handleOutput, 'full matrix of 1 constraints'));

    sparseOutput = captureSolverOutput(@() lobpcg(sparse(initialVectors), ...
        sparse(operatorA), sparse(2 * eye(C.matrixSize)), ...
        sparse(constraints), C.warningTolerance, C.shortMaxIterations, 1));
    assert(contains(sparseOutput, 'sparse initial guess'));
    assert(contains(sparseOutput, 'main operator is detected as a sparse matrix'));
    assert(contains(sparseOutput, 'sparse matrix of 1 constraints'));

    characterOutput = captureSolverOutput(@() lobpcg(initialVectors, ...
        'lobpcg_test_operator', 'lobpcg_test_b_operator', ...
        'lobpcg_test_preconditioner', constraints, C.warningTolerance, ...
        C.shortMaxIterations, 1));
    assert(contains(characterOutput, 'lobpcg_test_operator'));
    assert(contains(characterOutput, 'lobpcg_test_b_operator'));
    assert(contains(characterOutput, 'lobpcg_test_preconditioner'));
end

function test_dense_standard_problem
    [initialVectors, operatorA] = makeProblemData();
    [eigenvectors, eigenvalues, failureFlag] = ...
        solveStandardProblem(initialVectors, operatorA);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
    assert_close(eigenvectors' * eigenvectors, eye(C.blockSize), C.standardTolerance);
end

function test_sparse_standard_problem
    [initialVectors, operatorA] = makeProblemData();
    [eigenvectors, eigenvalues, failureFlag] = ...
        solveStandardProblem(sparse(initialVectors), sparse(operatorA));
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, true, true);
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
end

function test_function_handle_a
    [initialVectors, ~] = makeProblemData();
    [eigenvectors, eigenvalues, failureFlag] = ...
        solveStandardProblem(initialVectors, @test_operator);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
end

function test_character_name_a
    [initialVectors, ~] = makeProblemData();
    [eigenvectors, eigenvalues, failureFlag] = ...
        solveStandardProblem(initialVectors, 'lobpcg_test_operator');
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
end

function test_numeric_b
    [initialVectors, operatorA] = makeProblemData();
    operatorB = 2 * eye(C.matrixSize);
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, operatorB, C.standardTolerance, ...
        C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestGeneralizedEigenvalues, ...
        C.valueComparisonTolerance);
    assert_close(eigenvectors' * operatorB * eigenvectors, ...
        eye(C.blockSize), C.standardTolerance);
end

function test_function_handle_b
    [initialVectors, operatorA] = makeProblemData();
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, @test_b_operator, ...
        C.standardTolerance, C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestGeneralizedEigenvalues, ...
        C.valueComparisonTolerance);
end

function test_character_name_b
    [initialVectors, operatorA] = makeProblemData();
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, 'lobpcg_test_b_operator', ...
        C.standardTolerance, C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestGeneralizedEigenvalues, ...
        C.valueComparisonTolerance);
end

function test_handle_preconditioner
    preconditionerCalls = 0;
    [initialVectors, operatorA] = makeProblemData();
    [~, eigenvalues, failureFlag] = lobpcg(initialVectors, operatorA, [], ...
        @counting_preconditioner, C.standardTolerance, ...
        C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assert(preconditionerCalls > 0);
    assert(numel(eigenvalues) == C.blockSize && isreal(eigenvalues));
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);

    function output = counting_preconditioner(input)
        preconditionerCalls = preconditionerCalls + 1;
        output = diagonal_preconditioner(input);
    end
end

function test_character_preconditioner
    [initialVectors, operatorA] = makeProblemData();
    [~, eigenvalues, failureFlag] = lobpcg(initialVectors, operatorA, [], ...
        'lobpcg_test_preconditioner', C.standardTolerance, ...
        C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assert(numel(eigenvalues) == C.blockSize && isreal(eigenvalues));
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
end

function test_standard_constraints
    [initialVectors, operatorA] = makeProblemData();
    constraints = eye(C.matrixSize, 1);
    [eigenvectors, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, constraints, C.standardTolerance, ...
        C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.constrainedEigenvalues, C.valueComparisonTolerance);
    assert_close(constraints' * eigenvectors, zeros(1, C.blockSize), ...
        C.standardTolerance);
end

function test_generalized_constraints
    [initialVectors, operatorA] = makeProblemData();
    operatorB = 2 * eye(C.matrixSize);
    constraints = eye(C.matrixSize, 1);
    [eigenvectors, eigenvalues, failureFlag] = lobpcg(initialVectors, operatorA, ...
        operatorB, constraints, C.standardTolerance, ...
        C.standardMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.constrainedGeneralizedEigenvalues, ...
        C.valueComparisonTolerance);
    assert_close(constraints' * operatorB * eigenvectors, ...
        zeros(1, C.blockSize), C.standardTolerance);
end

function test_complex_problem
    [initialVectors, operatorA] = makeProblemData();
    operatorA(1, 2) = 1i;
    operatorA(2, 1) = -1i;
    initialVectors = initialVectors + 1i * randn(size(initialVectors));
    [eigenvectors, eigenvalues, failureFlag] = lobpcg(initialVectors, ...
        operatorA, C.complexTolerance, C.complexMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assert(max(abs(imag(eigenvalues))) < C.realEigenvalueTolerance);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, false);
    assert_close(eigenvectors' * eigenvectors, eye(C.blockSize), ...
        C.standardTolerance);
end

function test_single_precision_problem
    [initialVectors, operatorA] = makeProblemData();
    [eigenvectors, eigenvalues] = lobpcg(single(initialVectors), single(operatorA), ...
        single(C.singleTolerance), C.standardMaxIterations, C.quiet);
    assert(size(eigenvectors, 2) == C.blockSize);
    assert(isa(eigenvectors, 'single') && isreal(eigenvectors));
    assert(isa(eigenvalues, 'single'));
    assert(numel(eigenvalues) == C.blockSize);
    assert_close(double(eigenvalues), C.smallestEigenvalues, C.singleValueTolerance);
end

function test_partial_convergence
    [initialVectors, operatorA] = makeProblemData();
    initialVectors(:, 1) = eye(C.matrixSize, 1);
    [~, eigenvalues, failureFlag] = ...
        solveStandardProblem(initialVectors, operatorA);
    assert(failureFlag == 0);
    assert(numel(eigenvalues) == C.blockSize && isreal(eigenvalues));
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
end

function test_immediate_convergence
    [~, operatorA] = makeProblemData();
    initialVectors = eye(C.matrixSize, C.blockSize);
    [eigenvectors, eigenvalues, failureFlag, lambdaHistory, residualHistory] = ...
        lobpcg(initialVectors, operatorA, C.standardTolerance, ...
        C.immediateMaxIterations, C.quiet);
    assert(failureFlag == 0);
    assertEigenpairOutput(eigenvectors, eigenvalues, false, true);
    assert_close(eigenvalues, C.smallestEigenvalues, C.historyTolerance);
    assert(size(lambdaHistory, 2) == 1);
    assert(isempty(residualHistory));
end

function test_incomplete_convergence
    [initialVectors, operatorA] = makeProblemData();
    [~, ~, failureFlag] = lobpcg(initialVectors, operatorA, ...
        C.strictTolerance, C.shortMaxIterations, C.quiet);
    assert(failureFlag == 1);
end

function test_rank_deficient_residual
    [initialVectors, operatorA] = makeProblemData();
    % The preconditioner deliberately collapses all residual columns to rank one.
    assertWarningId(@() lobpcg(initialVectors, operatorA, [], ...
        @rank_deficient_preconditioner, C.standardTolerance, ...
        C.residualRankMaxIterations, C.quiet), ...
        'BLOPEX:lobpcg:ResidualNotFullRank');
end

function test_explicit_gram_path
    [~, operatorA] = makeProblemData();
    initialVectors = eye(C.matrixSize, C.blockSize) + ...
        C.nearExactPerturbation * randn(C.matrixSize, C.blockSize);
    [~, eigenvalues, failureFlag] = ...
        lobpcg(initialVectors, operatorA, C.strictTolerance, ...
        C.shortMaxIterations, C.quiet);
    assert(failureFlag == 1);
    assert_close(eigenvalues, C.smallestEigenvalues, C.valueComparisonTolerance);
end

function test_histories_and_diagnostics
    [initialVectors, operatorA] = makeProblemData();
    existingFigures = findall(0, 'Type', 'figure');
    previousFigureVisibility = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    figureCleanup = onCleanup(@() close_new_figures(existingFigures));
    visibilityCleanup = onCleanup(@() ...
        set(0, 'DefaultFigureVisible', previousFigureVisibility));
    % Verbosity level 2 creates two diagnostics; close only figures created here.
    [~, ~, ~, lambdaHistory, residualHistory] = ...
        solveWithDiagnostics(initialVectors, operatorA);
    currentFigures = findall(0, 'Type', 'figure');
    createdFigures = currentFigures(~ismember(currentFigures, existingFigures));
    assert(size(lambdaHistory, 1) == C.blockSize);
    assert(size(residualHistory, 1) == C.blockSize);
    assert(numel(createdFigures) == 2, ...
        'Expected verbosity diagnostics to create two figures.');
    clear visibilityCleanup figureCleanup
end

function test_sparsity_warning
    [initialVectors, operatorA] = makeProblemData();
    assertWarningId(@() lobpcg(sparse(initialVectors), operatorA, ...
        eye(C.matrixSize, 1), C.standardTolerance, C.shortMaxIterations, 1), ...
        'BLOPEX:lobpcg:SparsityInconsistent');
end

function [eigenvectors, eigenvalues, failureFlag, lambdaHistory, residualHistory] = ...
        solveStandardProblem(initialVectors, operatorA, tolerance, maxIterations)
    if nargin < 3
        tolerance = C.standardTolerance;
        maxIterations = C.standardMaxIterations;
    end
    [eigenvectors, eigenvalues, failureFlag, lambdaHistory, residualHistory] = ...
        lobpcg(initialVectors, operatorA, tolerance, maxIterations, C.quiet);
end

function assertEigenpairOutput(eigenvectors, eigenvalues, expectSparse, expectReal)
    assert(size(eigenvectors, 2) == C.blockSize);
    assert(numel(eigenvalues) == C.blockSize);
    assert(issparse(eigenvectors) == expectSparse);
    assert(isreal(eigenvectors) == expectReal);
    assert(isreal(eigenvalues));
end

function assertWarningId(action, expectedWarningId)
    % Preserve warning configuration and the caller's previous lastwarn state.
    warningState = warning;
    [previousWarningMessage, previousWarningId] = lastwarn;
    cleanup = onCleanup(@() restore_warning_state(warningState, ...
        previousWarningMessage, previousWarningId));
    warning('on', expectedWarningId);
    lastwarn('');
    action();
    [~, actualWarningId] = lastwarn;
    assert(strcmp(actualWarningId, expectedWarningId), ...
        'Expected warning %s, received %s.', expectedWarningId, actualWarningId);
    clear cleanup
end

function restore_warning_state(warningState, warningMessage, warningId)
    warning(warningState);
    lastwarn(warningMessage, warningId);
end

function output = captureSolverOutput(action)
    assert(isa(action, 'function_handle'));
    output = evalc('action();');
end

function [eigenvectors, eigenvalues, failureFlag, lambdaHistory, residualHistory] = ...
        solveWithDiagnostics(initialVectors, operatorA)
    eigenvectors = [];
    eigenvalues = [];
    failureFlag = [];
    lambdaHistory = [];
    residualHistory = [];
    capturedOutput = captureSolverOutput(@run_solver);
    assert(ischar(capturedOutput));

    function run_solver
        [eigenvectors, eigenvalues, failureFlag, lambdaHistory, residualHistory] = ...
            lobpcg(initialVectors, operatorA, C.standardTolerance, ...
            C.standardMaxIterations, 2);
    end
end

function close_new_figures(existingFigures)
    allFigures = findall(0, 'Type', 'figure');
    createdFigures = allFigures(~ismember(allFigures, existingFigures));
    createdFigures = createdFigures(isgraphics(createdFigures));
    if ~isempty(createdFigures)
        close(createdFigures);
    end
end

function [initialVectors, operatorA] = makeProblemData
    initialVectors = randn(C.matrixSize, C.blockSize);
    operatorA = diag(1:C.matrixSize);
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
    assert(isequal(size(actual), size(expected)), ...
        'Size mismatch: actual is %s; expected is %s.', ...
        mat2str(size(actual)), mat2str(size(expected)));
    difference = norm(actual - expected, 'fro');
    assert(difference <= tolerance, ...
        'Expected error <= %.3e, got %.3e.', tolerance, difference);
end
end