function test_laplacian_reference()
%TEST_LAPLACIAN_REFERENCE Verify LAPLACIAN against recorded JSON results.
%
%   TEST_LAPLACIAN_REFERENCE loads laplacian_reference.json from this
%   directory, executes each recorded LAPLACIAN call, and compares A,
%   lambda, and V with the pre-recorded values.  It also verifies that A
%   remains sparse and that all output dimensions match their contract.
%
%   The fixture contains 7,075 cases: 1D, 2D, and 3D grid-size vectors
%   composed from 1, 2, and 3; every combination of 'DD', 'DN', 'ND', 'NN',
%   and 'P' boundary conditions; M = 1; and M = 2 where the grid has at
%   least two points.  A is stored densely in JSON because JSON has no
%   sparse-matrix type, but the test separately requires sparse output.
%
%   Values are compared using an absolute tolerance of 1e-12 to allow for
%   floating-point serialization and trigonometric roundoff.  Run this
%   function after adding this directory to the MATLAB or Octave path:
%
%       test_laplacian_reference
%
%   Rebuild the fixture only with GENERATE_LAPLACIAN_REFERENCE when an
%   intentional behavior change requires new reference values.
testDirectory = fileparts(mfilename('fullpath'));
fixturePath = fullfile(testDirectory, 'laplacian_reference.json');
assert(exist(fixturePath, 'file') == 2, ...
    'Missing reference fixture. Run generate_laplacian_reference first.');

fixture = jsondecode(fileread(fixturePath));
assert(fixture.schemaVersion == 1, 'Unsupported fixture schema.');
assert(fixture.caseCount == 7075, 'Fixture must contain 7,075 cases.');
assert(numel(fixture.cases) == fixture.caseCount, ...
    'Fixture case count does not match its contents.');

tolerance = 1e-12;
for caseIndex = 1:fixture.caseCount
    reference = fixture.cases(caseIndex);
    N = double(reference.N(:).');
    B = normalize_boundary_conditions(reference.B);
    M = double(reference.M);
    order = prod(N);

    evalc('[A, lambda, V] = laplacian(N, B, M);');
    expectedA = reshape(double(reference.A), order, order);
    expectedLambda = reshape(double(reference.lambda), M, 1);
    expectedV = reshape(double(reference.V), order, M);
    caseLabel = sprintf('N=%s, B={%s}, M=%d', ...
        mat2str(N), strjoin(B, ', '), M);

    assert(issparse(A), 'A is not sparse for %s.', caseLabel);
    assert(isequal(size(A), [order, order]), ...
        'A has the wrong size for %s.', caseLabel);
    assert(isequal(size(lambda), [M, 1]), ...
        'lambda has the wrong size for %s.', caseLabel);
    assert(isequal(size(V), [order, M]), ...
        'V has the wrong size for %s.', caseLabel);
    assert(max(abs(full(A(:)) - expectedA(:))) <= tolerance, ...
        'A differs from the reference values for %s.', caseLabel);
    assert(max(abs(lambda(:) - expectedLambda(:))) <= tolerance, ...
        'lambda differs from the reference values for %s.', caseLabel);
    assert(max(abs(V(:) - expectedV(:))) <= tolerance, ...
        'V differs from the reference values for %s.', caseLabel);
end

fprintf('Laplacian reference test suite: %d cases passed.\n', fixture.caseCount);
end

function B = normalize_boundary_conditions(value)
%NORMALIZE_BOUNDARY_CONDITIONS Convert decoded JSON boundaries to a row cell.
if iscell(value)
    B = value(:).';
elseif ischar(value)
    B = cellstr(value);
else
    error('BLOPEX:laplacian:InvalidFixtureBoundaryConditions', ...
        'Fixture boundary conditions must be strings.');
end
end