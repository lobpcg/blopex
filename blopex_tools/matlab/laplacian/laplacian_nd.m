function [A, lambda, V] = laplacian_nd(varargin)
%LAPLACIAN_ND Sparse negative Laplacian in any number of dimensions.
%
%   A = LAPLACIAN_ND(N) creates a sparse negative Laplacian on a
%   tensor-product grid.  N is a nonempty row vector of positive grid sizes;
%   N(i) specifies the number of points along dimension i.  This form uses
%   Dirichlet conditions ('DD') on every dimension.
%
%   A = LAPLACIAN_ND(N, B) specifies one boundary condition string per grid
%   dimension in the cell row vector B.  Each entry must be one of:
%       'DD'  Dirichlet at both the low and high boundaries.
%       'DN'  Dirichlet low and Neumann high.
%       'ND'  Neumann low and Dirichlet high.
%       'NN'  Neumann at both boundaries.
%       'P'   Periodic at both boundaries.
%
%   [A, LAMBDA] = LAPLACIAN_ND(N, M) returns the M smallest exact
%   eigenvalues with default 'DD' boundaries.  [A, LAMBDA] =
%   LAPLACIAN_ND(N, B, M) uses the supplied boundaries.  In both forms M
%   must be an integer from 0 through PROD(N).  M = 0, the default when B
%   is provided without M, returns A and an empty LAMBDA.
%
%   [A, LAMBDA, V] = LAPLACIAN_ND(N, B, M) additionally returns the
%   orthonormal eigenvectors in V.  A is PROD(N)-by-PROD(N), LAMBDA is
%   M-by-1, and V is PROD(N)-by-M.  LAMBDA is sorted in ascending order and
%   A * V equals V * DIAG(LAMBDA), up to floating-point roundoff.
%
%   The grid point ordering follows MATLAB and Octave Kronecker products:
%   the first dimension in N varies fastest.  A is assembled as a sum of
%   one-dimensional sparse operators, so it remains sparse in all supported
%   dimensions.  Requesting V can consume substantial dense memory when
%   PROD(N) or M is large.
%
%   Examples
%   --------
%   Create the one-dimensional three-point Dirichlet operator and its two
%   smallest eigenpairs.  A is tridiagonal and LAMBDA is
%   [2-SQRT(2); 2].
%
%       [A, lambda, V] = laplacian_nd(3, 2);
%
%   Solve a two-dimensional problem with mixed and periodic boundaries.
%
%       [A, lambda, V] = laplacian_nd([3, 2], {'DN', 'P'}, 2);
%       residual = norm(A * V - V * diag(lambda), 'fro');
%
%   The same loop-based implementation also supports dimensions above three.
%
%       N = [2, 2, 2, 2];
%       B = {'DD', 'NN', 'P', 'DN'};
%       [A, lambda, V] = laplacian_nd(N, B, 2);
%       isOrthonormal = norm(V.' * V - eye(2), 'fro') < 1e-12;
%
%   Noninteger grid sizes are rounded with warning
%   BLOPEX:laplacian_nd:NonIntegerGridSize.  Invalid argument counts,
%   grid sizes, boundary conditions, and eigenpair counts raise errors with
%   BLOPEX:laplacian_nd:* identifiers.  A request for more than three outputs
%   can instead be rejected by MATLAB or Octave before this function begins.
%   A warning is issued when the Mth and (M+1)th eigenvalues are numerically
%   equal.
%
%   Unlike LAPLACIAN, this implementation accepts N and B of any positive
%   length.  Its one-, two-, and three-dimensional results are tested for
%   equality with LAPLACIAN.
%
%   A similar nD Python implementation can be found in
%   https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LaplacianNd.html
%   but is limited to pure Dirichlet, Neumann or Periodic boundary conditions
%   at both ends of each dimension.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   License:  MIT / Apache-2.0
%   Copyright (c) 2026 A.V. Knyazev, Andrew.Knyazev@ucdenver.edu
%   $Revision: 1.0 $  $Date: 1-November-2026
%   Tested in GNU Octave Version: 11.3.0 and
%   MATLAB 26.2 (R2026b)
if nargin > 3
    error('BLOPEX:laplacian_nd:TooManyInputs', ...
        '%s', 'Too many input arguments.');
elseif nargin == 0
    error('BLOPEX:laplacian_nd:NoInputArguments', ...
        '%s', 'Must have at least one input argument.');
end
if nargout > 3
    error('BLOPEX:laplacian_nd:TooManyOutputs', ...
        '%s', 'Maximum number of outputs is 3.');
end

u = validate_grid_sizes(varargin{1});
dimension = numel(u);
defaultBoundaries = repmat({'DD'}, 1, dimension);

if nargin == 3
    B = varargin{2};
    m = varargin{3};
elseif nargin == 2
    secondInput = varargin{2};
    if iscell(secondInput)
        B = secondInput;
        m = 0;
    elseif isnumeric(secondInput) && isscalar(secondInput)
        B = defaultBoundaries;
        m = secondInput;
    else
        error('BLOPEX:laplacian_nd:InvalidSecondInput', ...
            '%s', 'Second input must be a scalar number or a cell array.');
    end
else
    B = defaultBoundaries;
    m = 0;
end

[lowerBoundaries, upperBoundaries] = validate_boundaries(B, dimension);
order = prod(u);
m = validate_eigenvalue_count(m, order);

componentMatrices = cell(1, dimension);
identityMatrices = cell(1, dimension);
componentEigenvalues = cell(1, dimension);
for axisIndex = 1:dimension
    pointCount = u(axisIndex);
    componentMatrices{axisIndex} = component_matrix(pointCount, ...
        lowerBoundaries(axisIndex), upperBoundaries(axisIndex));
    identityMatrices{axisIndex} = speye(pointCount);
    componentEigenvalues{axisIndex} = component_eigenvalues(pointCount, ...
        lowerBoundaries(axisIndex), upperBoundaries(axisIndex));
end

A = sparse(order, order);
for activeAxis = 1:dimension
    term = sparse(1);
    for axisIndex = dimension:-1:1
        if axisIndex == activeAxis
            factor = componentMatrices{axisIndex};
        else
            factor = identityMatrices{axisIndex};
        end
        term = kron(term, factor);
    end
    A = A + term;
end

if m == 0
    lambda = [];
    if nargout == 3
        V = [];
    end
    return
end

[lambda, modes, nextEigenvalue] = combine_eigenvalues(componentEigenvalues, ...
    u, m);
if nargout == 3
    V = combine_eigenvectors(u, lowerBoundaries, upperBoundaries, modes);
end

if abs(lambda(end) - nextEigenvalue) < order * eps('double')
    warning('BLOPEX:laplacian_nd:RepeatedEigenvalue', ...
        'The (M+1)th eigenvalue is nearly equal to the Mth.');
end
end

function u = validate_grid_sizes(value)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
    size(value, 1) ~= 1 || isempty(value) || any(~isfinite(value))
    error('BLOPEX:laplacian_nd:WrongVectorOfGridPoints', ...
        '%s', 'Number of grid points must be a nonempty row vector.');
end
u = double(value);
roundedSizes = round(u);
if any(roundedSizes ~= u)
    warning('BLOPEX:laplacian_nd:NonIntegerGridSize', ...
        '%s', 'Grid sizes must be integers. Rounding...');
    u = roundedSizes;
end
if any(u <= 0)
    error('BLOPEX:laplacian_nd:NonPositiveGridSize', ...
        '%s', 'Grid sizes must be positive.');
end
end

function [lowerBoundaries, upperBoundaries] = validate_boundaries(B, dimension)
if ~iscell(B) || ~isequal(size(B), [1, dimension])
    error('BLOPEX:laplacian_nd:InvalidBdryConds', ...
        '%s', 'Boundary conditions must be a row cell array of length N.');
end

lowerBoundaries = zeros(1, dimension);
upperBoundaries = zeros(1, dimension);
for axisIndex = 1:dimension
    boundary = B{axisIndex};
    if ~ischar(boundary)
        error('BLOPEX:laplacian_nd:InvalidBdryConds', ...
            '%s', 'Boundary conditions must be character vectors.');
    elseif strcmp(boundary, 'P')
        lowerBoundaries(axisIndex) = 3;
        upperBoundaries(axisIndex) = 3;
    elseif numel(boundary) == 2 && ...
            any(boundary(1) == 'DN') && any(boundary(2) == 'DN')
        lowerBoundaries(axisIndex) = boundary_code(boundary(1));
        upperBoundaries(axisIndex) = boundary_code(boundary(2));
    else
        error('BLOPEX:laplacian_nd:InvalidBdryConds', ...
            '%s', ['Boundary conditions must use ''DD'', ''DN'', ''ND'', ', ...
            '''NN'', or ''P''.']);
    end
end
end

function code = boundary_code(boundary)
if boundary == 'D'
    code = 1;
else
    code = 2;
end
end

function m = validate_eigenvalue_count(value, order)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
    ~isfinite(value)
    error('BLOPEX:laplacian_nd:WrongNumberOfEigenvalues', ...
        '%s', 'The requested number of eigenvalues must be a scalar.');
end
m = round(double(value));
if m ~= value || m < 0 || m > order
    error('BLOPEX:laplacian_nd:InvalidNumberOfEigs', ...
        '%s', ['Number of eigenvalues must be a nonnegative integer no ', ...
        'bigger than the number of grid points.']);
end
end

function matrix = component_matrix(pointCount, lowerBoundary, upperBoundary)
edge = ones(pointCount, 1);
matrix = spdiags([-edge, 2 * edge, -edge], [-1, 0, 1], ...
    pointCount, pointCount);

if lowerBoundary == 2
    matrix(1, 1) = 1;
elseif lowerBoundary == 3
    matrix(1, pointCount) = matrix(1, pointCount) - 1;
    matrix(pointCount, 1) = matrix(pointCount, 1) - 1;
end
if upperBoundary == 2
    matrix(pointCount, pointCount) = 1;
end
end

function values = component_eigenvalues(pointCount, lowerBoundary, upperBoundary)
if lowerBoundary == 1 && upperBoundary == 1
    factor = pi / (2 * (pointCount + 1));
    modes = (1:pointCount).';
elseif lowerBoundary == 2 && upperBoundary == 2
    factor = pi / (2 * pointCount);
    modes = (0:pointCount - 1).';
elseif lowerBoundary == 3
    factor = pi / pointCount;
    modes = floor((1:pointCount) / 2).';
else
    factor = pi / (4 * (pointCount + 0.5));
    modes = 2 * (1:pointCount).' - 1;
end
values = 4 * sin(factor * modes) .^ 2;
end

function [lambda, modes, nextEigenvalue] = combine_eigenvalues(componentValues, u, m)
dimension = numel(u);
lambda = componentValues{1};
for axisIndex = 2:dimension
    previousOrder = numel(lambda);
    lambda = kron(ones(u(axisIndex), 1), lambda) + ...
        kron(componentValues{axisIndex}, ones(previousOrder, 1));
end
[lambda, permutation] = sort(lambda);
if m < numel(lambda)
    nextEigenvalue = lambda(m + 1);
else
    nextEigenvalue = inf;
end
lambda = lambda(1:m);
permutation = permutation(1:m).';

modes = zeros(dimension, m);
remaining = permutation - 1;
for axisIndex = 1:dimension
    modes(axisIndex, :) = mod(remaining, u(axisIndex)) + 1;
    remaining = floor(remaining / u(axisIndex));
end
end

function V = combine_eigenvectors(u, lowerBoundaries, upperBoundaries, modes)
dimension = numel(u);
eigenvectorParts = cell(1, dimension);
for axisIndex = 1:dimension
    eigenvectorParts{axisIndex} = component_eigenvectors(u(axisIndex), ...
        lowerBoundaries(axisIndex), upperBoundaries(axisIndex), ...
        modes(axisIndex, :));
end

V = eigenvectorParts{1};
previousOrder = u(1);
for axisIndex = 2:dimension
    V = kron(ones(u(axisIndex), 1), V) .* ...
        kron(eigenvectorParts{axisIndex}, ones(previousOrder, 1));
    previousOrder = previousOrder * u(axisIndex);
end
end

function vectors = component_eigenvectors(pointCount, lowerBoundary, ...
        upperBoundary, modes)
if lowerBoundary == 1 && upperBoundary == 1
    vectors = sin((1:pointCount).' * (pi / (pointCount + 1)) * modes) * ...
        sqrt(2 / (pointCount + 1));
elseif lowerBoundary == 2 && upperBoundary == 2
    vectors = cos((0.5:1:pointCount - 0.5).' * (pi / pointCount) * ...
        (modes - 1)) * sqrt(2 / pointCount);
    vectors(:, modes == 1) = 1 / sqrt(pointCount);
elseif lowerBoundary == 1
    vectors = sin((1:pointCount).' * (pi / (2 * (pointCount + 0.5))) * ...
        (2 * modes - 1)) * sqrt(2 / (pointCount + 0.5));
elseif lowerBoundary == 2
    vectors = cos((0.5:1:pointCount - 0.5).' * ...
        (pi / (2 * (pointCount + 0.5))) * (2 * modes - 1)) * ...
        sqrt(2 / (pointCount + 0.5));
else
    grid = (0.5:1:pointCount - 0.5).';
    vectors = zeros(pointCount, numel(modes));
    oddModes = mod(modes, 2) == 1;
    if any(oddModes)
        vectors(:, oddModes) = cos(grid * ...
            (pi / pointCount * (modes(oddModes) - 1))) * ...
            sqrt(2 / pointCount);
    end
    evenModes = ~oddModes;
    if any(evenModes)
        vectors(:, evenModes) = sin(grid * ...
            (pi / pointCount * modes(evenModes))) * sqrt(2 / pointCount);
    end
    vectors(:, modes == 1) = 1 / sqrt(pointCount);
    if mod(pointCount, 2) == 0
        vectors(:, modes == pointCount) = ...
            vectors(:, modes == pointCount) / sqrt(2);
    end
end
end