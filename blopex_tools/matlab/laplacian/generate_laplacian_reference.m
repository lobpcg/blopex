function generate_laplacian_reference()
%GENERATE_LAPLACIAN_REFERENCE Create the exhaustive Laplacian JSON fixture.
%
%   GENERATE_LAPLACIAN_REFERENCE evaluates LAPLACIAN for every requested
%   small problem and overwrites laplacian_reference.json beside this file.
%   Each JSON record has these fields:
%       N       Grid-size row vector.
%       B       Cell array of boundary-condition strings.
%       M       Number of requested eigenpairs.
%       A       Full numeric representation of the sparse output matrix.
%       lambda  Exact eigenvalue column vector.
%       V       Orthonormal eigenvector matrix.
%
%   Coverage includes one-, two-, and three-dimensional grids where every
%   dimension is 1, 2, or 3.  Every combination of 'DD', 'DN', 'ND', 'NN',
%   and 'P' boundary conditions is included.  M is 1 for all grids and 2
%   whenever the grid has at least two points, yielding 7,075 records.
%
%   Run this only when intentionally updating the recorded reference
%   values.  It requires JSONENCODE support and write permission for this
%   directory.  Validate a generated fixture with
%   TEST_LAPLACIAN_REFERENCE.
fixturePath = fullfile(fileparts(mfilename('fullpath')), ...
    'laplacian_reference.json');
boundaryValues = {'DD', 'DN', 'ND', 'NN', 'P'};
cases = struct('N', {}, 'B', {}, 'M', {}, 'A', {}, 'lambda', {}, 'V', {});
caseIndex = 0;

for dimension = 1:3
    gridSizes = enumerate_grid_sizes(dimension);
    boundarySets = enumerate_boundary_sets(boundaryValues, dimension);
    for gridIndex = 1:size(gridSizes, 1)
        N = gridSizes(gridIndex, :);
        for boundaryIndex = 1:size(boundarySets, 1)
            B = boundarySets(boundaryIndex, :);
            for M = 1:min(2, prod(N))
                caseIndex = caseIndex + 1;
                evalc('[A, lambda, V] = laplacian(N, B, M);');
                cases(caseIndex) = struct('N', N, 'B', {B}, 'M', M, ...
                    'A', full(A), 'lambda', lambda, 'V', V);
            end
        end
    end
end

assert(caseIndex == 7075, 'Expected 7,075 reference cases.');
fixture = struct('schemaVersion', 1, 'caseCount', caseIndex, 'cases', cases);
fileId = fopen(fixturePath, 'w');
assert(fileId ~= -1, 'Unable to write %s.', fixturePath);
cleanup = onCleanup(@() fclose(fileId));
fwrite(fileId, jsonencode(fixture, 'PrettyPrint', true), 'char');
fwrite(fileId, sprintf('\n'), 'char');
clear cleanup
fprintf('Wrote %d Laplacian reference cases to %s\n', caseIndex, fixturePath);
end

function gridSizes = enumerate_grid_sizes(dimension)
%ENUMERATE_GRID_SIZES Return all DIMENSION-element vectors with values 1:3.
axisSizes = 1:3;
gridSizeCount = numel(axisSizes) ^ dimension;
gridSizes = zeros(gridSizeCount, dimension);
for gridIndex = 0:gridSizeCount - 1
    remainder = gridIndex;
    for axisIndex = dimension:-1:1
        gridSizes(gridIndex + 1, axisIndex) = ...
            axisSizes(mod(remainder, numel(axisSizes)) + 1);
        remainder = floor(remainder / numel(axisSizes));
    end
end
end

function boundarySets = enumerate_boundary_sets(boundaryValues, dimension)
%ENUMERATE_BOUNDARY_SETS Return every boundary-value tuple for DIMENSION.
boundarySetCount = numel(boundaryValues) ^ dimension;
boundarySets = cell(boundarySetCount, dimension);
for boundarySetIndex = 0:boundarySetCount - 1
    remainder = boundarySetIndex;
    for axisIndex = dimension:-1:1
        boundarySets{boundarySetIndex + 1, axisIndex} = ...
            boundaryValues{mod(remainder, numel(boundaryValues)) + 1};
        remainder = floor(remainder / numel(boundaryValues));
    end
end
end