function matlab2hypreParVectors(V, num_procs, varargin)

% MATLAB2HYPREPARVECTORS Converts MATLAB Vectors to HYPRE IJ Format.
%
%    MATLAB2HYPREPARVECTORS(V, NUM_PROCS) converts a MATLAB multivector 
%    V to HYPRE format, saving the vectors in the local directory under the
%    filename prefix 'vectors'. The number of processors must be specified 
%    by num_procs.
%
%    MATLAB2HYPREPARVECTORS(V, NUM_PROCS, FILENAME) specify the 
%    path and filename prefix.
%
%    MATLAB2HYPREPARVECTORS(V, NUM_PROCS, FILENAME, PRECISION) specify the 
%    precision of the vectors' entries to be written to the Hypre files. 
%    Default precision is '16.15e'.
%
%
%    Examples:
%    matlab2hypreParVectors(V,16)
%    % Create Hypre files of the vectors that are the columns of V,
%    % partitioning across 16 processors.
%    matlab2hypreParVectors(V,16,'8.7e','c:\hypre\vectors')
%    % specify precision, path, and filename prefix.
%
%    See also matlabIJ2hypre, testIJmatlabhypre.m, 
%      hypreIJ2matlab.m, hypreParVectors2matlab.m


%    License: BSD
%    Copyright 2010 Bryan C. Smith, Diana Zakaryan, Andrew V. Knyazev
%    $Revision: 2.0 $ $Date: 21-Apr-2010
%    Tested in MATLAB 7.9.0.529 (R2009b) and Octave 3.2.3 and Hypre 2.6.0b.

% error handling
if ~isa(V,'numeric')
       error('BLOPEX:matlab2hypreParVectors:InvalidInput','%s',...
           'First input must be a numerical matrix.')
end


if (nargin > 2)
    filename = varargin{1};
else
    filename = 'vectors';
end

if (nargin > 3)
    precision = varargin{2};
else
    precision = '16.15e';
end
prt_format = strcat('%', precision);

[n, m] = size(V);

% Generate partitioning across processors.
% See Hypre getpart.c.
part_size = floor(n/num_procs);
rest = mod(n, num_procs);
part = [0 (rest + part_size):part_size:n];

% generate Hypre input files
Y = [n part(1:num_procs)]';
for i = 1:num_procs
    nrows = part(i+1) - part(i);
    for j = 1:m
        filename2 = [filename, '.', num2str(j-1) , '.', num2str(i-1)];
        fprintf('Generating file: %s\n', filename2);
        X = V((part(i) + 1):part(i+1), j);

        dlmwrite(filename2, nrows, 'precision', '%d');
        dlmwrite(filename2, X, '-append', 'precision', prt_format);

        % writing INFO file
        filename2 = [filename, '.', num2str(j-1), '.', 'INFO', '.',...
            num2str(i-1)];
        fprintf('Generating INFO file: %s\n', filename2);
        dlmwrite(filename2, Y, 'precision', '%d');
    end
end