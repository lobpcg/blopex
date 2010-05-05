function matlab2hypreIJ(A, num_procs, filename, varargin)

% MATLAB2HYPREIJ converts MATLAB matrix to HYPRE IJ format.
%
%    MATLAB2HYPREIJ(A, NUM_PROCS, FILENAME) generates a Hypre formatted
%    matrix from the Matlab sparse matrix A, and saves it in the current
%    directory with filename prefixes given by string FILENAME. The number
%    of processors must be specified by num_procs.
%
%    MATLAB2HYPREIJ(A, NUM_PROCS, FILENAME, PRECISION) specify the 
%    precision of the matrix entries to be written to the Hypre files.
%    Default is '16.15e'.
%
%    % Examples:
%    matlab2hypreIJ(A,16,'c:\hypre\matrixA')
%    matlab2hypreIJ(A,16,'c:\hypre\matrixA','8.7e')
%    % creates 16 files (matrixA.00000, ..., matrixA.00015) to use the 
%    % sparse matrix A in Hypre on a 16 processor system.
%    mpirun -np 2 ./ij -lobpcg -fromfile matrixA -vrand 20
%    % compute the 20 smallest eigenvalues of matrixA using Hypre.
%
%    See also testIJmatlabhypre.m, matlab2hypreParVectors.m, 
%    hypreIJ2matlab.m, hypreParVectors2matlab.m


%    License: BSD
%    Copyright 2010 Bryan C. Smith, Diana Zakaryan, Andrew V. Knyazev
%    $Revision: 2.0 $ $Date: 21-Apr-2010
%    Tested in MATLAB 7.9.0.529 (R2009b) and Octave 3.2.3 and Hypre 2.6.0b.


% Error handling.
if ~issparse(A)
       error('BLOPEX:matlab2hypreIJ:MatrixNotSparse','%s',...
           'Input matrix must be sparse format.')
end

[m, n] = size(A);
if m ~= n
    error('BLOPEX:matlab2hypreIJ:MatrixNotSquare','%s',...
        'Input matrix must be a square matrix.')
end

if ~ischar(filename)
    error('BLOPEX:matlab2hypreIJ:InvalidFilename','%s',...
        'Filename must be a string.')
end


% print format   
if (nargin > 3)
    precision = varargin{1};
else
    precision = '16.15e';
end
prt_format = strcat('%d %d %', precision, '\n');



[hypre_data(:,1), hypre_data(:,2), hypre_data(:,3)] = find(A);
hypre_data = sortrows(hypre_data);
nrows = size(hypre_data,1);
hypre_data(:,1:2) = hypre_data(:,1:2) - 1;

% generate partitioning across processes
% See Hypre getpart.c
part_size = floor(n/num_procs);
rest = mod(n, num_procs);
part = [0 (rest + part_size):part_size:n];

% generate Hypre input files
s0='00000';
index=1;
for i = 1:num_procs
    
    % generate filename suffix and full filename.
    s1 = int2str(i-1);
    ls = length(s1);   
    filename2 = [filename, '.', s0(1:(5-ls)), s1];
    fprintf('Generating file: %s\n', filename2);
    fid = fopen(filename2,'w');

    
    % find indices of rows contained in partition i.
    index1 = index;
    index = index + floor(nrows/num_procs);
    index = min(index, nrows);
    while hypre_data(index,1) >= part(i + 1)
        index = index - 1;
    end
    while (index <= nrows) && (hypre_data(index,1) < part(i + 1))
        index = index + 1;
    end
    index2 = index - 1;
    
    fprintf(fid,'%d %d %d %d\n',part(i),part(i+1)-1,part(i),part(i+1)-1);
    fprintf(fid,prt_format, hypre_data(index1:index2,:)');
    fclose(fid);
end

