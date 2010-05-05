function V = hypreParVectors2matlab(filename,varargin)
% HYPREPARVECTORS2MATLAB Converts Hypre IJ Vectors to Matlab.
%
%    V = HYPREPARVECTORS2MATLAB() converts Hypre IJ vectors with filename 
%    prefix 'vectors' stored in the local directory. Determines the length
%    of the vectors from the file vectors.0.INFO.0.
%
%    V = HYPREPARVECTORS2MATLAB(FILENAME) converts Hypre IJ vectors with 
%    filename prefix from input argument filename.
%
%    V = HYPREPARVECTORS2MATLAB(FILENAME, N) if info file filename.0.INFO.0
%    is not available, one can specify the length of vectors with input
%    argument n.
%
%    Example:
%    1) Create Hypre formatted vectors 'vectors.0.0', etc. in the specified
%       directory.
%    2) A = hypreIJ2matlab('c:\hypre\matrixA') 
%
%    See also matlabIJ2hypre, matlab2hypreParVectors.m, 
%      hypreIJ2matlab.m, testIJmatlabhypre.m


%    License: BSD
%    Copyright 2010 Bryan C. Smith, Diana Zakaryan, Andrew V. Knyazev
%    $Revision: 2.0 $ $Date: 21-Apr-2010
%    Tested in MATLAB 7.9.0.529 (R2009b) and Octave 3.2.3 and Hypre 2.6.0b.


[stat, files_temp] = fileattrib([filename, '*']); 
if ischar(files_temp)
    error('BLOPEX:hypreIJ2matlab:FilesNotFound','%s',...
           'No Hypre IJ matrix files were found with specified path ',...
           'and name.')
end
m = length(files_temp);
files = struct();
files.Name = {};
num_procs = 0;
k = 0;
for i = 1:m
    [pathstr, name, ext] = fileparts(files_temp(i).Name);
    kk = findstr('.INFO',name);
    if length(kk) == 1
        continue;
    elseif ~isempty(ext) && sum(double(ext) > 57 | double(ext) < 48) == 1
        k = k + 1;
        files(k).Name = files_temp(i).Name;
        a = length(ext);
        num_procs = max(str2double(ext(2:a)) + 1, num_procs);
    end
end
if isempty(files(1).Name)
    error('BLOPEX:hypreParVectors2matlab:FilesNotFound','%s',...
           'No Hypre IJ matrix files were found with specified path',...
           'and name.')
end

if nargin < 2
    n = dlmread([filename, '.0.INFO.0', '', [0,0,0,0]]);
else
    n = varargin(1);
end

% calls the actual program num_of_vectors time
% in order to convert all num_of_vectors hypre vectors to matlab format
length(files)
num_procs
num_vectors = length(files)/num_procs;
if num_vectors ~= floor(num_vectors)
        error('BLOPEX:hypreParVectors2matlab:MissingFiles','%s',...
           'Some files are missing.')
end

V = zeros(n,num_vectors);
for j = 1:num_vectors
    k = 1;
    for i = 1:num_procs
        filename2 = [filename, '.', num2str(j-1), '.', num2str(i-1)];
        fprintf('Reading file: %s\n', filename2);
        hypre_data = dlmread(filename2,';',1,0);
        m = length(hypre_data);
        V(k:(k+m-1),j) = hypre_data;
        k = k + m;
    end
end