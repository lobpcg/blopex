function A = hypreIJ2matlab(filename)

% HYPREIJ2MATLAB converts a HYPRE IJ matrix to MATLAB sparse format.
%
%    A = hypreIJ2matlab(filename) converts a HYPRE IJ formatted matrix 
%    in Matlab sparse format. Input argument filename is a string 
%    containing the path and filename prefix of the Hypre files.
%
%    Example:
%    1) Create Hypre formatted matrices 'matrixA.00000', 'matrixA.00001' in 
%       the specified directory.
%    2) A = hypreIJ2matlab ('c:\hypre\matrixA') 
%
%    See also matlabIJ2hypre, matlab2hypreParVectors.m, 
%      testIJmatlabhypre.m, hypreParVectors2matlab.m


%    License: BSD
%    Copyright 2010 Bryan C. Smith, Diana Zakaryan, Andrew V. Knyazev
%    $Revision: 2.0 $ $Date: 21-Apr-2010
%    Tested in MATLAB 7.9.0.529 (R2009b) and Octave 3.2.3 and Hypre 2.6.0b.


% read all the attributes of the specified filename
% in the specified directory
[stat, files_temp] = fileattrib([filename, '*']); 
if ischar(files_temp)
    error('BLOPEX:hypreIJ2matlab:FilesNotFound','%s',...
           'No Hypre IJ matrix files were found with specified path ',...
           'and name.')
end
m = length(files_temp);
files = struct();
files.Name = {};
k = 0;
for i = 1:m
    [pathstr, name, ext] = fileparts(files_temp(i).Name);
    if ~isempty(ext) && sum(double(ext) > 57 | double(ext) < 48)==1
        k = k + 1;
        files(k).Name = files_temp(i).Name;
    end
end
m = k;
if isempty(files(1).Name)
    error('BLOPEX:hypreIJ2matlab:FilesNotFound','%s',...
           'No Hypre IJ matrix files were found with specified path',...
           'and name.')
end


% find size of matrix by looking at last file.
filename = files(m).Name;
hypre_data = dlmread(filename,'',[0 0 0 1]);
imax = hypre_data(2) + 1;
A = spalloc(imax, imax, 1);

% fill the matrix
for i = 1:m
    filename = files(i).Name;
    fprintf('Reading file: %s\n', filename);
    hypre_data = dlmread(filename,'',1,0);
    
    size_hypredata = length(hypre_data);    
    if size_hypredata == 0
        continue;
    end
    A = A + sparse(hypre_data(:,1) + 1, hypre_data(:,2) + 1,...
        hypre_data(:,3), imax, imax, size_hypredata);
end

