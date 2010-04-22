function PetscBinaryWrite(filename,precision,varargin)
%
%  Writes in PETSc binary file sparse matrices and vectors
%  if the array is multidimensional and dense it is saved
%  as a one dimensional array
%
%  precision = 'int32' or 'int64' 
%              use 'int64' if Petsc configured --with-64-bit-indices
%
fd = fopen(filename,'w','ieee-be');
for l=1:nargin-2
  A = varargin{l};
  if issparse(A)
    % save sparse matrix in special Matlab format
    [m,n] = size(A);
    majic = 1.2345678910e-30;
    for i=1:min(m,n)
      if A(i,i) == 0
        A(i,i) = majic;
      end
    end
    if min(size(A)) == 1     %a one-rank matrix will be compressed to a
                             %scalar instead of a vectory by sum
      n_nz = full(A' ~= 0);
    else
      n_nz = full(sum(A' ~= 0));
    end
    nz   = sum(n_nz);
    header = fwrite(fd,[1211216,m,n,nz],precision);

    fwrite(fd,n_nz,precision);  %nonzeros per row
    [i,j,s] = find((A.' ~= 0).*(A.'));
    fwrite(fd,i-1,precision);
    for i=1:nz
      if s(i) == majic
        s(i) = 0;
      end
    end

    % modified by Lashuk to handle matrices with compex entries
    if isreal(s)
        if isreal(A)
          fwrite(fd,s,'double');
        else
          ss=nan(2*length(s),1);
          ss(1:2:end)=s;
          ss(2:2:end)=0;
          fwrite(fd,ss,'double');
        end
    else
        ss=nan(2*length(s),1);
        ss(1:2:end)=real(s);
        ss(2:2:end)=imag(s);
        fwrite(fd,ss,'double');
    end
  else
    [m,n] = size(A);
    fwrite(fd,[1211214,m*n],precision);
    if isreal(A)
       fwrite(fd,A,'double');
    else
       ss=nan(2*m,1);
       for i=1:n
          ss(1:2:end)=real(A(:,i));
          ss(2:2:end)=imag(A(:,i));
          fwrite(fd,ss,'double');
       end
    end
  end
end
fclose(fd);
