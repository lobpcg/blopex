function [blockVectorX,lambda,varargout] = blopex_matlab(blockVectorX,operatorA,varargin)
%LOBPCG solves Hermitian partial eigenproblems using preconditioning
%
%[blockVectorX,lambda]=lobpcg(blockVectorX,operatorA)
%
%outputs the array of algebraic smallest eigenvalues lambda and corresponding matrix of
%orthonormalized eigenvectors blockVectorX of the Hermitian (full or sparse) operator operatorA  
%using input matrix blockVectorX as an initial guess, without preconditioning, 
%somewhat similar to 
%
%opts.issym=1;opts.isreal=1;K=size(blockVectorX,2);[blockVectorX,lambda]=eigs(operatorA,K,'SA',opts);
%
%for real symmetric operator operatorA, or
%
%K=size(blockVectorX,2);[blockVectorX,lambda]=eigs(operatorA,K,'SR');
%for Hermitian operator operatorA. blockVectorX must be full rank.
%
%[blockVectorX,lambda,failureFlag]=lobpcg(blockVectorX,operatorA) also returns a convergence flag.  
%If failureFlag is 0 then all the eigenvalues converged; otherwise not all converged.
%
%[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory]=...
%lobpcg(blockVectorX,'operatorA','operatorB','operatorT',blockVectorY,...
%residualTolerance,maxIterations,verbosityLevel);
%
%computes smallest eigenvalues lambda and corresponding eigenvectors
%blockVectorX of the generalized eigenproblem Ax=lambda Bx, where 
%Hermitian operators operatorA and operatorB are given as functions, 
%as well as a preconditioner, operatorT. 
%The operators operatorB and operatorT must be in addition POSITIVE DEFINITE. 
%To compute the largest eigenpairs of operatorA, simply apply the code 
%to operatorA multiplied by -1. 
%The code does not involve ANY matrix factorizations of operratorA and operatorB, 
%thus, e.g., it preserves the sparsity and the structure of operatorA and operatorB. 
%
%residualTolerance and maxIterations control tolerance and max number of steps,
%and verbosityLevel = 0, 1, or 2 controls the amount of printed info.
%lambdaHistory is a matrix with all iterative lambdas, and
%residualNormsHistory are matrices of the history of 2-norms of residuals
%
%Required input: 
%   blockVectorX - initial approximation to eigenvectors, full or sparse matrix n-by-blockSize
%   operatorA - the operator of the problem, can be given as a matrix or as an M-file
%
%Optional function input:
%   operatorB - the second operator, if solving a generalized eigenproblem, 
%       can be given as a matrix or as an M-file; by default, or if empty, operatorB=I.
%   operatorT - preconditioner, must be given by an M-file; by default, operatorT=I.
%
%Optional constraints input: 
%   blockVectorY - a full or sparse n-by-sizeY matrix of constraints, sizeY < n. 
%   The iterations will be performed in the (operatorB-) orthogonal 
%   complement of the column-space of blockVectorY. blockVectorY must be full rank. 
%
%Optional scalar input parameters:
%   residualTolerance - tolerance, by default, residualTolerance=n*sqrt(eps)
%   maxIterations - max number of iterations, by default, maxIterations = min(n,20)
%   verbosityLevel - either 0 (no info), 1, or 2 (with pictures); by default, verbosityLevel = 0.
%
%Required output: blockVectorX and lambda are computed blockSize eigenpairs, 
%where blockSize=size(blockVectorX,2) for the initial guess blockVectorX if it is full rank.  
%
%Optional output: failureFlag, lambdaHistory and residualNormsHistory are described above.
%
%Functions operatorA(blockVectorX), operatorB(blockVectorX) and operatorT(blockVectorX) 
%must support blockVectorX being a matrix, not just a column vector.
%
%Every iteration involves one application of operatorA and operatorB, and one of operatorT. 
%
%Main memory requirements: 6 (9 if isempty(operatorB)=0) matrices of the same size
%as blockVectorX, 2 matrices of the same size as blockVectorY (if present), and 
%two square matrices of the size 3*blockSize. 
%
%The following
%Example:
%
%operatorA = 100.*delsq(numgrid('S',21)); [n,n]=size(operatorA);
%[blockVectorX,lambda,failureFlag]=lobpcg(randn(n,10),operatorA,1e-5,50,2);
%
%attempts to compute 10 first eigenpairs of the Poisson operator 
%in a 2x2 square with the mesh size 1/10 without preconditioning,
%but not all eigenpairs converge after 50 steps, so failureFlag=1.  
%
%The next 
%Example:
%
% operatorA = 100.*delsq(numgrid('S',21)); [n,n]=size(operatorA);
% blockVectorY=[];lambda_all=[];
% for j=1:5
%   [blockVectorX,lambda]=lobpcg(randn(n,2),operatorA,blockVectorY,1e-5,200,2);
%   blockVectorY=[blockVectorY,blockVectorX];
%   lambda_all=[lambda_all' lambda']';
% end
%
%attemps to compute the same eigenpairs by calling the code 5 times 
%with blockSize=2 using orthogonalization to the previously founded eigenvectors. 
%
%The following M-script:
%
%global R_cholinc
%operatorA = 100.*delsq(numgrid('S',21)); [n,n]=size(operatorA);
%R_cholinc=cholinc(operatorA,1e-3);
%[blockVectorX,lambda,failureFlag]=lobpcg(randn(n,10),operatorA,[],'precsol',1e-5,50,2);
%
%computes the same eigenpairs in less then 25 steps, so that failureFlag=0, 
%using the preconditioner function precsol:
%
%function blockVectorX=precsol(V)
%global R_cholinc 
%blockVectorX=R_cholinc\(R_cholinc'\V);
%In this example, operatorB=[] must be present in the input parameters. 
%[blockVectorX,lambda,failureFlag]=lobpcg(randn(n,10),operatorA,speye(n),'precsol',1e-5,50,2);
%
%produces similar answers, but is somewhat slower and needs more memory.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This main function LOBPCG is a version of 
%the preconditioned conjugate gradient method (Algorithm 5.1) described in
%A. V. Knyazev, Toward the Optimal Preconditioned Eigensolver:
%Locally Optimal Block Preconditioned Conjugate Gradient Method,
%SIAM Journal on Scientific Computing 23 (2001), no. 2, pp. 517-541. 
%http://epubs.siam.org/sam-bin/getfile/SISC/articles/36612.pdf

%Begin

%Initial settings

failureFlag = 1;
if nargin < 2
    error('There must be at least 2 input agruments: blockVectorX and operatorA')
end
if nargin > 8 
    warning('MATLAB:lobpcg_input_count','There must be at most 8 input agruments unless arguments are passed to a function')
end 

if ischar(blockVectorX)
    error('The first input argument blockVectorX cannot be a string')
end
[n,blockSize]=size(blockVectorX);
if blockSize > n
    error('The first input argument blockVectorX must be tall, not fat')
end
if n < 2
    error('The code does not work for 1x1 matrices')
end

if ~ischar(operatorA)
    [nA,nA] = size(operatorA);
    if any(size(operatorA) ~= nA)
        error('operatorA must be a square matrix or a string.')
    end
    if size(operatorA) ~= n
        error(['The size ' int2str(size(operatorA))...
                ' of operatorA is not the same as ' int2str(n)...
                ' - the number of rows of blockVectorX'])
    end
end

count_string = 0;

operatorT = [];
operatorB = [];
residualTolerance = [];
maxIterations = [];
verbosityLevel = [];
blockVectorY = []; sizeY = 0;
for j = 1:nargin-2
    if isequal(size(varargin{j}),[n,n]) 
        if isempty(operatorB)   
            operatorB = varargin{j};
        else
            error('Too many matrix input arguments. Preconditioner operatorT must be an M-function.')  
        end  
    elseif isequal(size(varargin{j},1),n) && size(varargin{j},2) < n
        if isempty(blockVectorY)   
            blockVectorY = varargin{j};
            sizeY=size(blockVectorY,2);
        else
            error('Something wrong with blockVectorY input argument') 
        end  
    elseif ischar(varargin{j})
        if count_string == 0
            if isempty(operatorB)  
                operatorB = varargin{j};
                count_string = count_string + 1;
            else
                operatorT = varargin{j};
            end
        elseif count_string == 1
            operatorT = varargin{j};
        else
            error('Too many string input arguments')
        end
    elseif isequal(size(varargin{j}),[n,n])
        error('Preconditioner operatorT must be an M-function')
    elseif max(size(varargin{j})) == 1
        if isempty(residualTolerance) 
            residualTolerance = varargin{j};
        elseif isempty(maxIterations)
            maxIterations = varargin{j};       
        elseif isempty(verbosityLevel)
            verbosityLevel = varargin{j};           
        else           
            error('Too many scalar parameters, need only three')
        end
    elseif isempty(varargin{j})
        if isempty(operatorB) 
            count_string = count_string + 1;
        elseif ~isempty(operatorT)
            count_string = count_string + 1;
        elseif ~isempty(blockVectorY)
            error(['Unrecognized empty input argument number ' int2str(j+2)]) %%%%%%% 
        end  
    else
        error(['Input argument number ' int2str(j+2) ' not recognized.'])
    end
end

if verbosityLevel 
    if issparse(blockVectorX)
        fprintf('The sparse initial guess with %i rows and %i columns is detected  \n',n,blockSize) 
    else
        fprintf('The full initial guess with %i rows and %i columns is detected  \n',n,blockSize) 
    end
    if ischar(operatorA)
        fprintf('The main operator is detected as an M-function %s \n',operatorA) 
    elseif issparse(operatorA)
        fprintf('The main operator is detected as a sparse matrix \n') 
    else
        fprintf('The main operator is detected as a full matrix \n')    
    end
    if isempty(operatorB)  
        fprintf('Solving a standard eigenvalue problem, not generalized \n')
    elseif ischar(operatorB) 
        fprintf('The second operator of the generalized eigenproblem \n')
        fprintf('is detected as an M-function %s \n',operatorB) 
    elseif issparse(operatorB)
        fprintf('The second operator of the generalized eigenproblem \n')
        fprintf('is detected as a sparse matrix \n')    
    else
        fprintf('The second operator of the generalized eigenproblem \n')
        fprintf('is detected as a full matrix \n') 
    end
    if isempty(operatorT)  
        fprintf('No preconditioner is detected \n')
    else
        fprintf('The preconditioner is detected as an M-function %s \n',operatorT)
    end    
    if isempty(blockVectorY)  
        fprintf('No matrix of constraints is detected \n')
    elseif issparse(blockVectorY)
        fprintf('The sparse matrix of %i constraints is detected \n',sizeY)
    else
        fprintf('The full matrix of %i constraints is detected \n',sizeY)
    end    
    if issparse(blockVectorY) ~= issparse(blockVectorX)
        warning('MATLAB:lobpcg_sparsity','The sparsity formats of the initial guess and the constraints are inconsistent')
    end
end

% Set defaults

if isempty(residualTolerance)
    residualTolerance = sqrt(eps)*n;   
end
if isempty(maxIterations)
    maxIterations = min(n,20);
end
if isempty(verbosityLevel)
    verbosityLevel = 0;
end

if verbosityLevel 
    fprintf('Solving with tolerance %e and maximum number of iterations %i \n',...
        residualTolerance,maxIterations)
end

% do actual solve by calling C version of blopex
[blockVectorX,lambda,failureFlag,lambdaHistory,residualNormsHistory,...
  iterationNumber] = ...
  blopex_matlab_gateway(blockVectorX,operatorA,...
  operatorB,operatorT,blockVectorY,residualTolerance,maxIterations,...
  verbosityLevel);

varargout(1)={failureFlag};
varargout(2)={lambdaHistory(1:blockSize,1:iterationNumber+1)};
varargout(3)={residualNormsHistory(1:blockSize,1:iterationNumber+1)};
