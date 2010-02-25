%% Function description.
%  Discretize the following PDE using bilinear quadrilateral FE's 
%  on a uniform mesh:
%  
%        -Laplacian(u)=lambda*u 
%                   
% on a square (0,0)x(s,0)x(s,s)x(0,s). 
% Boundary condition is of Dirichlet type.
%
%  Input:  s - side length of a square
%          n - number of elements along each axis
%  Output: A - sparse stiffness matrix
%          B - sparse mass matrix 
%%
function [A,B] = FEMDiscretizeSq(n,s)

h = 0;          % mesh parameter
n_el = 0;       % total number of elements
N = 0;          % total number of grid nodes

% compute total number of elements
n_el = n^2;

% compute total number of grid nodes
N = (n+1)^2;

% compute mesh parameter h
h = s/n;

% compute system size
m = N - 4*n; 

% Pre-allocate assembling arrays
ID  = (m+1).*ones(1,N); % m+1 will correspond to nodes on the boundary
IEN = zeros(4,n_el);

% Create element stiffness and mass matrices
A_el = (1/6).*[4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4];
B_el = (h^2/9).*[1 1/2 1/4 1/2; 1/2 1 1/2 1/4; 1/4 1/2 1 1/2; 1/2 1/4 1/2 1];

% Set up ID array (number nodes from left to right from bottom to top)
k = 0;
for i = 2:n
    t = (i-1)*(n+1);
    for j = 2:n
        k = k + 1;
        ID(t+j) = k;
    end
end

% Set up IEN array
r = 0;
k = 0;
for i = 1:n_el
    k = k + 1;
    r = rem(k,n+1);
    if r == 0 k = k + 1; end
    IEN(1,i) =k;
    IEN(2,i) =k+1;
    IEN(3,i) =k+n+2;
    IEN(4,i) =k+n+1;
end

% Assembling of A and B
L = (4^2)*n_el;

iA = nan(1,L);
jA = nan(1,L);
sA = zeros(1,L);

iB = nan(1,L);
jB = nan(1,L);
sB = zeros(1,L);

k = 1;
for e= 1:n_el
    gc = [ID(IEN(1,e)) ID(IEN(2,e)) ID(IEN(3,e)) ID(IEN(4,e))]; % global equation numbers
    for i = 1:4
        for j = i:4
            if i == j
                iA(k) = gc(i);
                jA(k) = gc(i);
                sA(k) = A_el(i,i);
                
                iB(k) = gc(i);
                jB(k) = gc(i);
                sB(k) = B_el(i,i);
            else
                iA(k) = gc(i);
                jA(k) = gc(j);
                sA(k) = A_el(i,j);
                
                iB(k) = gc(i);
                jB(k) = gc(j);
                sB(k) = B_el(i,j);
                
                k = k + 1;
                
                % manage the symmetric part
                iA(k) = gc(j);
                jA(k) = gc(i);
                sA(k) = A_el(j,i);
                
                iB(k) = gc(j);
                jB(k) = gc(i);
                sB(k) = B_el(j,i);
                
            end; % end if else
            
            k = k + 1;
            
        end
    end
end

A = sparse(iA,jA,sA,m+1,m+1);
clear iA jA sA
B = sparse(iB,jB,sB,m+1,m+1);
clear iB jB sB
A = A(1:m,1:m);
B = B(1:m,1:m);

clear('ID','IEN');

return;

end
