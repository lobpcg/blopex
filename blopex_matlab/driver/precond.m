function y = precond(x)

n = length(x);
invA = spdiags(1./[1:n]',0,n,n);
y = invA*x;