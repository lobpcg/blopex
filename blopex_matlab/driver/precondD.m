function y = precondD(x)
global operatorA
% may fail when operatorA has a zero on the diagonal
% for the demonstrative purpose only 
n = length(x);
y = spdiags(1./diag(operatorA),0,n,n)*x;
