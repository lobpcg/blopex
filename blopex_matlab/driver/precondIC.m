function y = precondIC(x)
global R_cholinc 

y=R_cholinc\(R_cholinc'\x);
