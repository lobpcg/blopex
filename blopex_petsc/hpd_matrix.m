function [A,err]=hpd_matrix(n,tau);

err=0;
A=eye(n,n);
a=-1;
b=1;
for i=1:n-1
   for j=i+1:n
      x=a+(b-a)*rand;
      if abs(x)<=tau 
         A(i,j)=x;
         A(j,i)=x;
      end
  end
end

nz=0;
for i=1:n
   for j=1:n
      if A(i,j) ~= 0
         nz=nz+1;
      end
   end
end

sprintf('tau=%d, nz=%d, K(A)=%d',tau,nz,cond(A))

e=eig(A);
for i=1:n
   if e(i)<0 
     'matrix has neg eigenvalues'
     err = 1;
     return
   end
end
 


 