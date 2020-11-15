function mbernoa
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program computes Bernoulli number Bn using
%                subroutine BERNOA
%       Example: Compute Bernouli number Bn for n = 0,1,...,10
%                Computed results:
%                   n            Bn
%                 --------------------------
%                   0     .100000000000D+01
%                   1    -.500000000000D+00
%                   2     .166666666667D+00
%                   4    -.333333333333D-01
%                   6     .238095238095D-01
%                   8    -.333333333333D-01
%                  10     .757575757576D-01
%       ===========================================================
n=[];b=[];
  b=0;
 b=zeros(1,200+1);
fprintf(1,'%s \n','  please enter nmax');
%        READ(*,*)N
n=10;
[n,b]=bernoa(n,b);
fprintf(1,'%s \n','   n            bn');
fprintf(1,'%s \n',' --------------------------');
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],0,b(0+1));
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],1,b(1+1));
for  k=2:2:n;
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],k,b(k+1));
end;  k=n+1;
%format(2x,i3,d22.12);
end
function [n,bn]=bernoa(n,bn,varargin);
%       ======================================
%       Purpose: Compute Bernoulli number Bn
%       Input :  n --- Serial number
%       Output:  BN(n)--- Bn
%       ======================================
bn(0+1)=1.0d0;
bn(1+1)=-0.5d0;
for  m=2:n;
s=-(1.0d0./(m+1.0d0)-0.5d0);
for  k=2:m-1;
r=1.0d0;
for  j=2:k;
r=r.*(j+m-k)./j;
end;  j=k+1;
s=s-r.*bn(k+1);
end;  k=m-1+1;
bn(m+1)=s;
end;  m=fix(n)+1;
for  m=3:2:n;
bn(m+1)=0.0d0;
end;  m=fix(n)+1;
return;
end

