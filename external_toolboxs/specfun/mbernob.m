function mbernob
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
%                subroutine BERNOB
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
[n,b]=bernob(n,b);
fprintf(1,'%s \n','   n            bn');
fprintf(1,'%s \n',' --------------------------');
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],0,b(0+1));
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],1,b(1+1));
for  k=2:2:n;
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],k,b(k+1));
end;  k=n+1;
%format(2x,i3,d22.12);
end
function [n,bn]=bernob(n,bn,varargin);
%       ======================================
%       Purpose: Compute Bernoulli number Bn
%       Input :  n --- Serial number
%       Output:  BN(n)--- Bn
%       ======================================
tpi=6.283185307179586d0;
bn(0+1)=1.0d0;
bn(1+1)=-0.5d0;
bn(2+1)=1.0d0./6.0d0;
r1=(2.0d0./tpi).^2;
for  m=4:2:n;
r1=-r1.*(m-1).*m./(tpi.*tpi);
r2=1.0d0;
for  k=2:10000;
s=(1.0d0./k).^m;
r2=r2+s;
if(s < 1.0d-15)break; end;
end;
bn(m+1)=r1.*r2;
end;
return;
end

