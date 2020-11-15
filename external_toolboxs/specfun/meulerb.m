function meulerb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes Euler number En using
%                subroutine EULERB
%       Example: Compute Euler number En for n = 0,2,...,10
%                Computed results:
%                   n            En
%                 --------------------------
%                   0     .100000000000D+01
%                   2    -.100000000000D+01
%                   4     .500000000000D+01
%                   6    -.610000000000D+02
%                   8     .138500000000D+04
%                  10    -.505210000000D+05
%       ==========================================================
n=[];e=[];
  e=0;
 e=zeros(1,200+1);
fprintf(1,'%s \n','  please enter nmax ');
%        READ(*,*)N
n=10;
[n,e]=eulerb(n,e);
fprintf(1,'%s \n','   n            en');
fprintf(1,'%s \n',' --------------------------');
for  k=0:2:n;
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],k,e(k+1));
end;  k=n+1;
%format(2x,i3,d22.12);
end
function [n,en]=eulerb(n,en,varargin);
%       ======================================
%       Purpose: Compute Euler number En
%       Input :  n --- Serial number
%       Output:  EN(n)--- En
%       ======================================
hpi=2.0d0./3.141592653589793d0;
en(0+1)=1.0d0;
en(2+1)=-1.0d0;
r1=-4.0d0.*hpi.^3;
for  m=4:2:n;
r1=-r1.*(m-1).*m.*hpi.*hpi;
r2=1.0d0;
isgn=1.0d0;
for  k=3:2:1000;
isgn=-isgn;
s=(1.0d0./k).^(m+1);
r2=r2+isgn.*s;
if(s < 1.0d-15)break; end;
end;
en(m+1)=r1.*r2;
end;
return;
end

