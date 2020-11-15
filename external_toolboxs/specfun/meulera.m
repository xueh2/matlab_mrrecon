function meulera
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
%                subroutine EULERA
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
[n,e]=eulera(n,e);
fprintf(1,'%s \n','   n            en');
fprintf(1,'%s \n',' --------------------------');
for  k=0:2:n;
fprintf(1,[repmat(' ',1,2),'%3g','%22.12g' ' \n'],k,e(k+1));
end;  k=n+1;
%format(2x,i3,d22.12);
end
function [n,en]=eulera(n,en,varargin);
%       ======================================
%       Purpose: Compute Euler number En
%       Input :  n --- Serial number
%       Output:  EN(n)--- En
%       ======================================
en(0+1)=1.0d0;
for  m=1:n./2;
s=1.0d0;
for  k=1:m-1;
r=1.0d0;
for  j=1:2.*k;
r=r.*(2.0d0.*m-2.0d0.*k+j)./j;
end;  j=2.*k+1;
s=s+r.*en(2.*k+1);
end;  k=m-1+1;
en(2.*m+1)=-s;
end;  m=fix(n)./2+1;
return;
end

