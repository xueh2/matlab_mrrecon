function mlpn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes the Legendre polynomials
%                Pn(x)and their derivatives Pn'(x)using
%                subroutine LPN
%       Input :  x --- Argument of Pn(x)
%                n --- Degree of Pn(x,n = 0,1,...)
%       Output:  PN(n)--- Pn(x)
%                PD(n)--- Pn'(x)
%       Example:    x = 0.5
%                  n          Pn(x)Pn'(x)
%                ---------------------------------------
%                  0       1.00000000        .00000000
%                  1        .50000000       1.00000000
%                  2       -.12500000       1.50000000
%                  3       -.43750000        .37500000
%                  4       -.28906250      -1.56250000
%                  5        .08984375      -2.22656250
%       ========================================================
n=[];x=[];pn=[];pd=[];
  pn=0;
 pd=0;
 x=0;
 pn=zeros(1,100+1);
pd=zeros(1,100+1);
fprintf(1,'%s \n','  please enter nmax and x ');
%        READ(*,*)N,X
n=5;
x=.5;
fprintf(1,[repmat(' ',1,3),'x =','%5.1g' ' \n'],x);
fprintf(1,'%0.15g \n');
[n,x,pn,pd]=lpn(n,x,pn,pd);
fprintf(1,'%s \n','  n         pn(x)pn''(x)');
fprintf(1,'%s \n','---------------------------------------');
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,2) ' \n'],k,pn(k+1),pd(k+1));
end;  k=n+1;
%format(1x,i3,2e17.8);
%format(3x,',f5.1);
end
function [n,x,pn,pd]=lpn(n,x,pn,pd,varargin);
%       ===============================================
%       Purpose: Compute Legendre polynomials Pn(x)
%                and their derivatives Pn'(x)
%       Input :  x --- Argument of Pn(x)
%                n --- Degree of Pn(x,n = 0,1,...)
%       Output:  PN(n)--- Pn(x)
%                PD(n)--- Pn'(x)
%       ===============================================
pn(0+1)=1.0d0;
pn(1+1)=x;
pd(0+1)=0.0d0;
pd(1+1)=1.0d0;
p0=1.0d0;
p1=x;
for  k=2:n;
pf=(2.0d0.*k-1.0d0)./k.*x.*p1-(k-1.0d0)./k.*p0;
pn(k+1)=pf;
if(abs(x)== 1.0d0);
pd(k+1)=0.5d0.*x.^(k+1).*k.*(k+1.0d0);
else;
pd(k+1)=k.*(p1-x.*pf)./(1.0d0-x.*x);
end;
p0=p1;
p1=pf;
end;  k=fix(n)+1;
return;
end

