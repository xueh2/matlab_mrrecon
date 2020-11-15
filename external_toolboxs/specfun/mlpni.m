function mlpni
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ========================================================
%     Purpose: This program computes the Legendre polynomials
%     Pn(x), Pn'(x)and the integral of Pn(t)from 0
%     to x using subroutine LPNI
%     Input :  x --- Argument of Pn(x)
%     n --- Degree of Pn(x,n = 0,1,...)
%     Output:  PN(n)--- Pn(x)
%     PD(n)--- Pn'(x)
%     PL(n)--- Integral of Pn(t)from 0 to x
%     Example: x = 0.50
%     n       Pn(x)Pn'(x)Pn(t)dt
%     ---------------------------------------------
%     0    1.00000000     .00000000     .50000000
%     1     .50000000    1.00000000     .12500000
%     2    -.12500000    1.50000000    -.18750000
%     3    -.43750000     .37500000    -.14843750
%     4    -.28906250   -1.56250000     .05859375
%     5     .08984375   -2.22656250     .11816406
%     ========================================================
n=[];x=[];pn=[];pd=[];pl=[];
  pn=0;
 pd=0;
 pl=0;
 x=0;
 pn=zeros(1,100+1);
pd=zeros(1,100+1);
pl=zeros(1,100+1);
fprintf(1,'%s \n','  please enter nmax and x');
%     READ(*,*)N,X
n=5;
x=.5;
fprintf(1,[repmat(' ',1,3),'x =','%5.2g' ' \n'],x);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n        pn(x)pn''(x)pn(t)dt');
fprintf(1,'%s \n',' ---------------------------------------------------');
[n,x,pn,pd,pl]=lpni(n,x,pn,pd,pl);
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%16.8g',1,3) ' \n'],k,pn(k+1),pd(k+1),pl(k+1));
end;  k=n+1;
%format(1x,i3,3e16.8);
%format(3x,',f5.2);
end
function [n,x,pn,pd,pl]=lpni(n,x,pn,pd,pl,varargin);
%     =====================================================
%     Purpose: Compute Legendre polynomials Pn(x), Pn'(x)
%     and the integral of Pn(t)from 0 to x
%     Input :  x --- Argument of Pn(x)
%     n --- Degree of Pn(x,n = 0,1,...)
%     Output:  PN(n)--- Pn(x)
%     PD(n)--- Pn'(x)
%     PL(n)--- Integral of Pn(t)from 0 to x
%     =====================================================
pn(0+1)=1.0d0;
pn(1+1)=x;
pd(0+1)=0.0d0;
pd(1+1)=1.0d0;
pl(0+1)=x;
pl(1+1)=0.5d0.*x.*x;
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
pl(k+1)=(x.*pn(k+1)-pn(k-1+1))./(k+1.0d0);
p0=p1;
p1=pf;
if(~(k == 2.*fix(k./2)));
r=1.0d0./(k+1.0d0);
n1=(k-1)./2;
for  j=1:n1;
r=(0.5d0./j-1.0d0).*r;
end;  j=n1+1;
pl(k+1)=pl(k+1)+r;
end;
end;  k=fix(n)+1;
return;
end

