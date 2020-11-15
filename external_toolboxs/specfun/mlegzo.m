function mlegzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose : This program computes the zeros of Legendre
%                 polynomial Pn(x)in the interval[-1,1]and the
%                 corresponding weighting coefficients for Gauss-
%                 Legendre integration using subroutine LEGZO
%       Input :   n    --- Order of the Legendre polynomial
%       Output:   X(n)--- Zeros of the Legendre polynomial
%                 W(n)--- Corresponding weighting coefficients
%       ============================================================
n=[];x=[];w=[];
 x=zeros(1,120);
w=zeros(1,120);
fprintf(1,'%s \n','please enter the order of pn(x), n ');
%        READ(*,*)N
n=5;
fprintf(1,[repmat(' ',1,1),'n =','%3g' ' \n'],n);
[n,x,w]=legzo(n,x,w);
fprintf(1,'%s \n','  nodes and weights for gauss-legendre integration');
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  i              xi                   wi');
fprintf(1,'%s \n',' ------------------------------------------------');
for  i=1:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,1),'%22.13g','%22.13g' ' \n'],i,x(i),w(i));
end;  i=n+1;
%format(1x,',i3);
%format(1x,i3,1x,f22.13,d22.13);
end
function [n,x,w]=legzo(n,x,w,varargin);
%       =========================================================
%       Purpose : Compute the zeros of Legendre polynomial Pn(x)
%                 in the interval[-1,1], and the corresponding
%                 weighting coefficients for Gauss-Legendre
%                 integration
%       Input :   n    --- Order of the Legendre polynomial
%       Output:   X(n)--- Zeros of the Legendre polynomial
%                 W(n)--- Corresponding weighting coefficients
%       =========================================================
n0=(fix(n)+1)./2;
for  nr=1:n0;
z=cos(3.1415926d0.*(nr-0.25d0)./fix(n));
while (1);
z0=z;
p=1.0d0;
for  i=1:nr-1;
p=p.*(z-x(i));
end;  i=nr-1+1;
f0=1.0d0;
if(nr == n0&n ~= 2.*fix(n./2))z=0.0d0; end;
f1=z;
for  k=2:n;
pf=(2.0d0-1.0d0./k).*z.*f1-(1.0d0-1.0d0./k).*f0;
pd=k.*(f1-z.*pf)./(1.0d0-z.*z);
f0=f1;
f1=pf;
end;  k=fix(n)+1;
if(z == 0.0)break; end;
fd=pf./p;
q=0.0d0;
for  i=1:nr-1;
wp=1.0d0;
for  j=1:nr-1;
if(j ~= i)wp=wp.*(z-x(j)); end;
end;  j=nr-1+1;
q=q+wp;
end;  i=nr-1+1;
gd=(pd-q.*fd)./p;
z=z-fd./gd;
if(~(abs(z-z0)> abs(z).*1.0d-15))break; end;
end;
x(nr)=z;
x(n+1-nr)=-z;
w(nr)=2.0d0./((1.0d0-z.*z).*pd.*pd);
w(n+1-nr)=w(nr);
end;
return;
end

