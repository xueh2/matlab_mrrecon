function mlpmn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes the associated Legendre
%                functions Pmn(x)and their derivatives Pmn'(x)
%                using subroutine LPMN
%       Input :  x --- Argument of Pmn(x)
%                m --- Order of Pmn(x),  m = 0,1,2,...,n
%                n --- Degree of Pmn(x), n = 0,1,2,...,N
%       Output:  PM(m,n)--- Pmn(x)
%                PD(m,n)--- Pmn'(x)
%       Example: x = 0.50
%          Pmn(x):
%          m\n        1            2            3            4
%         --------------------------------------------------------
%           0      .500000     -.125000     -.437500     -.289063
%           1     -.866025    -1.299038     -.324760     1.353165
%           2      .000000     2.250000     5.625000     4.218750
%           3      .000000      .000000    -9.742786   -34.099750
%           4      .000000      .000000      .000000    59.062500
%          Pmn'(x):
%          m\n        1            2            3            4
%         --------------------------------------------------------
%           0     1.000000     1.500000      .375000    -1.562500
%           1      .577350    -1.732051    -6.278684    -5.773503
%           2      .000000    -3.000000     3.750000    33.750000
%           3      .000000      .000000    19.485572      .000000
%           4      .000000      .000000      .000000  -157.500000
%       ==========================================================
m=[];n=[];x=[];pm=[];pd=[];
 pm=zeros(100+1,100+1);
pd=zeros(100+1,100+1);
fprintf(1,'%s \n','  please enter m, n and x');
%        READ(*,*)M,N,X
m=2;
n=4;
x=.5;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  m     n      x          pmn(x)pmn''(x)');
fprintf(1,'%s \n',' ---------------------------------------------------');
[dumvar1,m,n,x,pm,pd]=lpmn(100,m,n,x,pm,pd);
for  j=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,3),'%3g',repmat(' ',1,3),'%5.1g',repmat('%17.8g',1,2) ' \n'],m,j,x,pm(m+1,j+1),pd(m+1,j+1));
end;  j=n+1;
%format(1x,i3,3x,i3,3x,f5.1,2e17.8);
end
function [mm,m,n,x,pm,pd]=lpmn(mm,m,n,x,pm,pd,varargin);
%       =====================================================
%       Purpose: Compute the associated Legendre functions
%                Pmn(x)and their derivatives Pmn'(x)
%       Input :  x  --- Argument of Pmn(x)
%                m  --- Order of Pmn(x),  m = 0,1,2,...,n
%                n  --- Degree of Pmn(x), n = 0,1,2,...,N
%                mm --- Physical dimension of PM and PD
%       Output:  PM(m,n)--- Pmn(x)
%                PD(m,n)--- Pmn'(x)
%       =====================================================
for  i=0:n;
for  j=0:m;
pm(j+1,i+1)=0.0d0;
pd(j+1,i+1)=0.0d0;
end;  j=fix(m)+1;
end;  i=fix(n)+1;
pm(0+1,0+1)=1.0d0;
if(abs(x)== 1.0d0);
for  i=1:n;
pm(0+1,i+1)=x.^i;
pd(0+1,i+1)=0.5d0.*i.*(i+1.0d0).*x.^(i+1);
end;  i=fix(n)+1;
for  j=1:n;
for  i=1:m;
if(i == 1);
pd(i+1,j+1)=1.0d+300;
elseif(i == 2);
pd(i+1,j+1)=-0.25d0.*(j+2).*(j+1).*j.*(j-1).*x.^(j+1);
end;
end;  i=fix(m)+1;
end;  j=fix(n)+1;
return;
end;
ls=1;
if(abs(x)> 1.0d0)ls=-1; end;
xq=sqrt(ls.*(1.0d0-x.*x));
xs=ls.*(1.0d0-x.*x);
for  i=1:m;
pm(i+1,i+1)=-ls.*(2.0d0.*i-1.0d0).*xq.*pm(i-1+1,i-1+1);
end;  i=fix(m)+1;
for  i=0:m;
pm(i+1,i+1+1)=(2.0d0.*i+1.0d0).*x.*pm(i+1,i+1);
end;  i=fix(m)+1;
for  i=0:m;
for  j=i+2:n;
pm(i+1,j+1)=((2.0d0.*j-1.0d0).*x.*pm(i+1,j-1+1)-(i+j-1.0d0).*pm(i+1,j-2+1))./(j-i);
end;  j=fix(n)+1;
end;  i=fix(m)+1;
pd(0+1,0+1)=0.0d0;
for  j=1:n;
pd(0+1,j+1)=ls.*j.*(pm(0+1,j-1+1)-x.*pm(0+1,j+1))./xs;
end;  j=fix(n)+1;
for  i=1:m;
for  j=i:n;
pd(i+1,j+1)=ls.*i.*x.*pm(i+1,j+1)./xs+(j+i).*(j-i+1.0d0)./xq.*pm(i-1+1,j+1);
end;  j=fix(n)+1;
end;  i=fix(m)+1;
return;
end

