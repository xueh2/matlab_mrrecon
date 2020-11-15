function mlpmns
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes the associated Legendre
%                functions Pmn(x)and their derivatives Pmn'(x)
%                for a given order using subroutine LPMNS
%       Input :  x --- Argument of Pmn(x)
%                m --- Order of Pmn(x),  m = 0,1,2,...,n
%                n --- Degree of Pmn(x), n = 0,1,2,...,N
%       Output:  PM(n)--- Pmn(x)
%                PD(n)--- Pmn'(x)
%       Examples:
%                m = 1,  N = 5,  x = .5
%                n        Pmn(x)Pmn'(x)
%               -------------------------------------
%                0    .00000000D+00    .00000000D+00
%                1    .86602540D+00   -.57735027D+00
%                2    .12990381D+01    .17320508D+01
%                3    .32475953D+00    .62786842D+01
%                4   -.13531647D+01    .57735027D+01
%                5   -.19282597D+01   -.43977853D+01
%                m = 2,  N = 6,  x = 2.5
%                n        Pmn(x)Pmn'(x)
%               -------------------------------------
%                0    .00000000D+00    .00000000D+00
%                1    .00000000D+00    .00000000D+00
%                2    .15750000D+02    .15000000D+02
%                3    .19687500D+03    .26625000D+03
%                4    .16832813D+04    .29812500D+04
%                5    .12230859D+05    .26876719D+05
%                6    .81141416D+05    .21319512D+06
%       =======================================================
m=[];n=[];x=[];pm=[];pd=[];
 pm=zeros(1,200+1);
pd=zeros(1,200+1);
fprintf(1,'%s \n','please enter m, n, and x ');
%        READ(*,*)M,N,X
m=1;
n=5;
x=.5;
fprintf(1,[repmat(' ',1,1),'m =','%2g',',  ','n =','%2g',',  ','x =','%5.1g' ' \n'],m,n,x);
[m,n,x,pm,pd]=lpmns(m,n,x,pm,pd);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n        pmn(x)pmn''(x)');
fprintf(1,'%s \n',' -------------------------------------');
for  j=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,2) ' \n'],j,pm(j+1),pd(j+1));
end;  j=n+1;
%format(1x,i3,2d17.8);
%format(1x,',i2,',  ',',i2,',  ',',f5.1);
end
function [m,n,x,pm,pd]=lpmns(m,n,x,pm,pd,varargin);
%       ========================================================
%       Purpose: Compute associated Legendre functions Pmn(x)
%                and Pmn'(x)for a given order
%       Input :  x --- Argument of Pmn(x)
%                m --- Order of Pmn(x),  m = 0,1,2,...,n
%                n --- Degree of Pmn(x), n = 0,1,2,...,N
%       Output:  PM(n)--- Pmn(x)
%                PD(n)--- Pmn'(x)
%       ========================================================
for  k=0:n;
pm(k+1)=0.0d0;
pd(k+1)=0.0d0;
end;  k=fix(n)+1;
if(abs(x)== 1.0d0);
for  k=0:n;
if(m == 0);
pm(k+1)=1.0d0;
pd(k+1)=0.5d0.*k.*(k+1.0);
if(x < 0.0);
pm(k+1)=(-1).^k.*pm(k+1);
pd(k+1)=(-1).^(k+1).*pd(k+1);
end;
elseif(m == 1);
pd(k+1)=1.0d+300;
elseif(m == 2);
pd(k+1)=-0.25d0.*(k+2.0).*(k+1.0).*k.*(k-1.0);
if(x < 0.0)pd(k+1)=(-1).^(k+1).*pd(k+1); end;
end;
end;  k=fix(n)+1;
return;
end;
x0=abs(1.0d0-x.*x);
pm0=1.0d0;
pmk=pm0;
for  k=1:m;
pmk=(2.0d0.*k-1.0d0).*sqrt(x0).*pm0;
pm0=pmk;
end;  k=fix(m)+1;
pm1=(2.0d0.*fix(m)+1.0d0).*x.*pm0;
pm(m+1)=pmk;
pm(m+1+1)=pm1;
for  k=m+2:n;
pm2=((2.0d0.*k-1.0d0).*x.*pm1-(k+fix(m)-1.0d0).*pmk)./(k-fix(m));
pm(k+1)=pm2;
pmk=pm1;
pm1=pm2;
end;  k=fix(n)+1;
pd(0+1)=((1.0d0-fix(m)).*pm(1+1)-x.*pm(0+1))./(x.*x-1.0);
for  k=1:n;
pd(k+1)=(k.*x.*pm(k+1)-(k+fix(m)).*pm(k-1+1))./(x.*x-1.0d0);
end;  k=fix(n)+1;
return;
end

