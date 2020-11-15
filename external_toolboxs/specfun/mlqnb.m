function mlqnb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===============================================================
%       Purpose: This program computes the Legendre functions Qn(x)
%                and Qn'(x)using subroutine LQNB
%       Input :  x  --- Argument of Qn(x)
%                n  --- Degree of Qn(x,n = 0,1,תתת)
%       Output:  QN(n)--- Qn(x)
%                QD(n)--- Qn'(x)
%       Examples:     x1 = 0.50,    x2 = 2.50
%       n      Qn(x1)Qn'(x1)Qn(x2)Qn'(x2)
%     ----------------------------------------------------------------
%       0     .54930614    1.33333333   .42364893D+00  -.19047619D+00
%       1    -.72534693    1.21597281   .59122325D-01  -.52541546D-01
%       2    -.81866327    -.84270745   .98842555D-02  -.13109214D-01
%       3    -.19865477   -2.87734353   .17695141D-02  -.31202687D-02
%       4     .44017453   -2.23329085   .32843271D-03  -.72261513D-03
%       5     .55508089    1.08422720   .62335892D-04  -.16437427D-03
%       ===============================================================
n=[];x=[];qn=[];qd=[];
 qn=zeros(1,100+1);
qd=zeros(1,100+1);
fprintf(1,'%s \n','please enter nmax and x ');
%        READ(*,*)N,X
n=5;
x=.5;
fprintf(1,[repmat(' ',1,3),'x =','%5.2g' ' \n'],x);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n          qn(x)qn''(x)');
fprintf(1,'%s \n','--------------------------------------');
[n,x,qn,qd]=lqnb(n,x,qn,qd);
for  k=0:n;
if(abs(x)< 1.0);
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,2) ' \n'],k,qn(k+1),qd(k+1));
else;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,2) ' \n'],k,qn(k+1),qd(k+1));
end;
end;  k=n+1;
%format(1x,i3,2f17.8);
%format(1x,i3,2d17.8);
%format(3x,',f5.2);
end
function [n,x,qn,qd]=lqnb(n,x,qn,qd,varargin);
%       ====================================================
%       Purpose: Compute Legendre functions Qn(x)& Qn'(x)
%       Input :  x  --- Argument of Qn(x)
%                n  --- Degree of Qn(x,n = 0,1,2,תתת)
%       Output:  QN(n)--- Qn(x)
%                QD(n)--- Qn'(x)
%       ====================================================
eps=1.0d-14;
if(abs(x)== 1.0d0);
for  k=0:n;
qn(k+1)=1.0d+300;
qd(k+1)=1.0d+300;
end;  k=fix(n)+1;
return;
end;
if(x <= 1.021d0);
x2=abs((1.0d0+x)./(1.0d0-x));
q0=0.5d0.*log(x2);
q1=x.*q0-1.0d0;
qn(0+1)=q0;
qn(1+1)=q1;
qd(0+1)=1.0d0./(1.0d0-x.*x);
qd(1+1)=qn(0+1)+x.*qd(0+1);
for  k=2:n;
qf=((2.0d0.*k-1.0d0).*x.*q1-(k-1.0d0).*q0)./k;
qn(k+1)=qf;
qd(k+1)=(qn(k-1+1)-x.*qf).*k./(1.0d0-x.*x);
q0=q1;
q1=qf;
end;  k=fix(n)+1;
else;
qc2=1.0d0./x;
for  j=1:n;
qc2=qc2.*j./((2.0.*j+1.0d0).*x);
if(j == n-1)qc1=qc2; end;
end;  j=fix(n)+1;
for  l=0:1;
nl=fix(n)+l;
qf=1.0d0;
qr=1.0d0;
for  k=1:500;
qr=qr.*(0.5d0.*nl+k-1.0d0).*(0.5d0.*(nl-1)+k)./((nl+k-0.5d0).*k.*x.*x);
qf=qf+qr;
if(abs(qr./qf)< eps)break; end;
end;
if(l == 0);
qn(n-1+1)=qf.*qc1;
else;
qn(n+1)=qf.*qc2;
end;
end;
qf2=qn(fix(n)+1);
qf1=qn(fix(n)-1+1);
for  k=n:-1:2;
qf0=((2.*k-1.0d0).*x.*qf1-k.*qf2)./(k-1.0d0);
qn(k-2+1)=qf0;
qf2=qf1;
qf1=qf0;
end;  k=2-1;
qd(0+1)=1.0d0./(1.0d0-x.*x);
for  k=1:n;
qd(k+1)=k.*(qn(k-1+1)-x.*qn(k+1))./(1.0d0-x.*x);
end;  k=fix(n)+1;
end;
return;
end

