function mlqna
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ======================================================
%       Purpose: This program computes the Legendre functions
%                Qn(x)and Qn'(x)using subroutine LQNA
%       Input :  x  --- Argument of Qn(x,-1 ף x ף 1)
%                n  --- Degree of Qn(x,n = 0,1,תתת)
%       Output:  QN(n)--- Qn(x)
%                QD(n)--- Qn'(x)
%       Example:  x = 0.50
%                 n        Qn(x)Qn'(x)
%                ---------------------------------
%                 0      .54930614     1.33333333
%                 1     -.72534693     1.21597281
%                 2     -.81866327     -.84270745
%                 3     -.19865477    -2.87734353
%                 4      .44017453    -2.23329085
%                 5      .55508089     1.08422720
%       ======================================================
n=[];x=[];qn=[];qd=[];
  qn=0;
 qd=0;
 x=0;
 qn=zeros(1,100+1);
qd=zeros(1,100+1);
fprintf(1,'%s \n','  please enter nmax and x');
%        READ(*,*)N,X
n=5;
x=.5;
fprintf(1,[repmat(' ',1,3),'x =','%5.2g' ' \n'],x);
fprintf(1,'%0.15g \n');
[n,x,qn,qd]=lqna(n,x,qn,qd);
fprintf(1,'%s \n','  n        qn(x)qn''(x)');
fprintf(1,'%s \n',' ---------------------------------');
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%15.8g',1,2) ' \n'],k,qn(k+1),qd(k+1));
end;  k=n+1;
%format(1x,i3,2f15.8);
%format(3x,',f5.2);
end
function [n,x,qn,qd]=lqna(n,x,qn,qd,varargin);
%       =====================================================
%       Purpose: Compute Legendre functions Qn(x)and Qn'(x)
%       Input :  x  --- Argument of Qn(x,-1 ף x ף 1)
%                n  --- Degree of Qn(x,n = 0,1,2,תתת)
%       Output:  QN(n)--- Qn(x)
%                QD(n)--- Qn'(x)
%(1.0D+300 stands for infinity)
%       =====================================================
if(abs(x)== 1.0d0);
for  k=0:n;
qn(k+1)=1.0d+300;
qd(k+1)=-1.0d+300;
end;  k=fix(n)+1;
elseif(abs(x)< 1.0d0);
q0=0.5d0.*log((1.0d0+x)./(1.0d0-x));
q1=x.*q0-1.0d0;
qn(0+1)=q0;
qn(1+1)=q1;
qd(0+1)=1.0d0./(1.0d0-x.*x);
qd(1+1)=qn(0+1)+x.*qd(0+1);
for  k=2:n;
qf=((2.*k-1).*x.*q1-(k-1).*q0)./k;
qn(k+1)=qf;
qd(k+1)=(qn(k-1+1)-x.*qf).*k./(1.0d0-x.*x);
q0=q1;
q1=qf;
end;  k=fix(n)+1;
end;
return;
end

