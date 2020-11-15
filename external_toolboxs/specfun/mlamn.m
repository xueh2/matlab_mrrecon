function mlamn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ====================================================
%       Purpose: This program computes the lambda functions
%                and their derivatives using subroutine
%                LAMN
%       Input:   x --- Argument of lambda function
%                n --- Order of lambda function
%(n = 0,1,..., n ó 250)
%       Output:  BL(n)--- Lambda function of order n
%                DL(n)--- Derivative of lambda function
%       Example: Nmax = 5,  x = 10.00
%                 n       lambda(x)lambda'(x)
%                ---------------------------------------
%                 0    -.24593576D+00    -.43472746D-01
%                 1     .86945492D-02    -.50926063D-01
%                 2     .20370425D-01    -.46703503D-02
%                 3     .28022102D-02     .10540929D-01
%                 4    -.84327431D-02     .89879627D-02
%                 5    -.89879627D-02     .55521954D-03
%       ====================================================
n=[];x=[];nm=[];bl=[];dl=[];
 bl=zeros(1,250+1);
dl=zeros(1,250+1);
fprintf(1,'%s \n','  please enter n,x = ?');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,1),'n =','%4g',repmat(' ',1,6),'x =','%8.2g' ' \n'],n,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,x,nm,bl,dl]=lamn(n,x,nm,bl,dl);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n       lambda(x)lambda''(x)');
fprintf(1,'%s \n',' ---------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%18.8g',1,2) ' \n'],k,bl(k+1),dl(k+1));
end;  k=nm+1;
%format(1x,,i4,6x,,f8.2);
%format(1x,i3,2d18.8);
end
function [n,x,nm,bl,dl]=lamn(n,x,nm,bl,dl,varargin);
%       =========================================================
%       Purpose: Compute lambda functions and their derivatives
%       Input:   x --- Argument of lambda function
%                n --- Order of lambda function
%       Output:  BL(n)--- Lambda function of order n
%                DL(n)--- Derivative of lambda function
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 for computing the start
%                point for backward recurrence
%       =========================================================
nm=fix(fix(n));
if(abs(x)< 1.0d-100);
for  k=0:n;
bl(k+1)=0.0d0;
dl(k+1)=0.0d0;
end;  k=fix(n)+1;
bl(0+1)=1.0d0;
dl(1+1)=0.5d0;
return;
end;
if(x <= 12.0d0);
x2=x.*x;
for  k=0:n;
bk=1.0d0;
r=1.0d0;
for  i=1:50;
r=-0.25d0.*r.*x2./(i.*(i+k));
bk=bk+r;
if(abs(r)< abs(bk).*1.0d-15)break; end;
end;
bl(k+1)=bk;
if(k >= 1)dl(k-1+1)=-0.5d0.*x./k.*bk; end;
end;
uk=1.0d0;
r=1.0d0;
for  i=1:50;
r=-0.25d0.*r.*x2./(i.*(i+fix(n)+1.0d0));
uk=uk+r;
if(abs(r)< abs(uk).*1.0d-15)break; end;
end;
dl(n+1)=-0.5d0.*x./(fix(n)+1.0d0).*uk;
return;
end;
if(n == 0)nm=1; end;
m=msta1(x,200);
if(m < nm);
nm=fix(m);
else;
m=msta2(x,fix(nm),15);
end;
bs=0.0d0;
f0=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0).*f1./x-f0;
if(k <= nm)bl(k+1)=f; end;
if(k == 2.*fix(k./2))bs=bs+2.0d0.*f; end;
f0=f1;
f1=f;
end;  k=0-1;
bg=bs-f;
for  k=0:nm;
bl(k+1)=bl(k+1)./bg;
end;  k=fix(nm)+1;
r0=1.0d0;
for  k=1:nm;
r0=2.0d0.*r0.*k./x;
bl(k+1)=r0.*bl(k+1);
end;  k=fix(nm)+1;
dl(0+1)=-0.5d0.*x.*bl(1+1);
for  k=1:nm;
dl(k+1)=2.0d0.*k./x.*(bl(k-1+1)-bl(k+1));
end;  k=fix(nm)+1;
return;
end
function [msta1Result]=msta1(x,mp,varargin);
%       ===================================================
%       Purpose: Determine the starting point for backward
%                recurrence such that the magnitude of
%                Jn(x)at that point is about 10^(-MP)
%       Input :  x     --- Argument of Jn(x)
%                MP    --- Value of magnitude
%       Output:  MSTA1 --- Starting point
%       ===================================================
a0=abs(x);
n0=fix(1.1.*a0)+1;
f0=envj(n0,a0)-fix(mp);
n1=n0+5;
f1=envj(n1,a0)-fix(mp);
for  it=1:20;
nn=n1-(n1-n0)./(1.0d0-f0./f1);
f=envj(nn,a0)-fix(mp);
if(abs(nn-n1)< 1)break; end;
n0=n1;
f0=f1;
n1=nn;
f1=f;
end;
msta1Result=fix(nn);
return;
end
function [msta2Result]=msta2(x,n,mp,varargin);
%       ===================================================
%       Purpose: Determine the starting point for backward
%                recurrence such that all Jn(x)has MP
%                significant digits
%       Input :  x  --- Argument of Jn(x)
%                n  --- Order of Jn(x)
%                MP --- Significant digit
%       Output:  MSTA2 --- Starting point
%       ===================================================
a0=abs(x);
hmp=0.5d0.*fix(mp);
ejn=envj(fix(n),a0);
if(ejn <= hmp);
obj=fix(mp);
n0=fix(1.1.*a0);
else;
obj=hmp+ejn;
n0=fix(n);
end;
f0=envj(n0,a0)-obj;
n1=n0+5;
f1=envj(n1,a0)-obj;
for  it=1:20;
nn=n1-(n1-n0)./(1.0d0-f0./f1);
f=envj(nn,a0)-obj;
if(abs(nn-n1)< 1)break; end;
n0=n1;
f0=f1;
n1=nn;
f1=f;
end;
msta2Result=fix(nn+10);
return;
end
function [envjResult]=envj(n,x,varargin);
envjResult=0.5d0.*log10(6.28d0.*fix(n))-fix(n).*log10(1.36d0.*x./fix(n));
return;
end

