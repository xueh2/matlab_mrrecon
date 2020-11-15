function msphi
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ======================================================
%       Purpose: This program computes the modified spherical
%                Bessel functions of the first kind in(x)and
%                in'(x)using subroutine SPHI
%       Input :  x --- Argument of in(x)
%                n --- Order of in(x,0 ó n ó 250)
%       Output:  SI(n)--- in(x)
%                DI(n)--- in'(x)
%       Example: x = 10.0
%                  n          in(x)in'(x)
%                --------------------------------------------
%                  0     .1101323287D+04     .9911909633D+03
%                  1     .9911909633D+03     .9030850948D+03
%                  2     .8039659985D+03     .7500011637D+03
%                  3     .5892079640D+03     .5682828129D+03
%                  4     .3915204237D+03     .3934477522D+03
%                  5     .2368395827D+03     .2494166741D+03
%       ======================================================
n=[];x=[];nm=[];si=[];di=[];
 si=zeros(1,250+1);
di=zeros(1,250+1);
fprintf(1,'%s \n','please enter n and x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',     ','x =','%6.1g' ' \n'],n,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,x,nm,si,di]=sphi(n,x,nm,si,di);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n          in(x)in''(x)');
fprintf(1,'%s \n','--------------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%20.10g',1,2) ' \n'],k,si(k+1),di(k+1));
end;  k=nm+1;
%format(1x,i3,2d20.10);
%format(3x,',i3,',     ',',f6.1);
end
function [n,x,nm,si,di]=sphi(n,x,nm,si,di,varargin);
%       ========================================================
%       Purpose: Compute modified spherical Bessel functions
%                of the first kind, in(x)and in'(x)
%       Input :  x --- Argument of in(x)
%                n --- Order of in(x,n = 0,1,2,...)
%       Output:  SI(n)--- in(x)
%                DI(n)--- in'(x)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       ========================================================
nm=fix(fix(n));
if(abs(x)< 1.0d-100);
for  k=0:n;
si(k+1)=0.0d0;
di(k+1)=0.0d0;
end;  k=fix(n)+1;
si(0+1)=1.0d0;
di(1+1)=0.333333333333333d0;
return;
end;
si(0+1)=sinh(x)./x;
si(1+1)=-(sinh(x)./x-cosh(x))./x;
si0=si(0+1);
if(n >= 2);
m=msta1(x,200);
if(m < n);
nm=fix(m);
else;
m=msta2(x,fix(n),15);
end;
f0=0.0d0;
f1=1.0d0-100;
for  k=m:-1:0;
f=(2.0d0.*k+3.0d0).*f1./x+f0;
if(k <= nm)si(k+1)=f; end;
f0=f1;
f1=f;
end;  k=0-1;
cs=si0./f;
for  k=0:nm;
si(k+1)=cs.*si(k+1);
end;  k=fix(nm)+1;
end;
di(0+1)=si(1+1);
for  k=1:nm;
di(k+1)=si(k-1+1)-(k+1.0d0)./x.*si(k+1);
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

