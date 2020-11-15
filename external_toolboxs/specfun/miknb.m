function miknb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes modified Bessel functions
%                In(x)and Kn(x), and their derivatives using
%                subroutine IKNB
%       Input:   x --- Argument of In(x)and Kn(x,0 ó x ó 700)
%                n --- Order of In(x)and Kn(x)
%(n = 0,1,..., n ó 250)
%       Output:  BI(n)--- In(x)
%                DI(n)--- In'(x)
%                BK(n)--- Kn(x)
%                DK(n)--- Kn'(x)
%       Example: Nmax = 5,    x = 10.0
%     n      In(x)In'(x)Kn(x)Kn'(x)
%    ---------------------------------------------------------------
%     0   .2815717D+04   .2670988D+04   .1778006D-04  -.1864877D-04
%     1   .2670988D+04   .2548618D+04   .1864877D-04  -.1964494D-04
%     2   .2281519D+04   .2214685D+04   .2150982D-04  -.2295074D-04
%     3   .1758381D+04   .1754005D+04   .2725270D-04  -.2968563D-04
%     4   .1226491D+04   .1267785D+04   .3786144D-04  -.4239728D-04
%     5   .7771883D+03   .8378964D+03   .5754185D-04  -.6663236D-04
%       =============================================================
n=[];x=[];nm=[];bi=[];di=[];bk=[];dk=[];
 bi=zeros(1,250+1);
di=zeros(1,250+1);
bk=zeros(1,250+1);
dk=zeros(1,250+1);
fprintf(1,'%s \n','  please enter n, x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',    ','x =','%5.1g' ' \n'],n,x);
fprintf(1,'%0.15g \n');
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,x,nm,bi,di,bk,dk]=iknb(n,x,nm,bi,di,bk,dk);
fprintf(1,'%s ','  n      in(x)in''(x)');fprintf(1,'%s \n', '        kn(x)kn''(x)');
fprintf(1,'%s ',' -------------------------------');fprintf(1,'%s \n','--------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%15.7g',1,4) ' \n'],k,bi(k+1),di(k+1),bk(k+1),dk(k+1));
end;  k=nm+1;
%format(1x,i3,4d15.7);
%format(3x,,i3,',    ',,f5.1);
end
function [n,x,nm,bi,di,bk,dk]=iknb(n,x,nm,bi,di,bk,dk,varargin);
%       ============================================================
%       Purpose: Compute modified Bessel functions In(x)and Kn(x),
%                and their derivatives
%       Input:   x --- Argument of In(x)and Kn(x,0 ó x ó 700)
%                n --- Order of In(x)and Kn(x)
%       Output:  BI(n)--- In(x)
%                DI(n)--- In'(x)
%                BK(n)--- Kn(x)
%                DK(n)--- Kn'(x)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 for computing the starting point
%                for backward recurrence
%       ===========================================================
pi=3.141592653589793d0;
el=0.5772156649015329d0;
nm=fix(fix(n));
if(x <= 1.0d-100);
for  k=0:n;
bi(k+1)=0.0d0;
di(k+1)=0.0d0;
bk(k+1)=1.0d+300;
dk(k+1)=-1.0d+300;
end;  k=fix(n)+1;
bi(0+1)=1.0d0;
di(1+1)=0.5d0;
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
sk0=0.0d0;
f0=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0)./x.*f1+f0;
if(k <= nm)bi(k+1)=f; end;
if(k ~= 0&k == 2.*fix(k./2))sk0=sk0+4.0d0.*f./k; end;
bs=bs+2.0d0.*f;
f0=f1;
f1=f;
end;  k=0-1;
s0=exp(x)./(bs-f);
for  k=0:nm;
bi(k+1)=s0.*bi(k+1);
end;  k=fix(nm)+1;
if(x <= 8.0d0);
bk(0+1)=-(log(0.5d0.*x)+el).*bi(0+1)+s0.*sk0;
bk(1+1)=(1.0d0./x-bi(1+1).*bk(0+1))./bi(0+1);
else;
a0=sqrt(pi./(2.0d0.*x)).*exp(-x);
k0=16;
if(x >= 25.0)k0=10; end;
if(x >= 80.0)k0=8; end;
if(x >= 200.0)k0=6; end;
for  l=0:1;
bkl=1.0d0;
vt=4.0d0.*l;
r=1.0d0;
for  k=1:k0;
r=0.125d0.*r.*(vt-(2.0.*k-1.0).^2)./(k.*x);
bkl=bkl+r;
end;  k=k0+1;
bk(l+1)=a0.*bkl;
end;  l=1+1;
end;
g0=bk(0+1);
g1=bk(1+1);
for  k=2:nm;
g=2.0d0.*(k-1.0d0)./x.*g1+g0;
bk(k+1)=g;
g0=g1;
g1=g;
end;  k=fix(nm)+1;
di(0+1)=bi(1+1);
dk(0+1)=-bk(1+1);
for  k=1:nm;
di(k+1)=bi(k-1+1)-k./x.*bi(k+1);
dk(k+1)=-bk(k-1+1)-k./x.*bk(k+1);
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
%       Input :  x     --- Argument of Jn(x)
%                n     --- Order of Jn(x)
%                MP    --- Significant digit
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

