function mikna
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
%                subroutine IKNA
%       Input:   x --- Argument of In(x)and Kn(x,x ע 0)
%                n --- Order of In(x)and Kn(x)
%(n = 0,1,תתת, n ף 250)
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
[n,x,nm,bi,di,bk,dk]=ikna(n,x,nm,bi,di,bk,dk);
fprintf(1,'%s ','  n      in(x)in''(x)');fprintf(1,'%s \n', '        kn(x)kn''(x)');
fprintf(1,'%s ',' -------------------------------');fprintf(1,'%s \n','--------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%15.7g',1,4) ' \n'],k,bi(k+1),di(k+1),bk(k+1),dk(k+1));
end;  k=nm+1;
%format(1x,i3,4d15.7);
%format(3x,,i3,',    ',,f5.1);
end
function [n,x,nm,bi,di,bk,dk]=ikna(n,x,nm,bi,di,bk,dk,varargin);
%       ========================================================
%       Purpose: Compute modified Bessel functions In(x)and
%                Kn(x), and their derivatives
%       Input:   x --- Argument of In(x)and Kn(x,x ע 0)
%                n --- Order of In(x)and Kn(x)
%       Output:  BI(n)--- In(x)
%                DI(n)--- In'(x)
%                BK(n)--- Kn(x)
%                DK(n)--- Kn'(x)
%                NM --- Highest order computed
%       Routines called:
%(1)IK01A for computing I0(x),I1(x),K0(x)& K1(x)
%(2)MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       ========================================================
bi0=[];di0=[];bi1=[];di1=[];bk0=[];dk0=[];bk1=[];dk1=[];
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
[x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1]=ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1);
bi(0+1)=bi0;
bi(1+1)=bi1;
bk(0+1)=bk0;
bk(1+1)=bk1;
di(0+1)=di0;
di(1+1)=di1;
dk(0+1)=dk0;
dk(1+1)=dk1;
if(n <= 1)return; end;
if(x > 40.0&n < fix(0.25.*x));
h0=bi0;
h1=bi1;
for  k=2:n;
h=-2.0d0.*(k-1.0d0)./x.*h1+h0;
bi(k+1)=h;
h0=h1;
h1=h;
end;  k=fix(n)+1;
else;
m=msta1(x,200);
if(m < n);
nm=fix(m);
else;
m=msta2(x,fix(n),15);
end;
f0=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0).*f1./x+f0;
if(k <= nm)bi(k+1)=f; end;
f0=f1;
f1=f;
end;  k=0-1;
s0=bi0./f;
for  k=0:nm;
bi(k+1)=s0.*bi(k+1);
end;  k=fix(nm)+1;
end;
g0=bk0;
g1=bk1;
for  k=2:nm;
g=2.0d0.*(k-1.0d0)./x.*g1+g0;
bk(k+1)=g;
g0=g1;
g1=g;
end;  k=fix(nm)+1;
for  k=2:nm;
di(k+1)=bi(k-1+1)-k./x.*bi(k+1);
dk(k+1)=-bk(k-1+1)-k./x.*bk(k+1);
end;  k=fix(nm)+1;
return;
end
function [x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1]=ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1,varargin);
%       =========================================================
%       Purpose: Compute modified Bessel functions I0(x), I1(1),
%                K0(x)and K1(x), and their derivatives
%       Input :  x   --- Argument(x ע 0)
%       Output:  BI0 --- I0(x)
%                DI0 --- I0'(x)
%                BI1 --- I1(x)
%                DI1 --- I1'(x)
%                BK0 --- K0(x)
%                DK0 --- K0'(x)
%                BK1 --- K1(x)
%                DK1 --- K1'(x)
%       =========================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,8);
pi=3.141592653589793d0;
el=0.5772156649015329d0;
x2=x.*x;
if(x == 0.0d0);
bi0=1.0d0;
bi1=0.0d0;
bk0=1.0d+300;
bk1=1.0d+300;
di0=0.0d0;
di1=0.5d0;
dk0=-1.0d+300;
dk1=-1.0d+300;
return;
elseif(x <= 18.0);
bi0=1.0d0;
r=1.0d0;
for  k=1:50;
r=0.25d0.*r.*x2./(k.*k);
bi0=bi0+r;
if(abs(r./bi0)< 1.0d-15)break; end;
end;
bi1=1.0d0;
r=1.0d0;
for  k=1:50;
r=0.25d0.*r.*x2./(k.*(k+1));
bi1=bi1+r;
if(abs(r./bi1)< 1.0d-15)break; end;
end;
bi1=0.5d0.*x.*bi1;
else;
a(:)=[0.125d0,7.03125d-2,7.32421875d-2,1.1215209960938d-1,2.2710800170898d-1,5.7250142097473d-1,1.7277275025845d0,6.0740420012735d0,2.4380529699556d01,1.1001714026925d02,5.5133589612202d02,3.0380905109224d03];
b(:)=[-0.375d0,-1.171875d-1,-1.025390625d-1,-1.4419555664063d-1,-2.7757644653320d-1,-6.7659258842468d-1,-1.9935317337513d0,-6.8839142681099d0,-2.7248827311269d01,-1.2159789187654d02,-6.0384407670507d02,-3.3022722944809d03];
k0=12;
if(x >= 35.0)k0=9; end;
if(x >= 50.0)k0=7; end;
ca=exp(x)./sqrt(2.0d0.*pi.*x);
bi0=1.0d0;
xr=1.0d0./x;
for  k=1:k0;
bi0=bi0+a(k).*xr.^k;
end;  k=k0+1;
bi0=ca.*bi0;
bi1=1.0d0;
for  k=1:k0;
bi1=bi1+b(k).*xr.^k;
end;  k=k0+1;
bi1=ca.*bi1;
end;
if(x <= 9.0d0);
ct=-(log(x./2.0d0)+el);
bk0=0.0d0;
w0=0.0d0;
r=1.0d0;
for  k=1:50;
w0=w0+1.0d0./k;
r=0.25d0.*r./(k.*k).*x2;
bk0=bk0+r.*(w0+ct);
if(abs((bk0-ww)./bk0)< 1.0d-15)break; end;
ww=bk0;
end;
bk0=bk0+ct;
else;
a1(:)=[0.125d0,0.2109375d0,1.0986328125d0,1.1775970458984d01,2.1461706161499d02,5.9511522710323d03,2.3347645606175d05,1.2312234987631d07];
cb=0.5d0./x;
xr2=1.0d0./x2;
bk0=1.0d0;
for  k=1:8;
bk0=bk0+a1(k).*xr2.^k;
end;  k=8+1;
bk0=cb.*bk0./bi0;
end;
bk1=(1.0d0./x-bi1.*bk0)./bi0;
di0=bi1;
di1=bi0-bi1./x;
dk0=-bk1;
dk1=-bk0-bk1./x;
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

