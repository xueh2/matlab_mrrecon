function mikv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes modified Bessel functions
%                Iv(x)and Kv(x)with an arbitrary order, and
%                their derivatives using subroutine IKV
%       Input :  x --- Argument(x ע 0)
%                v --- Order of Iv(x)and Kv(x)
%(v = n+v0, 0 ף n ף 250 , 0 ף v0 < 1)
%       Output:  BI(n)--- In+v0(x)
%                DI(n)--- In+v0'(x)
%                BK(n)--- Kn+v0(x)
%                DK(n)--- Kn+v0'(x)
%       Example: v = n+v0,   v0 = .25,   x = 10.0
%     n       Iv(x)Iv'(x)Kv(x)Kv'(x)
%    ----------------------------------------------------------------
%     0   .28064359D+04  .26631677D+04  .17833184D-04 -.18709581D-04
%     1   .25930068D+04  .24823101D+04  .19155411D-04 -.20227611D-04
%     2   .21581842D+04  .21074153D+04  .22622037D-04 -.24245369D-04
%     3   .16218239D+04  .16310915D+04  .29335327D-04 -.32156018D-04
%     4   .11039987D+04  .11526244D+04  .41690000D-04 -.47053577D-04
%     5   .68342498D+03  .74520058D+03  .64771827D-04 -.75695209D-04
%       =============================================================
v=[];x=[];vm=[];
 global bi;
global di;
global bk;
global dk;
fprintf(1,'%s \n','  please enter v, x ');
%        READ(*,*)V,X
v=5.25;
x=10.0;
n=fix(v);
v0=v-n;
fprintf(1,[repmat(' ',1,8),'v = n+v0',',   ','v0 =','%7.5g',',   ','x =','%5.1g' ' \n'],v0,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[v,x,vm,bi,di,bk,dk]=ikv(v,x,vm,bi,di,bk,dk);
nm=fix(vm);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','    n       iv(x)iv''(x)');fprintf(1,'%s \n','kv(x)kv''(x)');
fprintf(1,'%s ','  -------------------------------------');fprintf(1,'%s \n','--------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,2),'%4g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,bi(k+1),di(k+1),bk(k+1),dk(k+1));
end;  k=nm+1;
%format(2x,i4,1x,4d16.8);
%format(8x,'v = n+v0',',   ',',f7.5,',   ',',f5.1);
end
function [v,x,vm,bi,di,bk,dk]=ikv(v,x,vm,bi,di,bk,dk,varargin);
%       =======================================================
%       Purpose: Compute modified Bessel functions Iv(x)and
%                Kv(x), and their derivatives
%       Input :  x --- Argument(x ע 0)
%                v --- Order of Iv(x)and Kv(x)
%(v = n+v0, n = 0,1,2,..., 0 ף v0 < 1)
%       Output:  BI(n)--- In+v0(x)
%                DI(n)--- In+v0'(x)
%                BK(n)--- Kn+v0(x)
%                DK(n)--- Kn+v0'(x)
%                VM --- Highest order computed
%       Routines called:
%(1)GAMMA for computing the gamma function
%(2)MSTA1 and MSTA2 to compute the starting
%                point for backward recurrence
%       =======================================================
v0p=[];gap=[];v0n=[];gan=[];
pi=3.141592653589793d0;
x2=x.*x;
n=fix(v);
v0=v-n;
if(n == 0)n=1; end;
if(x < 1.0d-100);
for  k=0:n;
bi(k+1)=0.0d0;
di(k+1)=0.0d0;
bk(k+1)=-1.0d+300;
dk(k+1)=1.0d+300;
end;  k=n+1;
if(v == 0.0);
bi(0+1)=1.0d0;
di(1+1)=0.5d0;
end;
vm=v;
return;
end;
piv=pi.*v0;
vt=4.0d0.*v0.*v0;
if(v0 == 0.0d0);
a1=1.0d0;
else;
v0p=1.0d0+v0;
[v0p,gap]=gamma(v0p,gap);
a1=(0.5d0.*x).^v0./gap;
end;
k0=14;
if(x >= 35.0)k0=10; end;
if(x >= 50.0)k0=8; end;
if(x <= 18.0);
bi0=1.0d0;
r=1.0d0;
for  k=1:30;
r=0.25d0.*r.*x2./(k.*(k+v0));
bi0=bi0+r;
if(abs(r./bi0)< 1.0d-15)break; end;
end;
bi0=bi0.*a1;
else;
ca=exp(x)./sqrt(2.0d0.*pi.*x);
sum=1.0d0;
r=1.0d0;
for  k=1:k0;
r=-0.125d0.*r.*(vt-(2.0d0.*k-1.0d0).^2.0)./(k.*x);
sum=sum+r;
end;  k=k0+1;
bi0=ca.*sum;
end;
m=msta1(x,200);
if(m < n);
n=m;
else;
m=msta2(x,n,15);
end;
f2=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(v0+k+1.0d0)./x.*f1+f2;
if(k <= n)bi(k+1)=f; end;
f2=f1;
f1=f;
end;  k=0-1;
cs=bi0./f;
for  k=0:n;
bi(k+1)=cs.*bi(k+1);
end;  k=n+1;
di(0+1)=v0./x.*bi(0+1)+bi(1+1);
for  k=1:n;
di(k+1)=-(k+v0)./x.*bi(k+1)+bi(k-1+1);
end;  k=n+1;
if(x <= 9.0d0);
if(v0 == 0.0d0);
ct=-log(0.5d0.*x)-0.5772156649015329d0;
cs=0.0d0;
w0=0.0d0;
r=1.0d0;
for  k=1:50;
w0=w0+1.0d0./k;
r=0.25d0.*r./(k.*k).*x2;
cs=cs+r.*(w0+ct);
wa=abs(cs);
if(abs((wa-ww)./wa)< 1.0d-15)break; end;
ww=wa;
end;
bk0=ct+cs;
else;
v0n=1.0d0-v0;
[v0n,gan]=gamma(v0n,gan);
a2=1.0d0./(gan.*(0.5d0.*x).^v0);
a1=(0.5d0.*x).^v0./gap;
sum=a2-a1;
r1=1.0d0;
r2=1.0d0;
for  k=1:120;
r1=0.25d0.*r1.*x2./(k.*(k-v0));
r2=0.25d0.*r2.*x2./(k.*(k+v0));
sum=sum+a2.*r1-a1.*r2;
wa=abs(sum);
if(abs((wa-ww)./wa)< 1.0d-15)break; end;
ww=wa;
end;
bk0=0.5d0.*pi.*sum./sin(piv);
end;
else;
cb=exp(-x).*sqrt(0.5d0.*pi./x);
sum=1.0d0;
r=1.0d0;
for  k=1:k0;
r=0.125d0.*r.*(vt-(2.0.*k-1.0).^2.0)./(k.*x);
sum=sum+r;
end;  k=k0+1;
bk0=cb.*sum;
end;
bk1=(1.0d0./x-bi(1+1).*bk0)./bi(0+1);
bk(0+1)=bk0;
bk(1+1)=bk1;
for  k=2:n;
bk2=2.0d0.*(v0+k-1.0d0)./x.*bk1+bk0;
bk(k+1)=bk2;
bk0=bk1;
bk1=bk2;
end;  k=n+1;
dk(0+1)=v0./x.*bk(0+1)-bk(1+1);
for  k=1:n;
dk(k+1)=-(k+v0)./x.*bk(k+1)-bk(k-1+1);
end;  k=n+1;
vm=n+v0;
return;
end
function [x,ga]=gamma(x,ga,varargin);
%       ==================================================
%       Purpose: Compute gamma function ג(x)
%       Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%       Output:  GA --- ג(x)
%       ==================================================
 g=zeros(1,26);
pi=3.141592653589793d0;
if(x == fix(x));
if(x > 0.0d0);
ga=1.0d0;
m1=x-1;
for  k=2:m1;
ga=ga.*k;
end;  k=m1+1;
else;
ga=1.0d+300;
end;
else;
if(abs(x)> 1.0d0);
z=abs(x);
m=fix(z);
r=1.0d0;
for  k=1:m;
r=r.*(z-k);
end;  k=m+1;
z=z-m;
else;
z=x;
end;
g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
gr=g(26);
for  k=25:-1:1;
gr=gr.*z+g(k);
end;  k=1-1;
ga=1.0d0./(gr.*z);
if(abs(x)> 1.0d0);
ga=ga.*r;
if(x < 0.0d0)ga=-pi./(x.*ga.*sin(pi.*x)); end;
end;
end;
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

