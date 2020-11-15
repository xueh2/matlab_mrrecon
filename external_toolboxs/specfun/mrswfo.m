function mrswfo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes the radial oblate spheriodal
%       functions of the first and second kinds, and their
%       derivatives using subroutine RSWFO
%       Input :  m  --- Mode parameter, m = 0,1,2,...
%       n  --- Mode parameter, n = m,m+1,m+2,...
%       c  --- Spheroidal parameter
%       x  --- Argument(x ע 0)
%       cv --- Characteristic value
%       KF --- Function code
%       KF=1 for the first kind
%       KF=2 for the second kind
%       KF=3 for both the first and second kinds
%       Output:  R1F --- Radial function of the first kind
%       R1D --- Derivative of the radial function of
%       the first kind
%       R2F --- Radial function of the second kind
%       R2D --- Derivative of the radial function of
%       the second kind
%       Example:
%       KD=-1, m = 2, n = 3, c = 5.0 and cv = 2.10980581604
%       x  R23(1,-ic,ix)R23(1)'(-ic,ix)R23(2,-ic,ix)R23(2)'(-ic,ix)
%       ----------------------------------------------------------------
%       --
%       0.0      0.0000000(-1)4.9911346(-1)-4.0071049(-1)3
%       .3065682
%       0.5(-1)2.0215966(-1)1.9559661(-1)-1.5154835(-1)6
%       .4482526
%       1.0(-1)1.0526737(-1)-5.4010624(-1)1.3065731(-1)2
%       .7958491
%       1.5(-1)-1.1611246(-2)-8.9414690(-2)3.7405936(-1)-5
%       .0118500
%       5.0(-2)3.6463251(-2)-8.3283065(-2)1.5549615(-1)1
%       .7544481
%       ==============================================================
m=[];n=[];c=[];cv=[];eg=[];x=[];kf=[];r1f=[];r1d=[];r2f=[];r2d=[];
 eg=zeros(1,200);
fprintf(1,'%s \n','please enter kf');
%       READ(*,*)KF
kf=3;
fprintf(1,[repmat(' ',1,1),'kf=','%3g' ' \n'],kf);
fprintf(1,'%s \n','please enter m, n, c and x(x ע 0)');
%       READ(*,*)M,N,C,X
m=2;
n=3;
c=5.0;
x=5.0;
[m,n,c,dumvar4,cv,eg]=segv(m,n,c,-1,cv,eg);
fprintf(1,[repmat(' ',1,1),'m=','%2g',',   ','n=','%2g',',   ','c=','%5.1g',',    ','cv =','%18.10g',',  ','x=','%5.2g' ' \n'],m,n,c,cv,x);
fprintf(1,'%0.15g \n');
[m,n,c,x,cv,kf,r1f,r1d,r2f,r2d]=rswfo(m,n,c,x,cv,kf,r1f,r1d,r2f,r2d);
fprintf(1,'%s ','   x    rmn(1,-ic,ix)rmn''(1,-ic,ix)');fprintf(1,'%s \n',' rmn(2,-ic,ix)rmn''(2,-ic,ix)');
fprintf(1,'%s ','  ---------------------------------------------');fprintf(1,'%s \n','---------------------------');
if(kf == 1);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%17.8g',1,4) ' \n'],x,r1f,r1d);
elseif(kf == 2);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,34),repmat('%17.8g',1,4) ' \n'],x,r2f,r2d);
elseif(kf == 3);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%17.8g',1,4) ' \n'],x,r1f,r1d,r2f,r2d);
end;
if(kf == 3);
fprintf(1,[repmat(' ',1,1), '\n ' ,repmat(' ',1,1),'wronskian check:', '\n ' ,repmat(' ',1,1),'computed value =','%17.8g',repmat(' ',1,5),'exact value =','%17.8g' ' \n'],r1f.*r2d-r2f.*r1d,1.0d0./(c.*(x.*x+1.0d0)));
fprintf(1,[repmat(' ',1,1), '\n ' ,repmat(' ',1,1),'caution: this check is not accurate if it ','involves', '\n ' ,repmat(' ',1,1),'         the subtraction of two ','similar numbers' ' \n']);
end;
%format(1x,,i3);
%format(1x,,i2,',   ',,i2,',   ',,f5.1,',    ',,f18.10,',  ',,f5.2);
%format(1x,f5.1,4d17.8);
%format(1x,f5.1,34x,4d17.8);
%format(1x,./1x,['wronskian check:'],./1x,',d17.8,5x,',d17.8);
%format(1x,./1x,['caution: this check is not accurate if it '],'involves',./1x,'         the subtraction of two ', 'similar numbers');
end
function [m,n,c,x,cv,kf,r1f,r1d,r2f,r2d]=rswfo(m,n,c,x,cv,kf,r1f,r1d,r2f,r2d,varargin);
%       ==========================================================
%       Purpose: Compute oblate radial functions of the first
%       and second kinds, and their derivatives
%       Input :  m  --- Mode parameter,  m = 0,1,2,...
%       n  --- Mode parameter,  n = m,m+1,m+2,...
%       c  --- Spheroidal parameter
%       x  --- Argument(x ע 0)
%       cv --- Characteristic value
%       KF --- Function code
%       KF=1 for the first kind
%       KF=2 for the second kind
%       KF=3 for both the first and second kinds
%       Output:  R1F --- Radial function of the first kind
%       R1D --- Derivative of the radial function of
%       the first kind
%       R2F --- Radial function of the second kind
%       R2D --- Derivative of the radial function of
%       the second kind
%       Routines called:
%(1)SDMN for computing expansion coefficients dk
%(2)RMN1 for computing prolate or oblate radial
%       function of the first kind
%(3)RMN2L for computing prolate or oblate radial
%       function of the second kind for a large argument
%(4)RMN2SO for computing oblate radial functions of
%       the second kind for a small argument
%       ==========================================================
kd=[];df=[];id=[];
 df=zeros(1,200);
kd=-1;
[m,n,c,cv,kd,df]=sdmn(fix(m),fix(n),c,cv,kd,df);
if(kf ~= 2);
[m,n,c,x,df,kd,r1f,r1d]=rmn1(fix(m),fix(n),c,x,df,kd,r1f,r1d);
end;
if(kf > 1);
id=10;
if(x > 1.0d-8);
[m,n,c,x,df,kd,r2f,r2d,id]=rmn2l(fix(m),fix(n),c,x,df,kd,r2f,r2d,id);
end;
if(id > -1);
[m,n,c,x,cv,df,kd,r2f,r2d]=rmn2so(fix(m),fix(n),c,x,cv,df,kd,r2f,r2d);
end;
end;
return;
end
function [m,n,c,cv,kd,df]=sdmn(m,n,c,cv,kd,df,varargin);
%       =====================================================
%       Purpose: Compute the expansion coefficients of the
%       prolate and oblate spheroidal functions, dk
%       Input :  m  --- Mode parameter
%       n  --- Mode parameter
%       c  --- Spheroidal parameter
%       cv --- Characteristic value
%       KD --- Function code
%       KD=1 for prolate; KD=-1 for oblate
%       Output:  DF(k)--- Expansion coefficients dk;
%       DF(1), DF(2), ... correspond to
%       d0, d2, ... for even n-m and d1,
%       d3, ... for odd n-m
%       =====================================================
 a=zeros(1,200);
d=zeros(1,200);
g=zeros(1,200);
fl=0.0;
sw=0.0;
nm=25+fix(0.5.*(fix(n)-fix(m))+c);
if(c < 1.0d-10);
for  i=1:nm;
df(i)=0d0;
end;  i=nm+1;
df((n-m)./2+1)=1.0d0;
return;
end;
cs=c.*c.*fix(kd);
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
for  i=1:nm+2;
if(ip == 0)k=2.*(i-1); end;
if(ip == 1)k=2.*i-1; end;
dk0=fix(m)+k;
dk1=fix(m)+k+1;
dk2=2.*(fix(m)+k);
d2k=2.*fix(m)+k;
a(i)=(d2k+2.0).*(d2k+1.0)./((dk2+3.0).*(dk2+5.0)).*cs;
d(i)=dk0.*dk1+(2.0.*dk0.*dk1-2.0.*fix(m).*fix(m)-1.0)./((dk2-1.0).*(dk2+3.0)).*cs;
g(i)=k.*(k-1.0)./((dk2-3.0).*(dk2-1.0)).*cs;
end;  i=nm+2+1;
fs=1.0d0;
f1=0.0d0;
f0=1.0d-100;
kb=0;
df(nm+1)=0.0d0;
for  k=nm:-1:1;
f=-((d(k+1)-cv).*f0+a(k+1).*f1)./g(k+1);
if(abs(f)> abs(df(k+1)));
df(k)=f;
f1=f0;
f0=f;
if(abs(f)> 1.0d+100);
for  k1=k:nm;
df(k1)=df(k1).*1.0d-100;
end;  k1=nm+1;
f1=f1.*1.0d-100;
f0=f0.*1.0d-100;
end;
else;
kb=k;
fl=df(k+1);
f1=1.0d-100;
f2=-(d(1)-cv)./a(1).*f1;
df(1)=f1;
if(kb == 1);
fs=f2;
elseif(kb == 2);
df(2)=f2;
fs=-((d(2)-cv).*f2+g(2).*f1)./a(2);
else;
df(2)=f2;
for  j=3:kb+1;
f=-((d(j-1)-cv).*f2+g(j-1).*f1)./a(j-1);
if(j <= kb)df(j)=f; end;
if(abs(f)> 1.0d+100);
for  k1=1:j;
df(k1)=df(k1).*1.0d-100;
end;  k1=j+1;
f=f.*1.0d-100;
f2=f2.*1.0d-100;
end;
f1=f2;
f2=f;
end;  j=kb+1+1;
fs=f;
end;
break;
end;
end;
su1=0.0d0;
r1=1.0d0;
for  j=m+ip+1:2.*(m+ip);
r1=r1.*j;
end;  j=2.*(fix(m)+ip)+1;
su1=df(1).*r1;
for  k=2:kb;
r1=-r1.*(k+fix(m)+ip-1.5d0)./(k-1.0d0);
su1=su1+r1.*df(k);
end;  k=kb+1;
su2=0.0d0;
for  k=kb+1:nm;
if(k ~= 1)r1=-r1.*(k+m+ip-1.5d0)./(k-1.0d0); end;
su2=su2+r1.*df(k);
if(abs(sw-su2)< abs(su2).*1.0d-14)break; end;
sw=su2;
end;
r3=1.0d0;
for  j=1:(m+n+ip)./2;
r3=r3.*(j+0.5d0.*(fix(n)+fix(m)+ip));
end;  j=(fix(m)+fix(n)+ip)./2+1;
r4=1.0d0;
for  j=1:(n-m-ip)./2;
r4=-4.0d0.*r4.*j;
end;  j=(fix(n)-fix(m)-ip)./2+1;
s0=r3./(fl.*(su1./fs)+su2)./r4;
for  k=1:kb;
df(k)=fl./fs.*s0.*df(k);
end;  k=kb+1;
for  k=kb+1:nm;
df(k)=s0.*df(k);
end;  k=nm+1;
return;
end
function [m,n,c,x,df,kd,r1f,r1d]=rmn1(m,n,c,x,df,kd,r1f,r1d,varargin);
%       =======================================================
%       Purpose: Compute prolate and oblate spheroidal radial
%       functions of the first kind for given m, n,
%       c and x
%       Routines called:
%(1)SCKB for computing expansion coefficients c2k
%(2)SPHJ for computing the spherical Bessel
%       functions of the first kind
%       =======================================================
ck=[];nm2=[];cx=[];sj=[];dj=[];
 ck=zeros(1,200);
 sj=zeros(1,251+1);
dj=zeros(1,251+1);
sw=0.0;
eps=1.0d-14;
ip=1;
nm1=fix((fix(n)-fix(m))./2);
if(n-m == 2.*nm1)ip=0; end;
nm=25+nm1+fix(c);
reg=1.0d0;
if(m+nm > 80)reg=1.0d-200; end;
r0=reg;
for  j=1:2.*m+ip;
r0=r0.*j;
end;  j=2.*fix(m)+ip+1;
r=r0;
suc=r.*df(1);
for  k=2:nm;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
suc=suc+r.*df(k);
if(k > nm1&abs(suc-sw)< abs(suc).*eps)break; end;
sw=suc;
end;
if(x == 0.0);
[m,n,c,df,ck]=sckb(fix(m),fix(n),c,df,ck);
sum=0.0d0;
for  j=1:nm;
sum=sum+ck(j);
if(abs(sum-sw1)< abs(sum).*eps)break; end;
sw1=sum;
end;
r1=1.0d0;
for  j=1:(n+m+ip)./2;
r1=r1.*(j+0.5d0.*(fix(n)+fix(m)+ip));
end;  j=(fix(n)+fix(m)+ip)./2+1;
r2=1.0d0;
for  j=1:m;
r2=2.0d0.*c.*r2.*j;
end;  j=fix(m)+1;
r3=1.0d0;
for  j=1:(n-m-ip)./2;
r3=r3.*j;
end;  j=(fix(n)-fix(m)-ip)./2+1;
sa0=(2.0.*(fix(m)+ip)+1.0).*r1./(2.0.^fix(n).*c.^ip.*r2.*r3);
if(ip == 0);
r1f=sum./(sa0.*suc).*df(1).*reg;
r1d=0.0d0;
elseif(ip == 1);
r1f=0.0d0;
r1d=sum./(sa0.*suc).*df(1).*reg;
end;
return;
end;
cx=c.*x;
nm2=2.*nm+fix(m);
[nm2,cx,nm2,sj,dj]=sphj(nm2,cx,nm2,sj,dj);
a0=(1.0d0-fix(kd)./(x.*x)).^(0.5d0.*fix(m))./suc;
r1f=0.0d0;
for  k=1:nm;
l=2.*k+fix(m)-fix(n)-2+ip;
if(l == 4.*fix(l./4))lg=1; end;
if(l ~= 4.*fix(l./4))lg=-1; end;
if(k == 1);
r=r0;
else;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
end;
np=fix(m)+2.*k-2+ip;
r1f=r1f+lg.*r.*df(k).*sj(np+1);
if(k > nm1&abs(r1f-sw)< abs(r1f).*eps)break; end;
sw=r1f;
end;
r1f=r1f.*a0;
b0=fix(kd).*fix(m)./x.^3.0d0./(1.0-fix(kd)./(x.*x)).*r1f;
sud=0.0d0;
for  k=1:nm;
l=2.*k+fix(m)-fix(n)-2+ip;
if(l == 4.*fix(l./4))lg=1; end;
if(l ~= 4.*fix(l./4))lg=-1; end;
if(k == 1);
r=r0;
else;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
end;
np=fix(m)+2.*k-2+ip;
sud=sud+lg.*r.*df(k).*dj(np+1);
if(k > nm1&abs(sud-sw)< abs(sud).*eps)break; end;
sw=sud;
end;
r1d=b0+a0.*c.*sud;
return;
end
function [m,n,c,df,ck]=sckb(m,n,c,df,ck,varargin);
%       ======================================================
%       Purpose: Compute the expansion coefficients of the
%       prolate and oblate spheroidal functions, c2k
%       Input :  m  --- Mode parameter
%       n  --- Mode parameter
%       c  --- Spheroidal parameter
%       DF(k)--- Expansion coefficients dk
%       Output:  CK(k)--- Expansion coefficients ck;
%       CK(1), CK(2), ... correspond to
%       c0, c2, ...
%       ======================================================
if(c <= 1.0d-10)c=1.0d-10; end;
nm=25+fix(0.5.*(fix(n)-fix(m))+c);
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
reg=1.0d0;
if(m+nm > 80)reg=1.0d-200; end;
fac=-0.5d0.^fix(m);
for  k=0:nm-1;
fac=-fac;
i1=2.*k+ip+1;
r=reg;
for  i=i1:i1+2.*m-1;
r=r.*i;
end;  i=i1+2.*fix(m)-1+1;
i2=k+fix(m)+ip;
for  i=i2:i2+k-1;
r=r.*(i+0.5d0);
end;  i=i2+k-1+1;
sum=r.*df(k+1);
for  i=k+1:nm;
d1=2.0d0.*i+ip;
d2=2.0d0.*fix(m)+d1;
d3=i+fix(m)+ip-0.5d0;
r=r.*d2.*(d2-1.0d0).*i.*(d3+k)./(d1.*(d1-1.0d0).*(i-k).*d3);
sum=sum+r.*df(i+1);
if(abs(sw-sum)< abs(sum).*1.0d-14)break; end;
sw=sum;
end;
r1=reg;
for  i=2:m+k;
r1=r1.*i;
end;  i=fix(m)+k+1;
ck(k+1)=fac.*sum./r1;
end;
return;
end
function [n,x,nm,sj,dj]=sphj(n,x,nm,sj,dj,varargin);
%       =======================================================
%       Purpose: Compute spherical Bessel functions jn(x)and
%       their derivatives
%       Input :  x --- Argument of jn(x)
%       n --- Order of jn(x,n = 0,1,תתת)
%       Output:  SJ(n)--- jn(x)
%       DJ(n)--- jn'(x)
%       NM --- Highest order computed
%       Routines called:
%       MSTA1 and MSTA2 for computing the starting
%       point for backward recurrence
%       =======================================================
nm=fix(fix(n));
if(abs(x)== 1.0d-100);
for  k=0:n;
sj(k+1)=0.0d0;
dj(k+1)=0.0d0;
end;  k=fix(n)+1;
sj(0+1)=1.0d0;
dj(1+1)=.3333333333333333d0;
return;
end;
sj(0+1)=sin(x)./x;
sj(1+1)=(sj(0+1)-cos(x))./x;
if(n >= 2);
sa=sj(0+1);
sb=sj(1+1);
m=msta1(x,200);
if(m < n);
nm=fix(m);
else;
m=msta2(x,fix(n),15);
end;
f0=0.0d0;
f1=1.0d0-100;
for  k=m:-1:0;
f=(2.0d0.*k+3.0d0).*f1./x-f0;
if(k <= nm)sj(k+1)=f; end;
f0=f1;
f1=f;
end;  k=0-1;
if(abs(sa)> abs(sb))cs=sa./f; end;
if(abs(sa)<= abs(sb))cs=sb./f0; end;
for  k=0:nm;
sj(k+1)=cs.*sj(k+1);
end;  k=fix(nm)+1;
end;
dj(0+1)=(cos(x)-sin(x)./x)./x;
for  k=1:nm;
dj(k+1)=sj(k-1+1)-(k+1.0d0).*sj(k+1)./x;
end;  k=fix(nm)+1;
return;
end
function [msta1Result]=msta1(x,mp,varargin);
%       ===================================================
%       Purpose: Determine the starting point for backward
%       recurrence such that the magnitude of
%       Jn(x)at that point is about 10^(-MP)
%       Input :  x     --- Argument of Jn(x)
%       MP    --- Value of magnitude
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
%       recurrence such that all Jn(x)has MP
%       significant digits
%       Input :  x  --- Argument of Jn(x)
%       n  --- Order of Jn(x)
%       MP --- Significant digit
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
function [m,n,c,x,df,kd,r2f,r2d,id]=rmn2l(m,n,c,x,df,kd,r2f,r2d,id,varargin);
%       ========================================================
%       Purpose: Compute prolate and oblate spheroidal radial
%       functions of the second kind for given m, n,
%       c and a large cx
%       Routine called:
%       SPHY for computing the spherical Bessel
%       functions of the second kind
%       ========================================================
nm2=[];cx=[];sy=[];dy=[];
  sy=zeros(1,251+1);
dy=zeros(1,251+1);
sw=0.0;
eps=1.0d-14;
ip=1;
nm1=fix((fix(n)-fix(m))./2);
if(n-m == 2.*nm1)ip=0; end;
nm=25+nm1+fix(c);
reg=1.0d0;
if(m+nm > 80)reg=1.0d-200; end;
nm2=2.*nm+fix(m);
cx=c.*x;
[nm2,cx,nm2,sy,dy]=sphy(nm2,cx,nm2,sy,dy);
r0=reg;
for  j=1:2.*m+ip;
r0=r0.*j;
end;  j=2.*fix(m)+ip+1;
r=r0;
suc=r.*df(1);
for  k=2:nm;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
suc=suc+r.*df(k);
if(k > nm1&abs(suc-sw)< abs(suc).*eps)break; end;
sw=suc;
end;
a0=(1.0d0-fix(kd)./(x.*x)).^(0.5d0.*fix(m))./suc;
r2f=0.0;
for  k=1:nm;
l=2.*k+fix(m)-fix(n)-2+ip;
if(l == 4.*fix(l./4))lg=1; end;
if(l ~= 4.*fix(l./4))lg=-1; end;
if(k == 1);
r=r0;
else;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
end;
np=fix(m)+2.*k-2+ip;
r2f=r2f+lg.*r.*(df(k).*sy(np+1));
eps1=abs(r2f-sw);
if(k > nm1&eps1 < abs(r2f).*eps)break; end;
sw=r2f;
end;
id1=fix(log10(eps1./abs(r2f)+eps));
r2f=r2f.*a0;
if(np >= nm2);
id=10;
return;
end;
b0=fix(kd).*fix(m)./x.^3.0d0./(1.0-fix(kd)./(x.*x)).*r2f;
sud=0.0d0;
for  k=1:nm;
l=2.*k+fix(m)-fix(n)-2+ip;
if(l == 4.*fix(l./4))lg=1; end;
if(l ~= 4.*fix(l./4))lg=-1; end;
if(k == 1);
r=r0;
else;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
end;
np=fix(m)+2.*k-2+ip;
sud=sud+lg.*r.*(df(k).*dy(np+1));
eps2=abs(sud-sw);
if(k > nm1&eps2 < abs(sud).*eps)break; end;
sw=sud;
end;
r2d=b0+a0.*c.*sud;
id2=fix(log10(eps2./abs(sud)+eps));
id=fix(max([id1,id2]));
return;
end
function [n,x,nm,sy,dy]=sphy(n,x,nm,sy,dy,varargin);
%       ======================================================
%       Purpose: Compute spherical Bessel functions yn(x)and
%       their derivatives
%       Input :  x --- Argument of yn(x,x ע 0)
%       n --- Order of yn(x,n = 0,1,תתת)
%       Output:  SY(n)--- yn(x)
%       DY(n)--- yn'(x)
%       NM --- Highest order computed
%       ======================================================
nm=fix(fix(n));
if(x < 1.0d-60);
for  k=0:n;
sy(k+1)=-1.0d+300;
dy(k+1)=1.0d+300;
end;  k=fix(n)+1;
return;
end;
sy(0+1)=-cos(x)./x;
sy(1+1)=(sy(0+1)-sin(x))./x;
f0=sy(0+1);
f1=sy(1+1);
for  k=2:n;
f=(2.0d0.*k-1.0d0).*f1./x-f0;
sy(k+1)=f;
if(abs(f)>= 1.0d+300)break; end;
f0=f1;
f1=f;
end;
nm=fix(k-1);
dy(0+1)=(sin(x)+cos(x)./x)./x;
for  k=1:nm;
dy(k+1)=sy(k-1+1)-(k+1.0d0).*sy(k+1)./x;
end;  k=fix(nm)+1;
return;
end
function [m,n,c,x,cv,df,kd,r2f,r2d]=rmn2so(m,n,c,x,cv,df,kd,r2f,r2d,varargin);
%       =============================================================
%       Purpose: Compute oblate radial functions of the second kind
%       with a small argument, Rmn(-ic,ix)& Rmn'(-ic,ix)
%       Routines called:
%(1)SCKB for computing the expansion coefficients c2k
%(2)KMN for computing the joining factors
%(3)QSTAR for computing the factor defined in(15.7.3)
%(4)CBK for computing the the expansion coefficient
%       defined in(15.7.6)
%(5)GMN for computing the function defined in(15.7.4)
%(6)RMN1 for computing the radial function of the first
%       kind
%       =============================================================
ck=[];dn=[];ck1=[];ck2=[];qs=[];qt=[];bk=[];gf=[];gd=[];r1f=[];r1d=[];
 bk=zeros(1,200);
ck=zeros(1,200);
 dn=zeros(1,200);
if(abs(df(1))<= 1.0d-280);
r2f=1.0d+300;
r2d=1.0d+300;
return;
end;
eps=1.0d-14;
pi=3.141592653589793d0;
nm=25+fix((fix(n)-fix(m))./2+c);
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
[m,n,c,df,ck]=sckb(fix(m),fix(n),c,df,ck);
[m,n,c,cv,kd,df,dn,ck1,ck2]=kmn(fix(m),fix(n),c,cv,fix(kd),df,dn,ck1,ck2);
[m,n,c,ck,ck1,qs,qt]=qstar(fix(m),fix(n),c,ck,ck1,qs,qt);
[m,n,c,cv,qt,ck,bk]=cbk(fix(m),fix(n),c,cv,qt,ck,bk);
if(x == 0.0d0);
sum=0.0d0;
for  j=1:nm;
sum=sum+ck(j);
if(abs(sum-sw)< abs(sum).*eps)break; end;
sw=sum;
end;
if(ip == 0);
r1f=sum./ck1;
r2f=-0.5d0.*pi.*qs.*r1f;
r2d=qs.*r1f+bk(1);
elseif(ip == 1);
r1d=sum./ck1;
r2f=bk(1);
r2d=-0.5d0.*pi.*qs.*r1d;
end;
return;
else;
[m,n,c,x,bk,gf,gd]=gmn(fix(m),fix(n),c,x,bk,gf,gd);
[m,n,c,x,df,kd,r1f,r1d]=rmn1(fix(m),fix(n),c,x,df,fix(kd),r1f,r1d);
h0=atan(x)-0.5d0.*pi;
r2f=qs.*r1f.*h0+gf;
r2d=qs.*(r1d.*h0+r1f./(1.0d0+x.*x))+gd;
end;
return;
end
function [m,n,c,ck,ck1,qs,qt]=qstar(m,n,c,ck,ck1,qs,qt,varargin);
%       =========================================================
%       Purpose: Compute Q*mn(-ic)for oblate radial functions
%       with a small argument
%       =========================================================
 ap=zeros(1,200);
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
r=1.0d0./ck(1).^2;
ap(1)=r;
for  i=1:m;
s=0.0d0;
for  l=1:i;
sk=0.0d0;
for  k=0:l;
sk=sk+ck(k+1).*ck(l-k+1);
end;  k=l+1;
s=s+sk.*ap(i-l+1);
end;  l=i+1;
ap(i+1)=-r.*s;
end;  i=fix(m)+1;
qs0=ap(fix(m)+1);
for  l=1:m;
r=1.0d0;
for  k=1:l;
r=r.*(2.0d0.*k+ip).*(2.0d0.*k-1.0d0+ip)./(2.0d0.*k).^2;
end;  k=l+1;
qs0=qs0+ap(fix(m)-l+1).*r;
end;  l=fix(m)+1;
qs=(-1).^ip.*ck1.*(ck1.*qs0)./c;
qt=-2.0d0./ck1.*qs;
return;
end
function [m,n,c,cv,qt,ck,bk]=cbk(m,n,c,cv,qt,ck,bk,varargin);
%       =====================================================
%       Purpose: Compute coefficient Bk's for oblate radial
%       functions with a small argument
%       =====================================================
  u=zeros(1,200);
v=zeros(1,200);
w=zeros(1,200);
eps=1.0d-14;
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
nm=25+fix(0.5.*(fix(n)-fix(m))+c);
u(1)=0.0d0;
n2=nm-2;
for  j=2:n2;
u(j)=c.*c;
end;  j=n2+1;
for  j=1:n2;
v(j)=(2.0.*j-1.0-ip).*(2.0.*(j-fix(m))-ip)+fix(m).*(fix(m)-1.0)-cv;
end;  j=n2+1;
for  j=1:nm-1;
w(j)=(2.0.*j-ip).*(2.0.*j+1.0-ip);
end;  j=nm-1+1;
if(ip == 0);
for  k=0:n2-1;
s1=0.0d0;
i1=k-fix(m)+1;
for  i=i1:nm;
if(~(i < 0));
r1=1.0d0;
for  j=1:k;
r1=r1.*(i+fix(m)-j)./j;
end;  j=k+1;
s1=s1+ck(i+1).*(2.0.*i+fix(m)).*r1;
if(abs(s1-sw)< abs(s1).*eps)break; end;
sw=s1;
end;
end;
bk(k+1)=qt.*s1;
end;
elseif(ip == 1);
for  k=0:n2-1;
s1=0.0d0;
i1=k-fix(m)+1;
for  i=i1:nm;
if(~(i < 0));
r1=1.0d0;
for  j=1:k;
r1=r1.*(i+fix(m)-j)./j;
end;  j=k+1;
if(i > 0)s1=s1+ck(i).*(2.0.*i+m-1).*r1; end;
s1=s1-ck(i+1).*(2.0.*i+fix(m)).*r1;
if(abs(s1-sw)< abs(s1).*eps)break; end;
sw=s1;
end;
end;
bk(k+1)=qt.*s1;
end;
end;
w(1)=w(1)./v(1);
bk(1)=bk(1)./v(1);
for  k=2:n2;
t=v(k)-w(k-1).*u(k);
w(k)=w(k)./t;
bk(k)=(bk(k)-bk(k-1).*u(k))./t;
end;  k=n2+1;
for  k=n2-1:-1:1;
bk(k)=bk(k)-w(k).*bk(k+1);
end;  k=1-1;
return;
end
function [m,n,c,x,bk,gf,gd]=gmn(m,n,c,x,bk,gf,gd,varargin);
%       ===========================================================
%       Purpose: Compute gmn(-ic,ix)and its derivative for oblate
%       radial functions with a small argument
%       ===========================================================
eps=1.0d-14;
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
nm=25+fix(0.5.*(fix(n)-fix(m))+c);
xm=(1.0d0+x.*x).^(-0.5d0.*fix(m));
gf0=0.0d0;
for  k=1:nm;
gf0=gf0+bk(k).*x.^(2.0.*k-2.0);
if(abs((gf0-gw)./gf0)< eps&k >= 10)break; end;
gw=gf0;
end;
gf=xm.*gf0.*x.^(1-ip);
gd1=-fix(m).*x./(1.0d0+x.*x).*gf;
gd0=0.0d0;
for  k=1:nm;
if(ip == 0);
gd0=gd0+(2.0d0.*k-1.0).*bk(k).*x.^(2.0.*k-2.0);
else;
gd0=gd0+2.0d0.*k.*bk(k+1).*x.^(2.0.*k-1.0);
end;
if(abs((gd0-gw)./gd0)< eps&k >= 10)break; end;
gw=gd0;
end;
gd=gd1+xm.*gd0;
return;
end
function [m,n,c,cv,kd,df,dn,ck1,ck2]=kmn(m,n,c,cv,kd,df,dn,ck1,ck2,varargin);
%       ===================================================
%       Purpose: Compute the expansion coefficients of the
%       prolate and oblate spheroidal functions
%       and joining factors
%       ===================================================
 u=zeros(1,200);
v=zeros(1,200);
w=zeros(1,200);
 tp=zeros(1,200);
rk=zeros(1,200);
nm=25+fix(0.5.*(fix(n)-fix(m))+c);
nn=nm+fix(m);
cs=c.*c.*fix(kd);
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
for  i=1:nn+3;
if(ip == 0)k=-2.*(i-1); end;
if(ip == 1)k=-(2.*i-3); end;
gk0=2.0d0.*fix(m)+k;
gk1=(fix(m)+k).*(fix(m)+k+1.0d0);
gk2=2.0d0.*(fix(m)+k)-1.0d0;
gk3=2.0d0.*(fix(m)+k)+3.0d0;
u(i)=gk0.*(gk0-1.0d0).*cs./(gk2.*(gk2+2.0d0));
v(i)=gk1-cv+(2.0d0.*(gk1-fix(m).*fix(m))-1.0d0).*cs./(gk2.*gk3);
w(i)=(k+1.0d0).*(k+2.0d0).*cs./((gk2+2.0d0).*gk3);
end;  i=nn+3+1;
for  k=1:m;
t=v(fix(m)+1);
for  l=0:m-k-1;
t=v(fix(m)-l)-w(fix(m)-l+1).*u(fix(m)-l)./t;
end;  l=fix(m)-k-1+1;
rk(k)=-u(k)./t;
end;  k=fix(m)+1;
r=1.0d0;
for  k=1:m;
r=r.*rk(k);
dn(k)=df(1).*r;
end;  k=fix(m)+1;
tp(nn)=v(nn+1);
for  k=nn-1:-1:m+1;
tp(k)=v(k+1)-w(k+2).*u(k+1)./tp(k+1);
if(k > m+1)rk(k)=-u(k)./tp(k); end;
end;  k=fix(m)+1-1;
if(m == 0)dnp=df(1); end;
if(m ~= 0)dnp=dn(m); end;
dn(m+1)=(-1).^ip.*dnp.*cs./((2.0.*fix(m)-1.0).*(2.0.*fix(m)+1.0-4.0.*ip).*tp(fix(m)+1));
for  k=m+2:nn;
dn(k)=rk(k).*dn(k-1);
end;  k=nn+1;
r1=1.0d0;
for  j=1:(n+m+ip)./2;
r1=r1.*(j+0.5d0.*(fix(n)+fix(m)+ip));
end;  j=(fix(n)+fix(m)+ip)./2+1;
nm1=(fix(n)-fix(m))./2;
r=1.0d0;
for  j=1:2.*m+ip;
r=r.*j;
end;  j=2.*fix(m)+ip+1;
su0=r.*df(1);
for  k=2:nm;
r=r.*(fix(m)+k-1.0).*(fix(m)+k+ip-1.5d0)./(k-1.0d0)./(k+ip-1.5d0);
su0=su0+r.*df(k);
if(k > nm1&abs((su0-sw)./su0)< 1.0d-14)break; end;
sw=su0;
end;
if(~(kd == 1));
r2=1.0d0;
for  j=1:m;
r2=2.0d0.*c.*r2.*j;
end;  j=fix(m)+1;
r3=1.0d0;
for  j=1:(n-m-ip)./2;
r3=r3.*j;
end;  j=(fix(n)-fix(m)-ip)./2+1;
sa0=(2.0.*(fix(m)+ip)+1.0).*r1./(2.0.^fix(n).*c.^ip.*r2.*r3.*df(1));
ck1=sa0.*su0;
if(kd == -1)return; end;
end;
r4=1.0d0;
for  j=1:(n-m-ip)./2;
r4=4.0d0.*r4.*j;
end;  j=(fix(n)-fix(m)-ip)./2+1;
r5=1.0d0;
for  j=1:m;
r5=r5.*(j+fix(m))./c;
end;  j=fix(m)+1;
g0=dn(fix(m));
if(m == 0)g0=df(1); end;
sb0=(ip+1.0).*c.^(ip+1)./(2.0.*ip.*(fix(m)-2.0)+1.0)./(2.0.*fix(m)-1.0);
ck2=(-1).^ip.*sb0.*r4.*r5.*g0./r1.*su0;
return;
end
function [m,n,c,kd,cv,eg]=segv(m,n,c,kd,cv,eg,varargin);
%     =========================================================
%     Purpose: Compute the characteristic values of spheroidal
%     wave functions
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     KD --- Function code
%     KD=1 for Prolate; KD=-1 for Oblate
%     Output:  CV --- Characteristic value for given m, n and c
%     EG(L)--- Characteristic value for mode m and n'
%(L = n' - m + 1)
%     =========================================================
 b=zeros(1,100);
h=zeros(1,100);
d=zeros(1,300);
e=zeros(1,300);
f=zeros(1,300);
cv0=zeros(1,100);
a=zeros(1,300);
g=zeros(1,300);
if(c < 1.0d-10);
for  i=1:n;
eg(i)=(i+fix(m)).*(i+fix(m)-1.0d0);
end;  i=fix(n)+1;
else;
icm=fix((fix(n)-fix(m)+2)./2);
nm=10+fix(0.5.*(fix(n)-fix(m))+c);
cs=c.*c.*fix(kd);
for  l=0:1;
for  i=1:nm;
if(l == 0)k=2.*(i-1); end;
if(l == 1)k=2.*i-1; end;
dk0=fix(m)+k;
dk1=fix(m)+k+1;
dk2=2.*(fix(m)+k);
d2k=2.*fix(m)+k;
a(i)=(d2k+2.0).*(d2k+1.0)./((dk2+3.0).*(dk2+5.0)).*cs;
d(i)=dk0.*dk1+(2.0.*dk0.*dk1-2.0.*fix(m).*fix(m)-1.0)./((dk2-1.0).*(dk2+3.0)).*cs;
g(i)=k.*(k-1.0)./((dk2-3.0).*(dk2-1.0)).*cs;
end;  i=nm+1;
for  k=2:nm;
e(k)=sqrt(a(k-1).*g(k));
f(k)=e(k).*e(k);
end;  k=nm+1;
f(1)=0.0d0;
e(1)=0.0d0;
xa=d(nm)+abs(e(nm));
xb=d(nm)-abs(e(nm));
nm1=nm-1;
for  i=1:nm1;
t=abs(e(i))+abs(e(i+1));
t1=d(i)+t;
if(xa < t1)xa=t1; end;
t1=d(i)-t;
if(t1 < xb)xb=t1; end;
end;  i=nm1+1;
for  i=1:icm;
b(i)=xa;
h(i)=xb;
end;  i=icm+1;
for  k=1:icm;
for  k1=k:icm;
if(b(k1)< b(k));
b(k)=b(k1);
break;
end;
end;
if(k ~= 1&h(k)< h(k-1))h(k)=h(k-1); end;
while(1==1);
x1=(b(k)+h(k))./2.0d0;
cv0(k)=x1;
if(abs((b(k)-h(k))./x1)< 1.0d-14)break; end;
j=0;
s=1.0d0;
for  i=1:nm;
if(s == 0.0d0)s=s+1.0d-30; end;
t=f(i)./s;
s=d(i)-t-x1;
if(s < 0.0d0)j=j+1; end;
end;  i=nm+1;
if(j < k);
h(k)=x1;
else;
b(k)=x1;
if(j >= icm);
b(icm)=x1;
else;
if(h(j+1)< x1)h(j+1)=x1; end;
if(x1 < b(j))b(j)=x1; end;
end;
end;
end;
cv0(k)=x1;
if(l == 0)eg(2.*k-1)=cv0(k); end;
if(l == 1)eg(2.*k)=cv0(k); end;
end;
end;
end;
cv=eg(fix(n)-fix(m)+1);
return;
end

