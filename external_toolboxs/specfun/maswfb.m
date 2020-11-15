function maswfb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ============================================================
%     Purpose: This program computes the prolate and oblate
%     spheroidal angular functions of the first kind
%     and their derivatives using subroutine ASWFB
%     Input :  m  --- Mode parameter,  m = 0,1,2,...
%     n  --- Mode parameter,  n = m,m+1,...
%     c  --- Spheroidal parameter
%     x  --- Argument of angular function, |x| ó 1.0
%     KD --- Function code
%     KD=1 for prolate;  KD=-1 for oblate
%     cv --- Characteristic value
%     Output:  S1F --- Angular function of the first kind
%     S1D --- Derivative of the angular function of
%     the first kind
%     Examples:
%     KD = 1, m = 2, n = 3, c = 3.0 and cv = 14.8277782138
%     x         Smn(c,x)Smn'(c,x)
%     --------------------------------------------
%     0.2      .28261309D+01       .12418631D+02
%     0.5      .49938554D+01       .92761604D+00
%     0.8      .31693975D+01      -.12646552D+02
%     KD =-1, m = 2, n = 3, c = 3.0 and cv = 8.8093939208
%     x         Smn(-ic,x)Smn'(-ic,x)
%     --------------------------------------------
%     0.2      .29417848D+01       .14106305D+02
%     0.5      .64138827D+01       .76007194D+01
%     0.8      .60069873D+01      -.14387479D+02
%     ============================================================
m=[];n=[];c=[];kd=[];cv=[];eg=[];x=[];s1f=[];s1d=[];
 eg=zeros(1,200);
fprintf(1,'%s \n','please enter kd, m, n and c ');
%     READ(*,*)KD,M,N,C
m=2;
n=3;
c=3.0;
kd=1;
fprintf(1,[repmat(' ',1,1),'kd =','%2g',', ','m =','%2g',', ','n =','%2g',', ','c =','%5.1g' ' \n'],kd,m,n,c);
[m,n,c,kd,cv,eg]=segv(m,n,c,kd,cv,eg);
fprintf(1,[repmat(' ',1,1),' cv =','%18.10g' ' \n'],cv);
fprintf(1,'%0.15g \n');
if(kd == 1);
fprintf(1,'%s \n','    x         smn(c,x)smn''(c,x)');
elseif(kd == -1);
fprintf(1,'%s \n','    x         smn(-ic,x)smn''(-ic,x)');
end;
fprintf(1,'%s \n','  --------------------------------------------');
for  i=0:20;
x=-1.0d0+0.1d0.*i;
if(i == 0)x=0.0d0; end;
if(i == 20)x=1.0d0; end;
[m,n,c,x,kd,cv,s1f,s1d]=aswfb(m,n,c,x,kd,cv,s1f,s1d);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%20.8g',1,2) ' \n'],x,s1f,s1d);
end;  i=20+1;
%format(1x,'i2,', ',',i2,', ',',i2,', ',',f5.1);
%format(1x,',f18.10);
%format(1x,f5.1,2d20.8);
end
function [m,n,c,x,kd,cv,s1f,s1d]=aswfb(m,n,c,x,kd,cv,s1f,s1d,varargin);
%     ===========================================================
%     Purpose: Compute the prolate and oblate spheroidal angular
%     functions of the first kind and their derivatives
%     Input :  m  --- Mode parameter,  m = 0,1,2,...
%     n  --- Mode parameter,  n = m,m+1,...
%     c  --- Spheroidal parameter
%     x  --- Argument of angular function, |x| < 1.0
%     KD --- Function code
%     KD=1 for prolate;  KD=-1 for oblate
%     cv --- Characteristic value
%     Output:  S1F --- Angular function of the first kind
%     S1D --- Derivative of the angular function of
%     the first kind
%     Routines called:
%(1)SDMN for computing expansion coefficients dk
%(2)LPMNS for computing associated Legendre function
%     of the first kind Pmn(x)
%     ===========================================================
df=[];nm2=[];pm=[];pd=[];
 df=zeros(1,200);
pm=zeros(1,251+1);
pd=zeros(1,251+1);
eps=1.0d-14;
ip=1;
sw=0.0;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
nm=25+fix((fix(n)-fix(m))./2+c);
nm2=2.*nm+fix(m);
[m,n,c,cv,kd,df]=sdmn(fix(m),fix(n),c,cv,fix(kd),df);
[m,nm2,x,pm,pd]=lpmns(fix(m),nm2,x,pm,pd);
su1=0.0d0;
for  k=1:nm;
mk=fix(m)+2.*(k-1)+ip;
su1=su1+df(k).*pm(mk+1+1);
if(abs(sw-su1)< abs(su1).*eps)break; end;
sw=su1;
end;
s1f=(-1).^fix(m).*su1;
su1=0.0d0;
for  k=1:nm;
mk=fix(m)+2.*(k-1)+ip;
su1=su1+df(k).*pd(mk+1+1);
if(abs(sw-su1)< abs(su1).*eps)break; end;
sw=su1;
end;
s1d=(-1).^fix(m).*su1;
return;
end
function [m,n,c,cv,kd,df]=sdmn(m,n,c,cv,kd,df,varargin);
%     =====================================================
%     Purpose: Compute the expansion coefficients of the
%     prolate and oblate spheroidal functions, dk
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     cv --- Characteristic value
%     KD --- Function code
%     KD=1 for prolate; KD=-1 for oblate
%     Output:  DF(k)--- Expansion coefficients dk;
%     DF(1), DF(2), ... correspond to
%     d0, d2, ... for even n-m and d1,
%     d3, ... for odd n-m
%     =====================================================
 a=zeros(1,200);
d=zeros(1,200);
g=zeros(1,200);
sw=0.0;
fl=0.0;
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
function [m,n,x,pm,pd]=lpmns(m,n,x,pm,pd,varargin);
%     ========================================================
%     Purpose: Compute associated Legendre functions Pmn(x)
%     and Pmn'(x)for a given order
%     Input :  x --- Argument of Pmn(x)
%     m --- Order of Pmn(x),  m = 0,1,2,...,n
%     n --- Degree of Pmn(x), n = 0,1,2,...,N
%     Output:  PM(n)--- Pmn(x)
%     PD(n)--- Pmn'(x)
%     ========================================================
for  k=0:n;
pm(k+1+1)=0.0d0;
pd(k+1+1)=0.0d0;
end;  k=fix(n)+1;
if(abs(x)== 1.0d0);
for  k=0:n;
if(m == 0);
pm(k+1+1)=1.0d0;
pd(k+1+1)=0.5d0.*k.*(k+1.0);
if(x < 0.0);
pm(k+1+1)=(-1).^k.*pm(k+1+1);
pd(k+1+1)=(-1).^(k+1).*pd(k+1+1);
end;
elseif(m == 1);
pd(k+1+1)=1.0d+300;
elseif(m == 2);
pd(k+1+1)=-0.25d0.*(k+2.0).*(k+1.0).*k.*(k-1.0);
if(x < 0.0)pd(k+1+1)=(-1).^(k+1).*pd(k+1+1); end;
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
pm(m+1+1)=pmk;
pm(m+1+1+1)=pm1;
for  k=m+2:n;
pm2=((2.0d0.*k-1.0d0).*x.*pm1-(k+fix(m)-1.0d0).*pmk)./(k-fix(m));
pm(k+1+1)=pm2;
pmk=pm1;
pm1=pm2;
end;  k=fix(n)+1;
for  k=1:n;
pm(k+1)=(-1).^fix(m).*pm(k+1);
pd(k+1)=(-1).^fix(m).*pd(k+1);
end;  k=fix(n)+1;
pd(0+1+1)=((1.0d0-fix(m)).*pm(1+1+1)-x.*pm(0+1+1))./(x.*x-1.0);
for  k=1:n;
pd(k+1+1)=(k.*x.*pm(k+1+1)-(k+fix(m)).*pm(k-1+1+1))./(x.*x-1.0d0);
end;  k=fix(n)+1;
return;
end

