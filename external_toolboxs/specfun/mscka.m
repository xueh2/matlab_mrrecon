function mscka
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ============================================================
%     Purpose: This program computes the expansion coefficients
%     of the prolate and oblate spheroidal functions,
%     c2k, using subroutine SCKA
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     cv --- Characteristic value
%     KD --- Function code
%     KD=1 for prolate; KD=-1 for oblate
%     Output:  CK(k)--- Expansion coefficients ck;
%     CK(1), CK(2),... correspond to
%     c0, c2,...
%     Example: Compute the first 13 expansion coefficients C2k for
%     KD= 1, m=2, n=3, c=3.0 and cv=14.8277782138; and
%     KD=-1, m=2, n=3, c=3.0 and cv=8.80939392077
%     Coefficients of Prolate and oblate functions
%     k         C2k(c)C2k(-ic)
%     ---------------------------------------------
%     0     .9173213327D+01     .2489664942D+02
%     1     .4718258929D+01    -.1205287032D+02
%     2     .9841212916D+00     .2410564082D+01
%     3     .1151870224D+00    -.2735821590D+00
%     4     .8733916403D-02     .2026057157D-01
%     5     .4663888254D-03    -.1061946315D-02
%     6     .1853910398D-04     .4158091152D-04
%     7     .5708084895D-06    -.1264400411D-05
%     8     .1402786472D-07     .3074963448D-07
%     9     .2817194508D-09    -.6120579463D-09
%     10     .4712094447D-11     .1015900041D-10
%     11     .6667838485D-13    -.1427953361D-12
%     12     .8087995432D-15     .1721924955D-14
%     ============================================================
m=[];n=[];c=[];kd=[];cv=[];eg=[];ck=[];
 ck=zeros(1,200);
eg=zeros(1,200);
fprintf(1,'%s \n','please kd, m, n and c ');
%     READ(*,*)KD,M,N,C
kd=1;
m=2;
n=3;
c=3.0;
[m,n,c,kd,cv,eg]=segv(m,n,c,kd,cv,eg);
fprintf(1,[repmat(' ',1,1),'kd=','%3g',',  ','m=','%3g',',  ','n=','%3g',',  ','c=','%5.1g',',  ','cv =','%18.10g' ' \n'],kd,m,n,c,cv);
[m,n,c,cv,kd,ck]=scka(m,n,c,cv,kd,ck);
fprintf(1,'%0.15g \n');
if(kd == 1);
fprintf(1,'%s \n','coefficients of prolate function');
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   k            c2k(c)');
else;
fprintf(1,'%s \n','coefficients of oblate function');
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   k           c2k(-ic)');
end;
fprintf(1,'%s \n','----------------------------');
nm=25+fix((n-m)./2+c);
for  k=1:nm;
fprintf(1,[repmat(' ',1,2),'%3g',repmat(' ',1,4),'%18.10g' ' \n'],k-1,ck(k));
end;  k=nm+1;
%format(2x,i3,4x,d18.10);
%format(1x,,i3,',  ',,i3,',  ',,i3,',  ',,f5.1,',  ',,f18.10);
end
function [m,n,c,cv,kd,ck]=scka(m,n,c,cv,kd,ck,varargin);
%     ======================================================
%     Purpose: Compute the expansion coefficients of the
%     prolate and oblate spheroidal functions, c2k
%     Input :  m  --- Mode parameter
%     n  --- Mode parameter
%     c  --- Spheroidal parameter
%     cv --- Characteristic value
%     KD --- Function code
%     KD=1 for prolate; KD=-1 for oblate
%     Output:  CK(k)--- Expansion coefficients ck;
%     CK(1), CK(2),... correspond to
%     c0, c2,...
%     ======================================================
if(c <= 1.0d-10)c=1.0d-10; end;
nm=25+fix((fix(n)-fix(m))./2+c);
cs=c.*c.*fix(kd);
ip=1;
if(n-m == 2.*fix((n-m)./2))ip=0; end;
fs=1.0d0;
f1=0.0d0;
f0=1.0d-100;
kb=0;
ck(nm+1)=0.0d0;
for  k=nm:-1:1;
f=(((2.0d0.*k+fix(m)+ip).*(2.0d0.*k+fix(m)+1.0d0+ip)-cv+cs).*f0-4.0d0.*(k+1.0d0).*(k+fix(m)+1.0d0).*f1)./cs;
if(abs(f)> abs(ck(k+1)));
ck(k)=f;
f1=f0;
f0=f;
if(abs(f)> 1.0d+100);
for  k1=nm:-1:k;
ck(k1)=ck(k1).*1.0d-100;
end;  k1=k-1;
f1=f1.*1.0d-100;
f0=f0.*1.0d-100;
end;
else;
kb=k;
fl=ck(k+1);
f1=1.0d0;
f2=0.25d0.*((fix(m)+ip).*(fix(m)+ip+1.0)-cv+cs)./(fix(m)+1.0).*f1;
ck(1)=f1;
if(kb == 1);
fs=f2;
elseif(kb == 2);
ck(2)=f2;
fs=0.125d0.*(((fix(m)+ip+2.0).*(fix(m)+ip+3.0)-cv+cs).*f2 -cs.*f1)./(fix(m)+2.0);
else;
ck(2)=f2;
for  j=3:kb+1;
f=0.25d0.*(((2.0.*j+fix(m)+ip-4.0).*(2.0.*j+fix(m)+ip-3.0)-cv+cs).*f2-cs.*f1)./((j-1.0).*(j+fix(m)-1.0));
if(j <= kb)ck(j)=f; end;
f1=f2;
f2=f;
end;  j=kb+1+1;
fs=f;
end;
break;
end;
end;
su1=0.0d0;
for  k=1:kb;
su1=su1+ck(k);
end;  k=kb+1;
su2=0.0d0;
for  k=kb+1:nm;
su2=su2+ck(k);
end;  k=nm+1;
r1=1.0d0;
for  j=1:(n+m+ip)./2;
r1=r1.*(j+0.5d0.*(fix(n)+fix(m)+ip));
end;  j=(fix(n)+fix(m)+ip)./2+1;
r2=1.0d0;
for  j=1:(n-m-ip)./2;
r2=-r2.*j;
end;  j=(fix(n)-fix(m)-ip)./2+1;
if(kb == 0);
s0=r1./(2.0d0.^fix(n).*r2.*su2);
else;
s0=r1./(2.0d0.^fix(n).*r2.*(fl./fs.*su1+su2));
end;
for  k=1:kb;
ck(k)=fl./fs.*s0.*ck(k);
end;  k=kb+1;
for  k=kb+1:nm;
ck(k)=s0.*ck(k);
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

