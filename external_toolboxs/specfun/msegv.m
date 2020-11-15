function msegv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes a sequence of characteristic
%                values of spheroidal prolate and oblate functions
%                using subroutine SEGV
%       Input :  m  --- Mode parameter
%                n  --- Mode parameter
%                c  --- Spheroidal parameter
%                KD --- Function code
%                       KD=1 for Prolate; KD=-1 for Oblate
%       Output:  CV --- Characteristic value for given m, n and c
%                EG(L)--- Characteristic value for mode m and n'
%(L = n' - m + 1)
%       Examples:
%                Prolate:(KD = 1 , m = 1, n = 5, c = 5.0)
%                m      n       c        Lambda mn(c)
%              ---------------------------------------
%                1      1      5.0        5.35042230
%                1      2      5.0       14.64295624
%                1      3      5.0       23.39761312
%                1      4      5.0       32.42194359
%                1      5      5.0       42.65818215
%                Oblate:(KD = -1 , m = 1, n = 5, c = 5.0)
%                m      n       c      Lambda mn(-ic)
%               --------------------------------------
%                1      1      5.0       -7.49338828
%                1      2      5.0       -7.12783752
%                1      3      5.0        2.75036721
%                1      4      5.0        8.69495925
%                1      5      5.0       18.43931577
%       =========================================================
m=[];n=[];c=[];kd=[];cv=[];eg=[];
 eg=zeros(1,100);
fprintf(1,'%s \n','please enter kd, m, n and c ');
%        READ(*,*)KD,M,N,C
kd=1;
m=1;
n=5;
c=5.0;
fprintf(1,[repmat(' ',1,1),'kd =','%2g',',  m =','%3g',',   n =','%3g',',  c =','%5.1g' ' \n'],kd,m,n,c);
fprintf(1,'%0.15g \n');
[m,n,c,kd,cv,eg]=segv(m,n,c,kd,cv,eg);
if(kd == 1);
fprintf(1,'%s \n','  m      n       c       lambda mn(c)');
elseif(kd == -1);
fprintf(1,'%s \n','  m      n       c      lambda mn(-ic)');
end;
fprintf(1,'%s \n','---------------------------------------');
for  l=1:n-m+1;
n1=m+l-1;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,4),'%3g',repmat(' ',1,4),'%5.1g','%18.8g' ' \n'],m,n1,c,eg(l));
end;  l=n-m+1+1;
%format(1x,',i2,',i3,',i3,',f5.1);
%format(1x,i3,4x,i3,4x,f5.1,f18.8);
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

