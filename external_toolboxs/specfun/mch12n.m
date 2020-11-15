function mch12n
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     =====================================================
%     Purpose: This program computes Hankel functions of
%     the first and second kinds and their
%     derivatives for a complex argument using
%     subroutine CH12N
%     Input :  z --- Complex argument
%     n --- Order of Hn(1,z)and Hn(2,z)
%(n = 0,1,תתת, n ף 250)
%     Output:  CHF1(n)--- Hn(1,z)
%     CHD1(n)--- Hn(1)'(z)
%     CHF2(n)--- Hn(2,z)
%     CHD2(n)--- Hn(2)'(z)
%     =====================================================
n=[];z=[];nm=[];chf1=[];chd1=[];chf2=[];chd2=[];
 chf1=zeros(1,250+1);
chd1=zeros(1,250+1);
chf2=zeros(1,250+1);
chd2=zeros(1,250+1);
fprintf(1,'%s \n','  please enter n, x and y(z=x+iy)');
%      READ(*,*)N,X,Y
n=5;
x=1.5;
y=1.5;
nm=1;
fprintf(1,[repmat(' ',1,3),'z =','%8.3g',' + i ','%8.3g',' ,',repmat(' ',1,6),'nmax =','%4g' ' \n'],x,y,n);
z=complex(x,y);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%       READ(*,*)NS
ns=2;
end;
[n,z,nm,chf1,chd1,chf2,chd2]=ch12n(n,z,nm,chf1,chd1,chf2,chd2);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n     re[hn(1,z)]im[hn(1,z)]');fprintf(1,'%s \n','      re[hn(1)''(z)]im[hn(1)''(z)]');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','---------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat('%18.10g',1,4) ' \n'],k,chf1(k+1),chd1(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n     re[hn(2,z)]im[hn(2,z)]');fprintf(1,'%s \n','      re[hn(2)''(z)]im[hn(2)''(z)]');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','---------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat('%18.10g',1,4) ' \n'],k,chf2(k+1),chd2(k+1));
end;  k=nm+1;
%format(1x,i4,4d18.10);
%format(3x,,f8.3,' + i ',f8.3,' ,',6x,,i4);
end
function [n,z,nm,chf1,chd1,chf2,chd2]=ch12n(n,z,nm,chf1,chd1,chf2,chd2,varargin);
%     ====================================================
%     Purpose: Compute Hankel functions of the first and
%     second kinds and their derivatives for a
%     complex argument
%     Input :  z --- Complex argument
%     n --- Order of Hn(1,z)and Hn(2,z)
%     Output:  CHF1(n)--- Hn(1,z)
%     CHD1(n)--- Hn(1)'(z)
%     CHF2(n)--- Hn(2,z)
%     CHD2(n)--- Hn(2)'(z)
%     NM --- Highest order computed
%     Routines called:
%(1)CJYNB for computing Jn(z)and Yn(z)
%(2)CIKNB for computing In(z)and Kn(z)
%     ====================================================
cbj=[];cdj=[];cby=[];cdy=[];zi=[];cbi=[];cdi=[];cbk=[];cdk=[];
 cbj=zeros(1,250+1);
cdj=zeros(1,250+1);
cby=zeros(1,250+1);
cdy=zeros(1,250+1);
cbi=zeros(1,250+1);
cdi=zeros(1,250+1);
cbk=zeros(1,250+1);
cdk=zeros(1,250+1);
ci=complex(0.0d0,1.0d0);
pi=3.141592653589793d0;
if(imag(z)< 0.0d0);
[n,z,nm,cbj,cdj,cby,cdy]=cjynb(fix(n),z,fix(nm),cbj,cdj,cby,cdy);
for  k=0:nm;
chf1(k+1)=cbj(k+1)+ci.*cby(k+1);
chd1(k+1)=cdj(k+1)+ci.*cdy(k+1);
end;  k=fix(nm)+1;
zi=ci.*z;
[n,zi,nm,cbi,cdi,cbk,cdk]=ciknb(fix(n),zi,fix(nm),cbi,cdi,cbk,cdk);
cfac=-2.0d0./(pi.*ci);
for  k=0:nm;
chf2(k+1)=cfac.*cbk(k+1);
chd2(k+1)=cfac.*ci.*cdk(k+1);
cfac=cfac.*ci;
end;  k=fix(nm)+1;
elseif(imag(z)> 0.0d0);
zi=-ci.*z;
[n,zi,nm,cbi,cdi,cbk,cdk]=ciknb(fix(n),zi,fix(nm),cbi,cdi,cbk,cdk);
cf1=-ci;
cfac=2.0d0./(pi.*ci);
for  k=0:nm;
chf1(k+1)=cfac.*cbk(k+1);
chd1(k+1)=-cfac.*ci.*cdk(k+1);
cfac=cfac.*cf1;
end;  k=fix(nm)+1;
[n,z,nm,cbj,cdj,cby,cdy]=cjynb(fix(n),z,fix(nm),cbj,cdj,cby,cdy);
for  k=0:nm;
chf2(k+1)=cbj(k+1)-ci.*cby(k+1);
chd2(k+1)=cdj(k+1)-ci.*cdy(k+1);
end;  k=fix(nm)+1;
else;
[n,z,nm,cbj,cdj,cby,cdy]=cjynb(fix(n),z,fix(nm),cbj,cdj,cby,cdy);
for  k=0:nm;
chf1(k+1)=cbj(k+1)+ci.*cby(k+1);
chd1(k+1)=cdj(k+1)+ci.*cdy(k+1);
chf2(k+1)=cbj(k+1)-ci.*cby(k+1);
chd2(k+1)=cdj(k+1)-ci.*cdy(k+1);
end;  k=fix(nm)+1;
end;
return;
end
function [n,z,nm,cbj,cdj,cby,cdy]=cjynb(n,z,nm,cbj,cdj,cby,cdy,varargin);
%     =======================================================
%     Purpose: Compute Bessel functions Jn(z), Yn(z)and
%     their derivatives for a complex argument
%     Input :  z --- Complex argument of Jn(z)and Yn(z)
%     n --- Order of Jn(z)and Yn(z)
%     Output:  CBJ(n)--- Jn(z)
%     CDJ(n)--- Jn'(z)
%     CBY(n)--- Yn(z)
%     CDY(n)--- Yn'(z)
%     NM --- Highest order computed
%     Routines called:
%     MSTA1 and MSTA2 to calculate the starting
%     point for backward recurrence
%     =======================================================
  a=zeros(1,4);
b=zeros(1,4);
a1=zeros(1,4);
b1=zeros(1,4);
el=0.5772156649015329d0;
pi=3.141592653589793d0;
r2p=.63661977236758d0;
y0=abs(imag(z));
a0=abs(z);
nm=fix(fix(n));
if(a0 < 1.0d-100);
for  k=0:n;
cbj(k+1)=complex(0.0d0,0.0d0);
cdj(k+1)=complex(0.0d0,0.0d0);
cby(k+1)=-complex(1.0d+300,0.0d0);
cdy(k+1)=complex(1.0d+300,0.0d0);
end;  k=fix(n)+1;
cbj(0+1)=complex(1.0d0,0.0d0);
cdj(1+1)=complex(0.5d0,0.0d0);
return;
end;
if(a0 <= 300.d0|n > 80);
if(n == 0)nm=1; end;
m=fix(msta1(a0,200));
if(m < nm);
nm=fix(m);
else;
m=fix(msta2(a0,fix(nm),15));
end;
cbs=complex(0.0d0,0.0d0);
csu=complex(0.0d0,0.0d0);
csv=complex(0.0d0,0.0d0);
cf2=complex(0.0d0,0.0d0);
cf1=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=2.0d0.*(k+1.0d0)./z.*cf1-cf2;
if(k <= nm)cbj(k+1)=cf; end;
if(k == 2.*fix(k./2)&k ~= 0);
if(y0 <= 1.0d0);
cbs=cbs+2.0d0.*cf;
else;
cbs=cbs+(-1).^(k./2).*2.0d0.*cf;
end;
csu=csu+(-1).^(k./2).*cf./k;
elseif(k > 1);
csv=csv+(-1).^(k./2).*k./(k.*k-1.0d0).*cf;
end;
cf2=cf1;
cf1=cf;
end;  k=0-1;
if(y0 <= 1.0d0);
cs0=cbs+cf;
else;
cs0=(cbs+cf)./cos(z);
end;
for  k=0:nm;
cbj(k+1)=cbj(k+1)./cs0;
end;  k=fix(nm)+1;
ce=log(z./2.0d0)+el;
cby(0+1)=r2p.*(ce.*cbj(0+1)-4.0d0.*csu./cs0);
cby(1+1)=r2p.*(-cbj(0+1)./z+(ce-1.0d0).*cbj(1+1)-4.0d0.*csv./cs0);
else;
a(:)=[-.7031250000000000d-01,.1121520996093750d+00,-.5725014209747314d+00,.6074042001273483d+01];
b(:)=[.7324218750000000d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02];
a1(:)=[.1171875000000000d+00,-.1441955566406250d+00,.6765925884246826d+00,-.6883914268109947d+01];
b1(:)=[-.1025390625000000d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02];
ct1=z-0.25d0.*pi;
cp0=complex(1.0d0,0.0d0);
for  k=1:4;
cp0=cp0+a(k).*z.^(-2.*k);
end;  k=4+1;
cq0=-0.125d0./z;
for  k=1:4;
cq0=cq0+b(k).*z.^(-2.*k-1);
end;  k=4+1;
cu=sqrt(r2p./z);
cbj0=cu.*(cp0.*cos(ct1)-cq0.*sin(ct1));
cby0=cu.*(cp0.*sin(ct1)+cq0.*cos(ct1));
cbj(0+1)=cbj0;
cby(0+1)=cby0;
ct2=z-0.75d0.*pi;
cp1=complex(1.0d0,0.0d0);
for  k=1:4;
cp1=cp1+a1(k).*z.^(-2.*k);
end;  k=4+1;
cq1=0.375d0./z;
for  k=1:4;
cq1=cq1+b1(k).*z.^(-2.*k-1);
end;  k=4+1;
cbj1=cu.*(cp1.*cos(ct2)-cq1.*sin(ct2));
cby1=cu.*(cp1.*sin(ct2)+cq1.*cos(ct2));
cbj(1+1)=cbj1;
cby(1+1)=cby1;
for  k=2:nm;
cbjk=2.0d0.*(k-1.0d0)./z.*cbj1-cbj0;
cbj(k+1)=cbjk;
cbj0=cbj1;
cbj1=cbjk;
end;  k=fix(nm)+1;
end;
cdj(0+1)=-cbj(1+1);
for  k=1:nm;
cdj(k+1)=cbj(k-1+1)-k./z.*cbj(k+1);
end;  k=fix(nm)+1;
if(abs(cbj(0+1))> 1.0d0);
cby(1+1)=(cbj(1+1).*cby(0+1)-2.0d0./(pi.*z))./cbj(0+1);
end;
for  k=2:nm;
if(abs(cbj(k-1+1))>= abs(cbj(k-2+1)));
cyy=(cbj(k+1).*cby(k-1+1)-2.0d0./(pi.*z))./cbj(k-1+1);
else;
cyy=(cbj(k+1).*cby(k-2+1)-4.0d0.*(k-1.0d0)./(pi.*z.*z))./cbj(k-2+1);
end;
cby(k+1)=cyy;
end;  k=fix(nm)+1;
cdy(0+1)=-cby(1+1);
for  k=1:nm;
cdy(k+1)=cby(k-1+1)-k./z.*cby(k+1);
end;  k=fix(nm)+1;
return;
end
function [n,z,nm,cbi,cdi,cbk,cdk]=ciknb(n,z,nm,cbi,cdi,cbk,cdk,varargin);
%     ============================================================
%     Purpose: Compute modified Bessel functions In(z)and Kn(z),
%     and their derivatives for a complex argument
%     Input:   z --- Complex argument
%     n --- Order of In(z)and Kn(z)
%     Output:  CBI(n)--- In(z)
%     CDI(n)--- In'(z)
%     CBK(n)--- Kn(z)
%     CDK(n)--- Kn'(z)
%     NM --- Highest order computed
%     Routones called:
%     MSTA1 and MSTA2 to compute the starting point for
%     backward recurrence
%     ===========================================================
pi=3.141592653589793d0;
el=0.57721566490153d0;
a0=abs(z);
nm=fix(fix(n));
if(a0 < 1.0d-100);
for  k=0:n;
cbi(k+1)=complex(0.0d0,0.0d0);
cbk(k+1)=complex(1.0d+300,0.0d0);
cdi(k+1)=complex(0.0d0,0.0d0);
cdk(k+1)=-complex(1.0d+300,0.0d0);
end;  k=fix(n)+1;
cbi(0+1)=complex(1.0d0,0.0d0);
cdi(1+1)=complex(0.5d0,0.0d0);
return;
end;
z1=z;
ci=complex(0.0d0,1.0d0);
if(real(z)< 0.0)z1=-z; end;
if(n == 0)nm=1; end;
m=fix(msta1(a0,200));
if(m < nm);
nm=fix(m);
else;
m=fix(msta2(a0,fix(nm),15));
end;
cbs=0.0d0;
csk0=0.0d0;
cf0=0.0d0;
cf1=1.0d-100;
for  k=m:-1:0;
cf=2.0d0.*(k+1.0d0).*cf1./z1+cf0;
if(k <= nm)cbi(k+1)=cf; end;
if(k ~= 0&k == 2.*fix(k./2))csk0=csk0+4.0d0.*cf./k; end;
cbs=cbs+2.0d0.*cf;
cf0=cf1;
cf1=cf;
end;  k=0-1;
cs0=exp(z1)./(cbs-cf);
for  k=0:nm;
cbi(k+1)=cs0.*cbi(k+1);
end;  k=fix(nm)+1;
if(a0 <= 9.0);
cbk(0+1)=-(log(0.5d0.*z1)+el).*cbi(0+1)+cs0.*csk0;
cbk(1+1)=(1.0d0./z1-cbi(1+1).*cbk(0+1))./cbi(0+1);
else;
ca0=sqrt(pi./(2.0d0.*z1)).*exp(-z1);
k0=16;
if(a0 >= 25.0)k0=10; end;
if(a0 >= 80.0)k0=8; end;
if(a0 >= 200.0)k0=6; end;
for  l=0:1;
cbkl=1.0d0;
vt=4.0d0.*l;
cr=complex(1.0d0,0.0d0);
for  k=1:k0;
cr=0.125d0.*cr.*(vt-(2.0.*k-1.0).^2)./(k.*z1);
cbkl=cbkl+cr;
end;  k=k0+1;
cbk(l+1)=ca0.*cbkl;
end;  l=1+1;
end;
cg0=cbk(0+1);
cg1=cbk(1+1);
for  k=2:nm;
cg=2.0d0.*(k-1.0d0)./z1.*cg1+cg0;
cbk(k+1)=cg;
cg0=cg1;
cg1=cg;
end;  k=fix(nm)+1;
if(real(z)< 0.0);
fac=1.0d0;
for  k=0:nm;
if(imag(z)< 0.0);
cbk(k+1)=fac.*cbk(k+1)+ci.*pi.*cbi(k+1);
else;
cbk(k+1)=fac.*cbk(k+1)-ci.*pi.*cbi(k+1);
end;
cbi(k+1)=fac.*cbi(k+1);
fac=-fac;
end;  k=fix(nm)+1;
end;
cdi(0+1)=cbi(1+1);
cdk(0+1)=-cbk(1+1);
for  k=1:nm;
cdi(k+1)=cbi(k-1+1)-k./z.*cbi(k+1);
cdk(k+1)=-cbk(k-1+1)-k./z.*cbk(k+1);
end;  k=fix(nm)+1;
return;
end
function [msta1Result]=msta1(x,mp,varargin);
%     ===================================================
%     Purpose: Determine the starting point for backward
%     recurrence such that the magnitude of
%     Jn(x)at that point is about 10^(-MP)
%     Input :  x     --- Argument of Jn(x)
%     MP    --- Value of magnitude
%     Output:  MSTA1 --- Starting point
%     ===================================================
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
msta1Result=nn;
return;
end
function [msta2Result]=msta2(x,n,mp,varargin);
%     ===================================================
%     Purpose: Determine the starting point for backward
%     recurrence such that all Jn(x)has MP
%     significant digits
%     Input :  x  --- Argument of Jn(x)
%     n  --- Order of Jn(x)
%     MP --- Significant digit
%     Output:  MSTA2 --- Starting point
%     ===================================================
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
msta2Result=nn+10;
return;
end
function [envjResult]=envj(n,x,varargin);
envjResult=0.5d0.*log10(6.28d0.*fix(n))-fix(n).*log10(1.36d0.*x./fix(n));
return;
end

