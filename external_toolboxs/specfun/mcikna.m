function mcikna
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes the modified Bessel functions
%                In(z)and Kn(z), and their derivatives for a
%                complex argument using subroutine CIKNA
%       Input :  z --- Complex argument of In(z)and Kn(z)
%                n --- Order of In(z)and Kn(z)
%(n = 0,1,תתת, n ף 250)
%       Output:  CBI(n)--- In(z)
%                CDI(n)--- In'(z)
%                CBK(n)--- Kn(z)
%                CDK(n)--- Kn'(z)
%       Example: z = 4.0 + i 2.0 ,      Nmax = 5
%     n     Re[In(z)]Im[In(z)]Re[In'(z)]Im[In'(z)]
%   -----------------------------------------------------------------
%     0  -.19056142D+01  .10403505D+02 -.23059657D+01  .92222463D+01
%     1  -.23059657D+01  .92222463D+01 -.23666457D+01  .83284588D+01
%     2  -.28276772D+01  .62534130D+01 -.24255774D+01  .61553456D+01
%     3  -.25451891D+01  .30884450D+01 -.22270972D+01  .36367893D+01
%     4  -.16265172D+01  .10201656D+01 -.16520416D+01  .16217056D+01
%     5  -.75889410D+00  .15496632D+00 -.94510625D+00  .48575220D+00
%     n     Re[Kn(z)]Im[Kn(z)]Re[Kn'(z)]Im[Kn'(z)]
%   -----------------------------------------------------------------
%     0  -.64221754D-02 -.84393648D-02  .74307276D-02  .89585853D-02
%     1  -.74307276D-02 -.89585853D-02  .88041795D-02  .94880091D-02
%     2  -.11186184D-01 -.10536653D-01  .14012532D-01  .10936010D-01
%     3  -.20594336D-01 -.12913435D-01  .27416815D-01  .12106413D-01
%     4  -.43647447D-01 -.13676173D-01  .60982763D-01  .63953943D-02
%     5  -.10137119D+00  .12264588D-03  .14495731D+00 -.37132068D-01
%       =============================================================
n=[];z=[];nm=[];
 global cbi;
global cdi;
global cbk;
global cdk;
fprintf(1,'%s \n','  please input n, x,y(z=x+iy)=?');
%        READ(*,*)N,X,Y
n=5;
x=4.0;
y=2.0;
z=complex(x,y);
fprintf(1,[repmat(' ',1,3),'z =','%7.1g',' + i','%7.1g',' ,',repmat(' ',1,6),'nmax =','%4g' ' \n'],x,y,n);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n',' please enter order step ns');
%         READ(*,*)NS
ns=1;
end;
[n,z,nm,cbi,cdi,cbk,cdk]=cikna(n,z,nm,cbi,cdi,cbk,cdk);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n      re[in(z)]im[in(z)]');fprintf(1,'%s \n','       re[in''(z)]im[in''(z)]');
fprintf(1,'%s ',' -----------------------------------');fprintf(1,'%s \n','----------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,cbi(k+1),cdi(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n      re[kn(z)]im[kn(z)]');fprintf(1,'%s \n','       re[kn''(z)]im[kn''(z)]');
fprintf(1,'%s ',' -----------------------------------');fprintf(1,'%s \n','----------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,cbk(k+1),cdk(k+1));
end;  k=nm+1;
%format(1x,i4,1x,4d16.8);
%format(3x,,f7.1,' + i',f7.1,' ,',6x,,i4);
end
function [n,z,nm,cbi,cdi,cbk,cdk]=cikna(n,z,nm,cbi,cdi,cbk,cdk,varargin);
%       ========================================================
%       Purpose: Compute modified Bessel functions In(z), Kn(x)
%                and their derivatives for a complex argument
%       Input :  z --- Complex argument of In(z)and Kn(z)
%                n --- Order of In(z)and Kn(z)
%       Output:  CBI(n)--- In(z)
%                CDI(n)--- In'(z)
%                CBK(n)--- Kn(z)
%                CDK(n)--- Kn'(z)
%                NM --- Highest order computed
%       Routines called:
%(1)CIK01 to compute I0(z), I1(z)K0(z)& K1(z)
%(2)MSTA1 and MSTA2 to compute the starting
%                 point for backward recurrence
%       ========================================================
cbi0=[];cdi0=[];cbi1=[];cdi1=[];cbk0=[];cdk0=[];cbk1=[];cdk1=[];
a0=abs(z);
nm=fix(fix(n));
if(a0 < 1.0d-100);
for  k=0:n;
cbi(k+1)=complex(0.0d0,0.0d0);
cdi(k+1)=complex(0.0d0,0.0d0);
cbk(k+1)=-complex(1.0d+300,0.0d0);
cdk(k+1)=complex(1.0d+300,0.0d0);
end;  k=fix(n)+1;
cbi(0+1)=complex(1.0d0,0.0d0);
cdi(1+1)=complex(0.5d0,0.0d0);
return;
end;
[z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1]=cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1);
cbi(0+1)=cbi0;
cbi(1+1)=cbi1;
cbk(0+1)=cbk0;
cbk(1+1)=cbk1;
cdi(0+1)=cdi0;
cdi(1+1)=cdi1;
cdk(0+1)=cdk0;
cdk(1+1)=cdk1;
if(n <= 1)return; end;
m=msta1(a0,200);
if(m < n);
nm=fix(m);
else;
m=msta2(a0,fix(n),15);
end;
cf2=complex(0.0d0,0.0d0);
cf1=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=2.0d0.*(k+1.0d0)./z.*cf1+cf2;
if(k <= nm)cbi(k+1)=cf; end;
cf2=cf1;
cf1=cf;
end;  k=0-1;
cs=cbi0./cf;
for  k=0:nm;
cbi(k+1)=cs.*cbi(k+1);
end;  k=fix(nm)+1;
for  k=2:nm;
if(abs(cbi(k-1+1))> abs(cbi(k-2+1)));
ckk=(1.0d0./z-cbi(k+1).*cbk(k-1+1))./cbi(k-1+1);
else;
ckk=(cbi(k+1).*cbk(k-2+1)+2.0d0.*(k-1.0d0)./(z.*z))./cbi(k-2+1);
end;
cbk(k+1)=ckk;
end;  k=fix(nm)+1;
for  k=2:nm;
cdi(k+1)=cbi(k-1+1)-k./z.*cbi(k+1);
cdk(k+1)=-cbk(k-1+1)-k./z.*cbk(k+1);
end;  k=fix(nm)+1;
return;
end
function [z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1]=cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1,varargin);
%       ==========================================================
%       Purpose: Compute modified complex Bessel functions I0(z),
%                I1(z), K0(z), K1(z), and their derivatives
%       Input :  z --- Complex argument
%       Output:  CBI0 --- I0(z)
%                CDI0 --- I0'(z)
%                CBI1 --- I1(z)
%                CDI1 --- I1'(z)
%                CBK0 --- K0(z)
%                CDK0 --- K0'(z)
%                CBK1 --- K1(z)
%                CDK1 --- K1'(z)
%       ==========================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,10);
cw=0.0;
pi=3.141592653589793d0;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z2=z.*z;
z1=z;
if(a0 == 0.0d0);
cbi0=complex(1.0d0,0.0d0);
cbi1=complex(0.0d0,0.0d0);
cdi0=complex(0.0d0,0.0d0);
cdi1=complex(0.5d0,0.0d0);
cbk0=complex(1.0d+300,0.0d0);
cbk1=complex(1.0d+300,0.0d0);
cdk0=-complex(1.0d+300,0.0d0);
cdk1=-complex(1.0d+300,0.0d0);
return;
end;
if(real(z)< 0.0)z1=-z; end;
if(a0 <= 18.0);
cbi0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:50;
cr=0.25d0.*cr.*z2./(k.*k);
cbi0=cbi0+cr;
if(abs(cr./cbi0)< 1.0d-15)break; end;
end;
cbi1=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:50;
cr=0.25d0.*cr.*z2./(k.*(k));
cbi1=cbi1+cr;
if(abs(cr./cbi1)< 1.0d-15)break; end;
end;
cbi1=0.5d0.*z1.*cbi1;
else;
a(:)=[0.125d0,7.03125d-2,7.32421875d-2,1.1215209960938d-1,2.2710800170898d-1,5.7250142097473d-1,1.7277275025845d0,6.0740420012735d0,2.4380529699556d01,1.1001714026925d02,5.5133589612202d02,3.0380905109224d03];
b(:)=[-0.375d0,-1.171875d-1,-1.025390625d-1,-1.4419555664063d-1,-2.7757644653320d-1,-6.7659258842468d-1,-1.9935317337513d0,-6.8839142681099d0,-2.7248827311269d01,-1.2159789187654d02,-6.0384407670507d02,-3.3022722944809d03];
k0=12;
if(a0 >= 35.0)k0=9; end;
if(a0 >= 50.0)k0=7; end;
ca=exp(z1)./sqrt(2.0d0.*pi.*z1);
cbi0=complex(1.0d0,0.0d0);
zr=1.0d0./z1;
for  k=1:k0;
cbi0=cbi0+a(k).*zr.^k;
end;  k=k0+1;
cbi0=ca.*cbi0;
cbi1=complex(1.0d0,0.0d0);
for  k=1:k0;
cbi1=cbi1+b(k).*zr.^k;
end;  k=k0+1;
cbi1=ca.*cbi1;
end;
if(a0 <= 9.0);
cs=complex(0.0d0,0.0d0);
ct=-log(0.5d0.*z1)-0.5772156649015329d0;
w0=0.0d0;
cr=complex(1.0d0,0.0d0);
for  k=1:50;
w0=w0+1.0d0./k;
cr=0.25d0.*cr./(k.*k).*z2;
cs=cs+cr.*(w0+ct);
if(abs((cs-cw)./cs)< 1.0d-15)break; end;
cw=cs;
end;
cbk0=ct+cs;
else;
a1(:)=[0.125d0,0.2109375d0,1.0986328125d0,1.1775970458984d01,2.1461706161499d02,5.9511522710323d03,2.3347645606175d05,1.2312234987631d07,8.401390346421d08,7.2031420482627d10];
cb=0.5d0./z1;
zr2=1.0d0./z2;
cbk0=complex(1.0d0,0.0d0);
for  k=1:10;
cbk0=cbk0+a1(k).*zr2.^k;
end;  k=10+1;
cbk0=cb.*cbk0./cbi0;
end;
cbk1=(1.0d0./z1-cbi1.*cbk0)./cbi0;
if(real(z)< 0.0);
if(imag(z)< 0.0)cbk0=cbk0+ci.*pi.*cbi0; end;
if(imag(z)> 0.0)cbk0=cbk0-ci.*pi.*cbi0; end;
if(imag(z)< 0.0)cbk1=-cbk1+ci.*pi.*cbi1; end;
if(imag(z)> 0.0)cbk1=-cbk1-ci.*pi.*cbi1; end;
cbi1=-cbi1;
end;
cdi0=cbi1;
cdi1=cbi0-1.0d0./z.*cbi1;
cdk0=-cbk1;
cdk1=-cbk0-1.0d0./z.*cbk1;
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

