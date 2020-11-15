function mcjyna
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ================================================================
%     Purpose: This program computes Bessel functions Jn(z), Yn(z)
%     and their derivatives for a complex argument using
%     subroutine CJYNA
%     Input :  z --- Complex argument of Jn(z)and Yn(z)
%     n --- Order of Jn(z)and Yn(z)
%(n = 0,1,תתת, n ף 250)
%     Output:  CBJ(n)--- Jn(z)
%     CDJ(n)--- Jn'(z)
%     CBY(n)--- Yn(z)
%     CDY(n)--- Yn'(z)
%     Eaxmple: z = 4.0 + i 2.0
%     n     Re[Jn(z)]Im[Jn(z)]Re[Jn'(z)]Im[Jn'(z)]
%     ------------------------------------------------------------------
%     -
%     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
%     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
%     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
%     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
%     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
%     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00
%     n     Re[Yn(z)]Im[Yn(z)]Re[Yn'(z)]Im[Yn'(z)]
%     ------------------------------------------------------------------
%     --
%     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
%     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
%     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
%     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
%     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
%     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01
%     z = 20.0 + i 10.0 ,      Nmax = 5
%     n     Re[Jn(z)]Im[Jn(z)]Re[Jn'(z)]Im[Jn'(z)]
%     ------------------------------------------------------------------
%     --
%     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
%     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
%     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
%     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
%     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
%     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04
%     n     Re[Yn(z)]Im[Yn(z)]Re[Yn'(z)]Im[Yn'(z)]
%     ------------------------------------------------------------------
%     --
%     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
%     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
%     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
%     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
%     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
%     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
%     ================================================================
n=[];z=[];nm=[];
 global cbj;
global cdj;
global cby;
global cdy;
fprintf(1,'%s \n','  please enter n, x,y(z=x+iy)');
%     READ(*,*)N,X,Y
n=5;
x=4.0;
y=2.0;
z=complex(x,y);
fprintf(1,[repmat(' ',1,3),'z =','%5.1g',' + i ','%5.1g',' ,',repmat(' ',1,6),'nmax =','%3g' ' \n'],x,y,n);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%       READ(*,*)NS
ns=1;
end;
[n,z,nm,cbj,cdj,cby,cdy]=cjyna(n,z,nm,cbj,cdj,cby,cdy);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n     re[jn(z)]im[jn(z)]');fprintf(1,'%s \n','       re[jn''(z)]im[jn''(z)]');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','-------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat('%16.8g',1,4) ' \n'],k,cbj(k+1),cdj(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n     re[yn(z)]im[yn(z)]');fprintf(1,'%s \n','       re[yn''(z)]im[yn''(z)]');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','-------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat('%16.8g',1,4) ' \n'],k,cby(k+1),cdy(k+1));
end;  k=nm+1;
%format(1x,i4,4d16.8);
%format(3x,,f5.1,' + i ',f5.1,' ,',6x,,i3);
end
function [n,z,nm,cbj,cdj,cby,cdy]=cjyna(n,z,nm,cbj,cdj,cby,cdy,varargin);
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
%     Rouitines called:
%(1)CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
%(2)MSTA1 and MSTA2 to calculate the starting
%     point for backward recurrence
%     =======================================================
cbj0=[];cdj0=[];cbj1=[];cdj1=[];cby0=[];cdy0=[];cby1=[];cdy1=[];
lb0=0.0;
pi=3.141592653589793d0;
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
[z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1]=cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1);
cbj(0+1)=cbj0;
cbj(1+1)=cbj1;
cby(0+1)=cby0;
cby(1+1)=cby1;
cdj(0+1)=cdj0;
cdj(1+1)=cdj1;
cdy(0+1)=cdy0;
cdy(1+1)=cdy1;
if(n <= 1)return; end;
if(n < fix(0.25.*a0));
cj0=cbj0;
cj1=cbj1;
for  k=2:n;
cjk=2.0d0.*(k-1.0d0)./z.*cj1-cj0;
cbj(k+1)=cjk;
cj0=cj1;
cj1=cjk;
end;  k=fix(n)+1;
else;
m=msta1(a0,200);
if(m < n);
nm=fix(m);
else;
m=msta2(a0,fix(n),15);
end;
cf2=complex(0.0d0,0.0d0);
cf1=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=2.0d0.*(k+1.0d0)./z.*cf1-cf2;
if(k <= nm)cbj(k+1)=cf; end;
cf2=cf1;
cf1=cf;
end;  k=0-1;
if(abs(cbj0)> abs(cbj1));
cs=cbj0./cf;
else;
cs=cbj1./cf2;
end;
for  k=0:nm;
cbj(k+1)=cs.*cbj(k+1);
end;  k=fix(nm)+1;
end;
for  k=2:nm;
cdj(k+1)=cbj(k-1+1)-k./z.*cbj(k+1);
end;  k=fix(nm)+1;
ya0=abs(cby0);
lb=0;
cg0=cby0;
cg1=cby1;
for  k=2:nm;
cyk=2.0d0.*(k-1.0d0)./z.*cg1-cg0;
if(~(abs(cyk)> 1.0d+290));
yak=abs(cyk);
ya1=abs(cg0);
if(yak < ya0&yak < ya1)lb=k; end;
cby(k+1)=cyk;
cg0=cg1;
cg1=cyk;
end;
end;  k=fix(nm)+1;
while (1);
if(lb <= 4|imag(z)== 0.0d0)break; end;
if(lb == lb0)break; end;
ch2=complex(1.0d0,0.0d0);
ch1=complex(0.0d0,0.0d0);
lb0=lb;
for  k=lb:-1:1;
ch0=2.0d0.*k./z.*ch1-ch2;
ch2=ch1;
ch1=ch0;
end;  k=1-1;
cp12=ch0;
cp22=ch2;
ch2=complex(0.0d0,0.0d0);
ch1=complex(1.0d0,0.0d0);
for  k=lb:-1:1;
ch0=2.0d0.*k./z.*ch1-ch2;
ch2=ch1;
ch1=ch0;
end;  k=1-1;
cp11=ch0;
cp21=ch2;
if(lb == nm)cbj(lb+1+1)=2.0d0.*lb./z.*cbj(lb+1)-cbj(lb-1+1); end;
if(abs(cbj(0+1))> abs(cbj(1+1)));
cby(lb+1+1)=(cbj(lb+1+1).*cby0-2.0d0.*cp11./(pi.*z))./cbj(0+1);
cby(lb+1)=(cbj(lb+1).*cby0+2.0d0.*cp12./(pi.*z))./cbj(0+1);
else;
cby(lb+1+1)=(cbj(lb+1+1).*cby1-2.0d0.*cp21./(pi.*z))./cbj(1+1);
cby(lb+1)=(cbj(lb+1).*cby1+2.0d0.*cp22./(pi.*z))./cbj(1+1);
end;
cyl2=cby(lb+1+1);
cyl1=cby(lb+1);
for  k=lb-1:-1:0;
cylk=2.0d0.*(k+1.0d0)./z.*cyl1-cyl2;
cby(k+1)=cylk;
cyl2=cyl1;
cyl1=cylk;
end;  k=0-1;
cyl1=cby(lb+1);
cyl2=cby(lb+1+1);
for  k=lb+1:nm-1;
cylk=2.0d0.*k./z.*cyl2-cyl1;
cby(k+1+1)=cylk;
cyl1=cyl2;
cyl2=cylk;
end;  k=fix(nm)-1+1;
for  k=2:nm;
wa=abs(cby(k+1));
if(wa < abs(cby(k-1+1)))lb=k; end;
end;  k=fix(nm)+1;
end;
for  k=2:nm;
cdy(k+1)=cby(k-1+1)-k./z.*cby(k+1);
end;  k=fix(nm)+1;
return;
end
function [z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1]=cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1,varargin);
%     ===========================================================
%     Purpose: Compute complex Bessel functions J0(z), J1(z)
%     Y0(z), Y1(z), and their derivatives
%     Input :  z --- Complex argument
%     Output:  CBJ0 --- J0(z)
%     CDJ0 --- J0'(z)
%     CBJ1 --- J1(z)
%     CDJ1 --- J1'(z)
%     CBY0 --- Y0(z)
%     CDY0 --- Y0'(z)
%     CBY1 --- Y1(z)
%     CDY1 --- Y1'(z)
%     ===========================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,12);
b1=zeros(1,12);
pi=3.141592653589793d0;
el=0.5772156649015329d0;
rp2=2.0d0./pi;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z2=z.*z;
z1=z;
if(a0 == 0.0d0);
cbj0=complex(1.0d0,0.0d0);
cbj1=complex(0.0d0,0.0d0);
cdj0=complex(0.0d0,0.0d0);
cdj1=complex(0.5d0,0.0d0);
cby0=-complex(1.0d300,0.0d0);
cby1=-complex(1.0d300,0.0d0);
cdy0=complex(1.0d300,0.0d0);
cdy1=complex(1.0d300,0.0d0);
return;
end;
if(real(z)< 0.0)z1=-z; end;
if(a0 <= 12.0);
cbj0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*k);
cbj0=cbj0+cr;
if(abs(cr./cbj0)< 1.0d-15)break; end;
end;
cbj1=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*(k+1.0d0));
cbj1=cbj1+cr;
if(abs(cr./cbj1)< 1.0d-15)break; end;
end;
cbj1=0.5d0.*z1.*cbj1;
w0=0.0d0;
cr=complex(1.0d0,0.0d0);
cs=complex(0.0d0,0.0d0);
for  k=1:40;
w0=w0+1.0d0./k;
cr=-0.25d0.*cr./(k.*k).*z2;
cp=cr.*w0;
cs=cs+cp;
if(abs(cp./cs)< 1.0d-15)break; end;
end;
cby0=rp2.*(log(z1./2.0d0)+el).*cbj0-rp2.*cs;
w1=0.0d0;
cr=complex(1.0d0,0.0d0);
cs=complex(1.0d0,0.0d0);
for  k=1:40;
w1=w1+1.0d0./k;
cr=-0.25d0.*cr./(k.*(k)).*z2;
cp=cr.*(2.0d0.*w1+1.0d0./(k+1.0d0));
cs=cs+cp;
if(abs(cp./cs)< 1.0d-15)break; end;
end;
cby1=rp2.*((log(z1./2.0d0)+el).*cbj1-1.0d0./z1-.25d0.*z1.*cs);
else;
a(:)=[-.703125d-01,.112152099609375d+00,-.5725014209747314d+00,.6074042001273483d+01,-.1100171402692467d+03,.3038090510922384d+04,-.1188384262567832d+06,.6252951493434797d+07,-.4259392165047669d+09,.3646840080706556d+11,-.3833534661393944d+13,.4854014686852901d+15];
b(:)=[.732421875d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02,.5513358961220206d+03,-.1825775547429318d+05,.8328593040162893d+06,-.5006958953198893d+08,.3836255180230433d+10,-.3649010818849833d+12,.4218971570284096d+14,-.5827244631566907d+16];
a1(:)=[.1171875d+00,-.144195556640625d+00,.6765925884246826d+00,-.6883914268109947d+01,.1215978918765359d+03,-.3302272294480852d+04,.1276412726461746d+06,-.6656367718817688d+07,.4502786003050393d+09,-.3833857520742790d+11,.4011838599133198d+13,-.5060568503314727d+15];
b1(:)=[-.1025390625d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02,-.6038440767050702d+03,.1971837591223663d+05,-.8902978767070678d+06,.5310411010968522d+08,-.4043620325107754d+10,.3827011346598605d+12,-.4406481417852278d+14,.6065091351222699d+16];
k0=12;
if(a0 >= 35.0)k0=10; end;
if(a0 >= 50.0)k0=8; end;
ct1=z1-0.25d0.*pi;
cp0=complex(1.0d0,0.0d0);
for  k=1:k0;
cp0=cp0+a(k).*z1.^(-2.*k);
end;  k=k0+1;
cq0=-0.125d0./z1;
for  k=1:k0;
cq0=cq0+b(k).*z1.^(-2.*k-1);
end;  k=k0+1;
cu=sqrt(rp2./z1);
cbj0=cu.*(cp0.*cos(ct1)-cq0.*sin(ct1));
cby0=cu.*(cp0.*sin(ct1)+cq0.*cos(ct1));
ct2=z1-0.75d0.*pi;
cp1=complex(1.0d0,0.0d0);
for  k=1:k0;
cp1=cp1+a1(k).*z1.^(-2.*k);
end;  k=k0+1;
cq1=0.375d0./z1;
for  k=1:k0;
cq1=cq1+b1(k).*z1.^(-2.*k-1);
end;  k=k0+1;
cbj1=cu.*(cp1.*cos(ct2)-cq1.*sin(ct2));
cby1=cu.*(cp1.*sin(ct2)+cq1.*cos(ct2));
end;
if(real(z)< 0.0);
if(imag(z)< 0.0)cby0=cby0-2.0d0.*ci.*cbj0; end;
if(imag(z)> 0.0)cby0=cby0+2.0d0.*ci.*cbj0; end;
if(imag(z)< 0.0)cby1=-(cby1-2.0d0.*ci.*cbj1); end;
if(imag(z)> 0.0)cby1=-(cby1+2.0d0.*ci.*cbj1); end;
cbj1=-cbj1;
end;
cdj0=-cbj1;
cdj1=cbj0-1.0d0./z.*cbj1;
cdy0=-cby1;
cdy1=cby0-1.0d0./z.*cby1;
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
msta1Result=fix(nn);
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
msta2Result=fix(nn+10);
return;
end
function [envjResult]=envj(n,x,varargin);
envjResult=0.5d0.*log10(6.28d0.*fix(n))-fix(n).*log10(1.36d0.*x./fix(n));
return;
end

