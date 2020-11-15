function mcjynb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ================================================================
%       Purpose: This program computes Bessel functions Jn(z), Yn(z)
%                and their derivatives for a complex argument using
%                subroutine CJYNB
%       Input :  z --- Complex argument of Jn(z)and Yn(z)
%                n --- Order of Jn(z)and Yn(z)
%(n = 0,1,תתת, n ף 250)
%       Output:  CBJ(n)--- Jn(z)
%                CDJ(n)--- Jn'(z)
%                CBY(n)--- Yn(z)
%                CDY(n)--- Yn'(z)
%       Eaxmple: z = 4.0 + i 2.0
%     n     Re[Jn(z)]Im[Jn(z)]Re[Jn'(z)]Im[Jn'(z)]
%    -------------------------------------------------------------------
%     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
%     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
%     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
%     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
%     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
%     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00
%     n     Re[Yn(z)]Im[Yn(z)]Re[Yn'(z)]Im[Yn'(z)]
%   --------------------------------------------------------------------
%     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
%     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
%     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
%     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
%     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
%     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01
%                z = 20.0 + i  10.0 ,      Nmax =  5
%     n     Re[Jn(z)]Im[Jn(z)]Re[Jn'(z)]Im[Jn'(z)]
%   --------------------------------------------------------------------
%     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
%     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
%     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
%     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
%     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
%     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04
%     n     Re[Yn(z)]Im[Yn(z)]Re[Yn'(z)]Im[Yn'(z)]
%   --------------------------------------------------------------------
%     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
%     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
%     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
%     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
%     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
%     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
%       ================================================================
n=[];z=[];nm=[];
 global cbj;
global cdj;
global cby;
global cdy;
fprintf(1,'%s \n','  please enter n, x,y(z=x+iy)');
%        READ(*,*)N,X,Y
n=3;
x=4.0;
y=2.0;
z=complex(x,y);
fprintf(1,[repmat(' ',1,3),'z =','%5.1g',' + i ','%5.1g',' ,',repmat(' ',1,6),'nmax =','%3g' ' \n'],x,y,n);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,z,nm,cbj,cdj,cby,cdy]=cjynb(n,z,nm,cbj,cdj,cby,cdy);
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
function [n,z,nm,cbj,cdj,cby,cdy]=cjynb(n,z,nm,cbj,cdj,cby,cdy,varargin);
%       =======================================================
%       Purpose: Compute Bessel functions Jn(z), Yn(z)and
%                their derivatives for a complex argument
%       Input :  z --- Complex argument of Jn(z)and Yn(z)
%                n --- Order of Jn(z)and Yn(z)
%       Output:  CBJ(n)--- Jn(z)
%                CDJ(n)--- Jn'(z)
%                CBY(n)--- Yn(z)
%                CDY(n)--- Yn'(z)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 to calculate the starting
%                point for backward recurrence
%       =======================================================
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
if(a0 <= 300.d0|n > fix(0.25.*a0));
if(n == 0)nm=1; end;
m=msta1(a0,200);
if(m < nm);
nm=fix(m);
else;
m=msta2(a0,fix(nm),15);
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
cyy=(cbj(k+1).*cby(k-2+1)-4.0d0.*(k-1.0d0)./(pi.*z.*z))./cbj(k-2+1+1);
end;
cby(k+1)=cyy;
end;  k=fix(nm)+1;
cdy(0+1)=-cby(1+1);
for  k=1:nm;
cdy(k+1)=cby(k-1+1)-k./z.*cby(k+1);
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

