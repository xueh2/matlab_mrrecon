function mcsphik
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes the modified spherical Bessel
%                functions and their derivatives for a complex
%                argument using subroutine CSPHIK
%       Input :  z --- Complex argument
%                n --- Order of in(z)& kn(z,0 ó n ó 250)
%       Output:  CSI(n)--- in(z)
%                CDI(n)--- in'(z)
%                CSK(n)--- kn(z)
%                CDK(n)--- kn'(z)
%       Example: z =4.0+i 2.0
%     n     Re[in(z)]Im[in(z)]Re[in'(z)]Im[in'(z)]
%    ---------------------------------------------------------------
%     0   .2118080D+00   .6101922D+01  -.4439356D+00   .4900150D+01
%     1  -.4439356D+00   .4900150D+01  -.5906477D+00   .4053075D+01
%     2  -.9918756D+00   .3028652D+01  -.7574058D+00   .2785396D+01
%     3  -.9663859D+00   .1375561D+01  -.7689911D+00   .1541649D+01
%     4  -.6018277D+00   .4263967D+00  -.5777565D+00   .6482500D+00
%     5  -.2668530D+00   .6640148D-01  -.3214450D+00   .1866032D+00
%     n     Re[kn(z)]Im[kn(z)]Re[kn'(z)]Im[kn'(z)]
%    ---------------------------------------------------------------
%     0  -.5010582D-02  -.4034862D-02   .6416184D-02   .4340777D-02
%     1  -.6416184D-02  -.4340777D-02   .8445211D-02   .4487936D-02
%     2  -.1016253D-01  -.4714473D-02   .1392804D-01   .4120703D-02
%     3  -.1893595D-01  -.3973987D-02   .2690088D-01   .3192843D-03
%     4  -.3945464D-01   .2977107D-02   .5690203D-01  -.1873044D-01
%     5  -.8727490D-01   .3689398D-01   .1220481D+00  -.9961483D-01
%       =============================================================
n=[];z=[];nm=[];csi=[];cdi=[];csk=[];cdk=[];
  x=0;
 y=0;
 csi=zeros(1,250+1);
cdi=zeros(1,250+1);
csk=zeros(1,250+1);
cdk=zeros(1,250+1);
fprintf(1,'%s \n','please enter n,x,y(z=x+iy)');
%        READ(*,*)N,X,Y
n=5;
x=4.0;
y=2.0;
fprintf(1,[repmat(' ',1,3),'nmaz =','%3g',',     ','z = ','%8.1g','+ i','%8.1g' ' \n'],n,x,y);
z=complex(x,y);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','please enter order step ns ');
%           READ(*,*)NS
ns=1;
end;
[n,z,nm,csi,cdi,csk,cdk]=csphik(n,z,nm,csi,cdi,csk,cdk);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[in(z)]im[in(z)]');fprintf(1,'%s \n','        re[in''(z)]im[in''(z)]');
fprintf(1,'%s ','--------------------------------------------');fprintf(1,'%s \n','----------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,4) ' \n'],k,csi(k+1),cdi(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[kn(z)]im[kn(z)]');fprintf(1,'%s \n','        re[kn''(z)]im[kn''(z)]');
fprintf(1,'%s ','--------------------------------------------');fprintf(1,'%s \n','----------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,4) ' \n'],k,csk(k+1),cdk(k+1));
end;  k=nm+1;
%format(1x,i3,4d17.8);
%format(3x,',i3,',     ','z = ',f8.1,'+ i',f8.1);
end
function [n,z,nm,csi,cdi,csk,cdk]=csphik(n,z,nm,csi,cdi,csk,cdk,varargin);
%       =======================================================
%       Purpose: Compute modified spherical Bessel functions
%                and their derivatives for a complex argument
%       Input :  z --- Complex argument
%                n --- Order of in(z)& kn(z,n = 0,1,2,...)
%       Output:  CSI(n)--- in(z)
%                CDI(n)--- in'(z)
%                CSK(n)--- kn(z)
%                CDK(n)--- kn'(z)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       =======================================================
  a0=0;
 pi=0;
pi=3.141592653589793d0;
a0=abs(z);
nm=fix(fix(n));
if(a0 < 1.0d-60);
for  k=0:n;
csi(k+1)=0.0d0;
cdi(k+1)=0.0d0;
csk(k+1)=1.0d+300;
cdk(k+1)=-1.0d+300;
end;  k=fix(n)+1;
csi(0+1)=1.0d0;
cdi(1+1)=0.3333333333333333d0;
return;
end;
ci=complex(0.0d0,1.0d0);
csinh=sin(ci.*z)./ci;
ccosh=cos(ci.*z);
csi0=csinh./z;
csi1=(-csinh./z+ccosh)./z;
csi(0+1)=csi0;
csi(1+1)=csi1;
if(n >= 2);
m=msta1(a0,200);
if(m < n);
nm=fix(m);
else;
m=msta2(a0,fix(n),15);
end;
cf0=0.0d0;
cf1=1.0d0-100;
for  k=m:-1:0;
cf=(2.0d0.*k+3.0d0).*cf1./z+cf0;
if(k <= nm)csi(k+1)=cf; end;
cf0=cf1;
cf1=cf;
end;  k=0-1;
if(abs(csi0)> abs(csi1))cs=csi0./cf; end;
if(abs(csi0)<= abs(csi1))cs=csi1./cf0; end;
for  k=0:nm;
csi(k+1)=cs.*csi(k+1);
end;  k=fix(nm)+1;
end;
cdi(0+1)=csi(1+1);
for  k=1:nm;
cdi(k+1)=csi(k-1+1)-(k+1.0d0).*csi(k+1)./z;
end;  k=fix(nm)+1;
csk(0+1)=0.5d0.*pi./z.*exp(-z);
csk(1+1)=csk(0+1).*(1.0d0+1.0d0./z);
for  k=2:nm;
if(abs(csi(k-1+1))> abs(csi(k-2+1)));
csk(k+1)=(0.5d0.*pi./(z.*z)-csi(k+1).*csk(k-1+1))./csi(k-1+1);
else;
csk(k+1)=(csi(k+1).*csk(k-2+1)+(k-0.5d0).*pi./z.^3)./csi(k-2+1);
end;
end;  k=fix(nm)+1;
cdk(0+1)=-csk(1+1);
for  k=1:nm;
cdk(k+1)=-csk(k-1+1)-(k+1.0d0).*csk(k+1)./z;
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

