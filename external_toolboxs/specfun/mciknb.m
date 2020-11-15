function mciknb
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
%                complex argument using subroutine CIKNB
%       Input:   z --- Complex argument
%                n --- Order of In(z)and Kn(z)
%(n = 0,1,תתת, n ף 250)
%       Output:  CBI(n)--- In(z)
%                CDI(n)--- In'(z)
%                CBK(n)--- Kn(z)
%                CDK(n)--- Kn'(z)
%       Example: Nmax = 5,   z = 4.0 + i 2.0
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
n=[];z=[];nm=[];cbi=[];cdi=[];cbk=[];cdk=[];
 cbi=zeros(1,250+1);
cdi=zeros(1,250+1);
cbk=zeros(1,250+1);
cdk=zeros(1,250+1);
fprintf(1,'%s \n','  please enter n, x and y(z = x+iy)');
%        READ(*,*)N,X,Y
n=5;
x=4.0;
y=2.0;
fprintf(1,[repmat(' ',1,3),'nmaz =','%3g',',    ','z =','%6.1g',' + i','%6.1g' ' \n'],n,x,y);
z=complex(x,y);
fprintf(1,'%0.15g \n');
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,z,nm,cbi,cdi,cbk,cdk]=ciknb(n,z,nm,cbi,cdi,cbk,cdk);
fprintf(1,'%s ','   n      re[in(z)]im[in(z)]');fprintf(1,'%s \n','      re[in''(z)]im[in''(z)]');
fprintf(1,'%s ',' ---------------------------------');fprintf(1,'%s \n','------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),repmat(' ',1,1),'%3g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,cbi(k+1),cdi(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n      re[kn(z)]im[kn(z)]');fprintf(1,'%s \n','      re[kn''(z)]im[kn''(z)]');
fprintf(1,'%s ',' ---------------------------------');fprintf(1,'%s \n','------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),repmat(' ',1,1),'%3g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,cbk(k+1),cdk(k+1));
end;  k=nm+1;
%format(1x,1x,i3,1x,4d16.8);
%format(3x,,i3,',    ',', f6.1,' + i',f6.1);
end
function [n,z,nm,cbi,cdi,cbk,cdk]=ciknb(n,z,nm,cbi,cdi,cbk,cdk,varargin);
%       ============================================================
%       Purpose: Compute modified Bessel functions In(z)and Kn(z),
%                and their derivatives for a complex argument
%       Input:   z --- Complex argument
%                n --- Order of In(z)and Kn(z)
%       Output:  CBI(n)--- In(z)
%                CDI(n)--- In'(z)
%                CBK(n)--- Kn(z)
%                CDK(n)--- Kn'(z)
%                NM --- Highest order computed
%       Routones called:
%                MSTA1 and MSTA2 to compute the starting point for
%                backward recurrence
%       ===========================================================
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
m=msta1(a0,200);
if(m < nm);
nm=fix(m);
else;
m=msta2(a0,fix(nm),15);
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

