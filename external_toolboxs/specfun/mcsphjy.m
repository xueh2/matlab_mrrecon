function mcsphjy
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ================================================================
%       Purpose: This program computes the spherical Bessel functions
%                jn(z), yn(z), and their derivatives for a complex
%                argument using subroutine CSPHJY
%       Input :  z --- Complex argument
%                n --- Order of jn(z)& yn(z,0 ó n ó 250)
%       Output:  CSJ(n)--- jn(z)
%                CDJ(n)--- jn'(z)
%                CSY(n)--- yn(z)
%                CDY(n)--- yn'(z)
%       Example: z = 4.0+i 2.0
%     n     Re[jn(z)]Im[jn(z)]Re[jn'(z)]Im[jn'(z)]
%   --------------------------------------------------------------------
%     0  -.80651523D+00  -.18941093D+00  -.37101203D-01   .75210758D+00
%     1   .37101203D-01  -.75210758D+00  -.67093420D+00   .11885235D+00
%     2   .60314368D+00  -.27298399D+00  -.24288981D+00  -.40737409D+00
%     3   .42955048D+00   .17755176D+00   .18848259D+00  -.24320520D+00
%     4   .12251323D+00   .22087111D+00   .19660170D+00   .17937264D-01
%     5  -.10242676D-01   .10975433D+00   .68951842D-01   .83020305D-01
%     n     Re[yn(z)]Im[yn(z)]Re[yn'(z)]Im[yn'(z)]
%   --------------------------------------------------------------------
%     0   .21734534D+00  -.79487692D+00  -.77049661D+00  -.87010064D-02
%     1   .77049661D+00   .87010064D-02  -.92593503D-01  -.64425800D+00
%     2   .24756293D+00   .56894854D+00   .45127429D+00  -.25839924D+00
%     3  -.23845941D+00   .43646607D+00   .26374403D+00   .12439192D+00
%     4  -.27587985D+00   .20902555D+00  -.67092335D-01   .89500599D-01
%     5  -.70001327D-01   .18807178D+00  -.30472133D+00  -.58661384D-01
%       ================================================================
n=[];z=[];nm=[];csj=[];cdj=[];csy=[];cdy=[];
  x=0;
 y=0;
 csj=zeros(1,250+1);
cdj=zeros(1,250+1);
csy=zeros(1,250+1);
cdy=zeros(1,250+1);
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
[n,z,nm,csj,cdj,csy,cdy]=csphjy(n,z,nm,csj,cdj,csy,cdy);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[jn(z)]im[jn(z)]');fprintf(1,'%s \n','        re[jn''(z)]im[jn''(z)]');
fprintf(1,'%s ','--------------------------------------------');fprintf(1,'%s \n','----------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,4) ' \n'],k,csj(k+1),cdj(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[yn(z)]im[yn(z)]');fprintf(1,'%s \n','        re[yn''(z)]im[yn''(z)]');
fprintf(1,'%s ','--------------------------------------------');fprintf(1,'%s \n','----------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,4) ' \n'],k,csy(k+1),cdy(k+1));
end;  k=nm+1;
%format(1x,i3,4d17.8);
%format(3x,,i3,',     ','z = ',f8.1,'+ i',f8.1);
end
function [n,z,nm,csj,cdj,csy,cdy]=csphjy(n,z,nm,csj,cdj,csy,cdy,varargin);
%       ==========================================================
%       Purpose: Compute spherical Bessel functions jn(z)& yn(z)
%                and their derivatives for a complex argument
%       Input :  z --- Complex argument
%                n --- Order of jn(z)& yn(z,n = 0,1,2,...)
%       Output:  CSJ(n)--- jn(z)
%                CDJ(n)--- jn'(z)
%                CSY(n)--- yn(z)
%                CDY(n)--- yn'(z)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       ==========================================================
  a0=0;
a0=abs(z);
nm=fix(fix(n));
if(a0 < 1.0d-60);
for  k=0:n;
csj(k+1)=0.0d0;
cdj(k+1)=0.0d0;
csy(k+1)=-1.0d+300;
cdy(k+1)=1.0d+300;
end;  k=fix(n)+1;
csj(0+1)=complex(1.0d0,0.0d0);
cdj(1+1)=complex(.333333333333333d0,0.0d0);
return;
end;
csj(0+1)=sin(z)./z;
csj(1+1)=(csj(0+1)-cos(z))./z;
if(n >= 2);
csa=csj(0+1);
csb=csj(1+1);
m=msta1(a0,200);
if(m < n);
nm=fix(m);
else;
m=msta2(a0,fix(n),15);
end;
cf0=0.0d0;
cf1=1.0d0-100;
for  k=m:-1:0;
cf=(2.0d0.*k+3.0d0).*cf1./z-cf0;
if(k <= nm)csj(k+1)=cf; end;
cf0=cf1;
cf1=cf;
end;  k=0-1;
if(abs(csa)> abs(csb))cs=csa./cf; end;
if(abs(csa)<= abs(csb))cs=csb./cf0; end;
for  k=0:nm;
csj(k+1)=cs.*csj(k+1);
end;  k=fix(nm)+1;
end;
cdj(0+1)=(cos(z)-sin(z)./z)./z;
for  k=1:nm;
cdj(k+1)=csj(k-1+1)-(k+1.0d0).*csj(k+1)./z;
end;  k=fix(nm)+1;
csy(0+1)=-cos(z)./z;
csy(1+1)=(csy(0+1)-sin(z))./z;
cdy(0+1)=(sin(z)+cos(z)./z)./z;
cdy(1+1)=(2.0d0.*cdy(0+1)-cos(z))./z;
for  k=2:nm;
if(abs(csj(k-1+1))> abs(csj(k-2+1)));
csy(k+1)=(csj(k+1).*csy(k-1+1)-1.0d0./(z.*z))./csj(k-1+1);
else;
csy(k+1)=(csj(k+1).*csy(k-2+1)-(2.0d0.*k-1.0d0)./z.^3)./csj(k-2+1);
end;
end;  k=fix(nm)+1;
for  k=2:nm;
cdy(k+1)=csy(k-1+1)-(k+1.0d0).*csy(k+1)./z;
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

