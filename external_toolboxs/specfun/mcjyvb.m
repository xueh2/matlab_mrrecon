function mcjyvb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes Bessel functions Jv(z), Yv(z),
%                and their derivatives for a complex argument using
%                subroutine CJYVB
%       Input :  z --- Complex argument
%                v --- Order of Jv(z)and Yv(z)
%(v = n+v0, 0 ף n ף 250, 0 ף v0 < 1)
%       Output:  CBJ(n)--- Jn+v0(z)
%                CDJ(n)--- Jn+v0'(z)
%                CBY(n)--- Yn+v0(z)
%                CDY(n)--- Yn+v0'(z)
%       Example:
%                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0
%     n     Re[Jv(z)]Im[Jv(z)]Re[Jv'(z)]Im[Jv'(z)]
%    -------------------------------------------------------------------
%     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
%     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
%     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
%     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
%     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
%     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00
%     n     Re[Yv(z)]Im[Yv(z)]Re[Yv'(z)]Im[Yv'(z)]
%    -------------------------------------------------------------------
%     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
%     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
%     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
%     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
%     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
%     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
%       =============================================================
v=[];z=[];vm=[];cbj=[];cdj=[];cby=[];cdy=[];
 cbj=zeros(1,250+1);
cdj=zeros(1,250+1);
cby=zeros(1,250+1);
cdy=zeros(1,250+1);
fprintf(1,'%s \n','  please enter v, x and y(z=x+iy)');
%        READ(*,*)V,X,Y
v=5;
x=4.0;
y=2.0;
z=complex(x,y);
n=fix(v);
v0=v-n;
fprintf(1,[repmat(' ',1,8),'v = n+v0',',  v0 =','%5.2g',',  z =','%7.2g',' +','%7.2g','i' ' \n'],v0,x,y);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[v,z,vm,cbj,cdj,cby,cdy]=cjyvb(v,z,vm,cbj,cdj,cby,cdy);
nm=fix(vm);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n       re[jv(z)]im[jv(z)]');fprintf(1,'%s \n','       re[jv''(z)]im[jv''(z)]');
fprintf(1,'%s ',' ----------------------------------');fprintf(1,'%s \n','-----------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,2),repmat('%16.8g',1,4) ' \n'],k,cbj(k+1),cdj(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n       re[yv(z)]im[yv(z)]');fprintf(1,'%s \n','       re[yv''(z)]im[yv''(z)]');
fprintf(1,'%s ',' ----------------------------------');fprintf(1,'%s \n','-----------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,2),repmat('%16.8g',1,4) ' \n'],k,cby(k+1),cdy(k+1));
end;  k=nm+1;
%format(1x,i3,2x,4d16.8);
%format(8x,'v = n+v0',',f5.2,',f7.2,' +',f7.2,'i');
end
function [v,z,vm,cbj,cdj,cby,cdy]=cjyvb(v,z,vm,cbj,cdj,cby,cdy,varargin);
%       ===========================================================
%       Purpose: Compute Bessel functions Jv(z), Yv(z)and their
%                derivatives for a complex argument
%       Input :  z --- Complex argument
%                v --- Order of Jv(z)and Yv(z)
%(v = n+v0, n = 0,1,2,..., 0 ף v0 < 1)
%       Output:  CBJ(n)--- Jn+v0(z)
%                CDJ(n)--- Jn+v0'(z)
%                CBY(n)--- Yn+v0(z)
%                CDY(n)--- Yn+v0'(z)
%                VM --- Highest order computed
%       Routines called:
%(1)GAMMA for computing the gamma function
%(2)MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       ===========================================================
vg=[];ga=[];gb=[];
pi=3.141592653589793d0;
rp2=.63661977236758d0;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z1=z;
z2=z.*z;
n=fix(v);
v0=v-n;
pv0=pi.*v0;
if(a0 < 1.0d-100);
for  k=0:n;
cbj(k+1)=complex(0.0d0,0.0d0);
cdj(k+1)=complex(0.0d0,0.0d0);
cby(k+1)=-complex(1.0d+300,0.0d0);
cdy(k+1)=complex(1.0d+300,0.0d0);
end;  k=n+1;
if(v0 == 0.0);
cbj(0+1)=complex(1.0d0,0.0d0);
cdj(1+1)=complex(0.5d0,0.0d0);
else;
cdj(0+1)=complex(1.0d+300,0.0d0);
end;
vm=v;
return;
end;
if(real(z)< 0.0d0)z1=-z; end;
if(a0 <= 12.0);
cjv0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*(k+v0));
cjv0=cjv0+cr;
if(abs(cr)< abs(cjv0).*1.0d-15)break; end;
end;
vg=1.0d0+v0;
[vg,ga]=gamma(vg,ga);
ca=(0.5d0.*z1).^v0./ga;
cjv0=cjv0.*ca;
else;
k0=11;
if(a0 >= 35.0)k0=10; end;
if(a0 >= 50.0)k0=8; end;
vv=4.0d0.*v0.*v0;
cpz=complex(1.0d0,0.0d0);
crp=complex(1.0d0,0.0d0);
for  k=1:k0;
crp=-0.78125d-2.*crp.*(vv-(4.0.*k-3.0).^2.0).*(vv-(4.0.*k-1.0).^2.0)./(k.*(2.0.*k-1.0).*z2);
cpz=cpz+crp;
end;  k=k0+1;
cqz=complex(1.0d0,0.0d0);
crq=complex(1.0d0,0.0d0);
for  k=1:k0;
crq=-0.78125d-2.*crq.*(vv-(4.0.*k-1.0).^2.0).*(vv-(4.0.*k+1.0).^2.0)./(k.*(2.0.*k+1.0).*z2);
cqz=cqz+crq;
end;  k=k0+1;
cqz=0.125d0.*(vv-1.0).*cqz./z1;
zk=z1-(0.5d0.*v0+0.25d0).*pi;
ca0=sqrt(rp2./z1);
cck=cos(zk);
csk=sin(zk);
cjv0=ca0.*(cpz.*cck-cqz.*csk);
cyv0=ca0.*(cpz.*csk+cqz.*cck);
end;
if(a0 <= 12.0);
if(v0 ~= 0.0);
cjvn=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*(k-v0));
cjvn=cjvn+cr;
if(abs(cr)< abs(cjvn).*1.0d-15)break; end;
end;
vg=1.0d0-v0;
[vg,gb]=gamma(vg,gb);
cb=(2.0d0./z1).^v0./gb;
cju0=cjvn.*cb;
cyv0=(cjv0.*cos(pv0)-cju0)./sin(pv0);
else;
cec=log(z1./2.0d0)+.5772156649015329d0;
cs0=complex(0.0d0,0.0d0);
w0=0.0d0;
cr0=complex(1.0d0,0.0d0);
for  k=1:30;
w0=w0+1.0d0./k;
cr0=-0.25d0.*cr0./(k.*k).*z2;
cs0=cs0+cr0.*w0;
end;  k=30+1;
cyv0=rp2.*(cec.*cjv0-cs0);
end;
end;
if(n == 0)n=1; end;
m=msta1(a0,200);
if(m < n);
n=m;
else;
m=msta2(a0,n,15);
end;
cf2=complex(0.0d0,0.0d0);
cf1=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=2.0d0.*(v0+k+1.0d0)./z1.*cf1-cf2;
if(k <= n)cbj(k+1)=cf; end;
cf2=cf1;
cf1=cf;
end;  k=0-1;
cs=cjv0./cf;
for  k=0:n;
cbj(k+1)=cs.*cbj(k+1);
end;  k=n+1;
if(real(z)< 0.0d0);
cfac0=exp(pv0.*ci);
if(imag(z)< 0.0d0);
cyv0=cfac0.*cyv0-2.0d0.*ci.*cos(pv0).*cjv0;
elseif(imag(z)> 0.0d0);
cyv0=cyv0./cfac0+2.0d0.*ci.*cos(pv0).*cjv0;
end;
for  k=0:n;
if(imag(z)< 0.0d0);
cbj(k+1)=exp(-pi.*(k+v0).*ci).*cbj(k+1);
elseif(imag(z)> 0.0d0);
cbj(k+1)=exp(pi.*(k+v0).*ci).*cbj(k+1);
end;
end;  k=n+1;
z1=z1;
end;
cby(0+1)=cyv0;
for  k=1:n;
cyy=(cbj(k+1).*cby(k-1+1)-2.0d0./(pi.*z))./cbj(k-1+1);
cby(k+1)=cyy;
end;  k=n+1;
cdj(0+1)=v0./z.*cbj(0+1)-cbj(1+1);
for  k=1:n;
cdj(k+1)=-(k+v0)./z.*cbj(k+1)+cbj(k-1+1);
end;  k=n+1;
cdy(0+1)=v0./z.*cby(0+1)-cby(1+1);
for  k=1:n;
cdy(k+1)=cby(k-1+1)-(k+v0)./z.*cby(k+1);
end;  k=n+1;
vm=n+v0;
return;
end
function [x,ga]=gamma(x,ga,varargin);
%       ==================================================
%       Purpose: Compute gamma function ג(x)
%       Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%       Output:  GA --- ג(x)
%       ==================================================
 g=zeros(1,26);
pi=3.141592653589793d0;
if(x == fix(x));
if(x > 0.0d0);
ga=1.0d0;
m1=x-1;
for  k=2:m1;
ga=ga.*k;
end;  k=m1+1;
else;
ga=1.0d+300;
end;
else;
if(abs(x)> 1.0d0);
z=abs(x);
m=fix(z);
r=1.0d0;
for  k=1:m;
r=r.*(z-k);
end;  k=m+1;
z=z-m;
else;
z=x;
end;
g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
gr=g(26);
for  k=25:-1:1;
gr=gr.*z+g(k);
end;  k=1-1;
ga=1.0d0./(gr.*z);
if(abs(x)> 1.0d0);
ga=ga.*r;
if(x < 0.0d0)ga=-pi./(x.*ga.*sin(pi.*x)); end;
end;
end;
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

