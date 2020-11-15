function mcikva
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==============================================================
%       Purpose: This program computes the modified Bessel functions
%                Iv(z), Kv(z)and their derivatives for an arbitrary
%                order and complex argument using subroutine CIKVA
%       Input :  z --- Complex argument z
%                v --- Real order of Iv(z)and Kv(z)
%(v =n+v0,  0 ף n ף 250, 0 ף v0 < 1)
%       Output:  CBI(n)--- In+v0(z)
%                CDI(n)--- In+v0'(z)
%                CBK(n)--- Kn+v0(z)
%                CDK(n)--- Kn+v0'(z)
%       Example: Compute Iv(z), Kv(z)and their derivatives for
%                v =n+v0, v0=0.25, n =0(1)5, and z =4.0 +i 2.0
%                Computation results:
%                v= n+v0,   v0 = .25,   z =  4.0+ i  2.0
%      n     Re[Iv(z)]Im[Iv(z)]Re[Iv'(z)]Im[Iv'(z)]
%    -----------------------------------------------------------------
%      0  -.19336550D+01  .10328998D+02 -.23119621D+01  .91612230D+01
%      1  -.24735044D+01  .85964317D+01 -.23898329D+01  .78707023D+01
%      2  -.28460107D+01  .54124063D+01 -.24105909D+01  .55204965D+01
%      3  -.23476775D+01  .24445612D+01 -.21145027D+01  .30604463D+01
%      4  -.13829947D+01  .70848630D+00 -.14732387D+01  .12545751D+01
%      5  -.59879982D+00  .64588999D-01 -.78816416D+00  .32629794D+00
%      n     Re[Kv(z)]Im[Kv(z)]Re[Kv'(z)]Im[Kv'(z)]
%     ----------------------------------------------------------------
%      0  -.64820386D-02 -.84715754D-02  .75118612D-02  .89920077D-02
%      1  -.80477525D-02 -.92535355D-02  .96506687D-02  .97789903D-02
%      2  -.12819299D-01 -.11086405D-01  .16310878D-01  .11358076D-01
%      3  -.24574004D-01 -.13462616D-01  .33167751D-01  .11850554D-01
%      4  -.53516204D-01 -.12614703D-01  .75424026D-01  .14407268D-02
%      5  -.12627405D+00  .10581162D-01  .18054884D+00 -.64789392D-01
%       ==============================================================
v=[];z=[];vm=[];
  global cbi;
global cdi;
global cbk;
global cdk;
fprintf(1,'%s \n','  please enter v, x,y(z=x+iy)');
%        READ(*,*)V,X,Y
v=1.25;
x=4.0;
y=2.0;
z=complex(x,y);
n=fix(v);
v0=v-n;
fprintf(1,[repmat(' ',1,8),'v= n+v0',',   ','v0 =','%7.2g',',   ','z =','%6.1g','+ i','%6.1g' ' \n'],v0,x,y);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[v,z,vm,cbi,cdi,cbk,cdk]=cikva(v,z,vm,cbi,cdi,cbk,cdk);
nm=fix(vm);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n      re[iv(z)]im[iv(z)]');fprintf(1,'%s \n','     re[iv''(z)]im[iv''(z)]');
fprintf(1,'%s ',' ---------------------------------------');fprintf(1,'%s \n','------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,cbi(k+1),cdi(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   n      re[kv(z)]im[kv(z)]');fprintf(1,'%s \n','     re[kv''(z)]im[kv''(z)]');
fprintf(1,'%s ',' ---------------------------------------');fprintf(1,'%s \n','------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%4g',repmat(' ',1,1),repmat('%16.8g',1,4) ' \n'],k,cbk(k+1),cdk(k+1));
end;  k=nm+1;
%format(1x,i4,1x,4d16.8);
%format(8x,'v= n+v0',',   ',',f7.2,',   ',',f6.1, '+ i',f6.1);
end
function [v,z,vm,cbi,cdi,cbk,cdk]=cikva(v,z,vm,cbi,cdi,cbk,cdk,varargin);
%       ============================================================
%       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
%                and their derivatives for an arbitrary order and
%                complex argument
%       Input :  z --- Complex argument
%                v --- Real order of Iv(z)and Kv(z)
%(v = n+v0, n = 0,1,2,תתת, 0 ף v0 < 1)
%       Output:  CBI(n)--- In+v0(z)
%                CDI(n)--- In+v0'(z)
%                CBK(n)--- Kn+v0(z)
%                CDK(n)--- Kn+v0'(z)
%                VM --- Highest order computed
%       Routines called:
%(1)GAMMA for computing the gamma function
%(2)MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       ============================================================
v0p=[];gap=[];v0n=[];gan=[];
pi=3.141592653589793d0;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z1=z;
z2=z.*z;
n=fix(v);
v0=v-n;
piv=pi.*v0;
vt=4.0d0.*v0.*v0;
if(n == 0)n=1; end;
if(a0 < 1.0d-100);
for  k=0:n;
cbi(k+1)=0.0d0;
cdi(k+1)=0.0d0;
cbk(k+1)=-1.0d+300;
cdk(k+1)=1.0d+300;
end;  k=n+1;
if(v0 == 0.0);
cbi(0+1)=complex(1.0d0,0.0d0);
cdi(1+1)=complex(0.5d0,0.0d0);
end;
vm=v;
return;
end;
k0=14;
if(a0 >= 35.0)k0=10; end;
if(a0 >= 50.0)k0=8; end;
if(real(z)< 0.0)z1=-z; end;
if(a0 < 18.0);
if(v0 == 0.0);
ca1=complex(1.0d0,0.0d0);
else;
v0p=1.0d0+v0;
[v0p,gap]=gamma(v0p,gap);
ca1=(0.5d0.*z1).^v0./gap;
end;
ci0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:50;
cr=0.25d0.*cr.*z2./(k.*(k+v0));
ci0=ci0+cr;
if(abs(cr)< abs(ci0).*1.0d-15)break; end;
end;
cbi0=ci0.*ca1;
else;
ca=exp(z1)./sqrt(2.0d0.*pi.*z1);
cs=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:k0;
cr=-0.125d0.*cr.*(vt-(2.0d0.*k-1.0d0).^2.0)./(k.*z1);
cs=cs+cr;
end;  k=k0+1;
cbi0=ca.*cs;
end;
m=msta1(a0,200);
if(m < n);
n=m;
else;
m=msta2(a0,n,15);
end;
cf2=complex(0.0d0,0.0d0);
cf1=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=2.0d0.*(v0+k+1.0d0)./z1.*cf1+cf2;
if(k <= n)cbi(k+1)=cf; end;
cf2=cf1;
cf1=cf;
end;  k=0-1;
cs=cbi0./cf;
for  k=0:n;
cbi(k+1)=cs.*cbi(k+1);
end;  k=n+1;
if(a0 <= 9.0);
if(v0 == 0.0);
ct=-log(0.5d0.*z1)-0.5772156649015329d0;
cs=complex(0.0d0,0.0d0);
w0=0.0d0;
cr=complex(1.0d0,0.0d0);
for  k=1:50;
w0=w0+1.0d0./k;
cr=0.25d0.*cr./(k.*k).*z2;
cp=cr.*(w0+ct);
cs=cs+cp;
if(k >= 10&abs(cp./cs)< 1.0d-15)break; end;
end;
cbk0=ct+cs;
else;
v0n=1.0d0-v0;
[v0n,gan]=gamma(v0n,gan);
ca2=1.0d0./(gan.*(0.5d0.*z1).^v0);
ca1=(0.5d0.*z1).^v0./gap;
csu=ca2-ca1;
cr1=complex(1.0d0,0.0d0);
cr2=complex(1.0d0,0.0d0);
for  k=1:50;
cr1=0.25d0.*cr1.*z2./(k.*(k-v0));
cr2=0.25d0.*cr2.*z2./(k.*(k+v0));
csu=csu+ca2.*cr1-ca1.*cr2;
ws=abs(csu);
if(k >= 10&abs(ws-ws0)./ws < 1.0d-15)break; end;
ws0=ws;
end;
cbk0=0.5d0.*pi.*csu./sin(piv);
end;
else;
cb=exp(-z1).*sqrt(0.5d0.*pi./z1);
cs=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:k0;
cr=0.125d0.*cr.*(vt-(2.0d0.*k-1.0d0).^2.0)./(k.*z1);
cs=cs+cr;
end;  k=k0+1;
cbk0=cb.*cs;
end;
cbk1=(1.0d0./z1-cbi(1+1).*cbk0)./cbi(0+1);
cbk(0+1)=cbk0;
cbk(1+1)=cbk1;
cg0=cbk0;
cg1=cbk1;
for  k=2:n;
cgk=2.0d0.*(v0+k-1.0d0)./z1.*cg1+cg0;
cbk(k+1)=cgk;
cg0=cg1;
cg1=cgk;
end;  k=n+1;
if(real(z)< 0.0);
for  k=0:n;
cvk=exp((k+v0).*pi.*ci);
if(imag(z)< 0.0d0);
cbk(k+1)=cvk.*cbk(k+1)+pi.*ci.*cbi(k+1);
cbi(k+1)=cbi(k+1)./cvk;
elseif(imag(z)> 0.0);
cbk(k+1)=cbk(k+1)./cvk-pi.*ci.*cbi(k+1);
cbi(k+1)=cvk.*cbi(k+1);
end;
end;  k=n+1;
end;
cdi(0+1)=v0./z.*cbi(0+1)+cbi(1+1);
cdk(0+1)=v0./z.*cbk(0+1)-cbk(1+1);
for  k=1:n;
cdi(k+1)=-(k+v0)./z.*cbi(k+1)+cbi(k-1+1);
cdk(k+1)=-(k+v0)./z.*cbk(k+1)-cbk(k-1+1);
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

