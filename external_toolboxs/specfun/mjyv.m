function mjyv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ============================================================
%     Purpose: This program computes Bessel functions Jv(x)and
%     Yv(x)and their derivatives using subroutine JYV
%     Input :  x --- Argument of Jv(x)and Yv(x)
%     v --- Order of Jv(x)and Yv(x)
%(v = n+v0,  0 ף n ף 250, 0 ף v0 < 1)
%     Output:  BJ(n)--- Jn+v0(x)
%     DJ(n)--- Jn+v0'(x)
%     BY(n)--- Yn+v0(x)
%     DY(n)--- Yn+v0'(x)
%     Example: Compute Jv(x)and Yv(x)and their derivatives
%     for v = 0.25(1.0)5.25 and x = 10.0
%     Computation results:
%     v =  5.25,      x = 10.00
%     v        Jv(x)Jv'(x)Yv(x)Yv'(x)
%     ------------------------------------------------------------
%     .25   -.20639379    -.13476340     .14493044    -.21381777
%     1.25    .12960355    -.22259423     .21744103     .11775031
%     2.25    .23879467     .07587475    -.09057018     .23781932
%     3.25   -.02214595     .24599211    -.25819761    -.00665596
%     4.25   -.25318954     .08545961    -.07725827    -.22536285
%     5.25   -.19306516    -.15183033     .19252809    -.17833551
%     ============================================================
v=[];x=[];vm=[];
 global bj;
global dj;
global by;
global dy;
fprintf(1,'%s \n','  please enter v, x ');
%     READ(*,*)V,X
v=5.25;
x=10.0;
fprintf(1,[repmat(' ',1,8),'v =','%6.2g',',    ','x =','%6.2g' ' \n'],v,x);
if(v <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%      READ(*,*)NS
ns=2;
end;
[v,x,vm,bj,dj,by,dy]=jyv(v,x,vm,bj,dj,by,dy);
nm=fix(vm);
v0=vm-nm;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   v         jv(x)jv''(x)');fprintf(1,'%s \n','          yv(x)yv''(x)');
fprintf(1,'%s ',' ---------------------------------------------');fprintf(1,'%s \n','------------------------');
for  k=0:ns:nm;
vk=k+v0;
fprintf(1,[repmat(' ',1,1),'%6.2g',repmat('%16.8g',1,4) ' \n'],vk,bj(k+1),dj(k+1),by(k+1),dy(k+1));
end;  k=nm+1;
%format(1x,f6.2,4d16.8);
%format(8x,,f6.2,',    ',,f6.2);
end
function [v,x,vm,bj,dj,by,dy]=jyv(v,x,vm,bj,dj,by,dy,varargin);
%     =======================================================
%     Purpose: Compute Bessel functions Jv(x)and Yv(x)
%     and their derivatives
%     Input :  x --- Argument of Jv(x)and Yv(x)
%     v --- Order of Jv(x)and Yv(x)
%(v = n+v0, 0 ף v0 < 1, n = 0,1,2,...)
%     Output:  BJ(n)--- Jn+v0(x)
%     DJ(n)--- Jn+v0'(x)
%     BY(n)--- Yn+v0(x)
%     DY(n)--- Yn+v0'(x)
%     VM --- Highest order computed
%     Routines called:
%(1)GAMMA for computing gamma function
%(2)MSTA1 and MSTA2 for computing the starting
%     point for backward recurrence
%     =======================================================
vg=[];ga=[];gb=[];
el=.5772156649015329d0;
pi=3.141592653589793d0;
rp2=.63661977236758d0;
x2=x.*x;
n=fix(v);
v0=v-n;
if(x < 1.0d-100);
for  k=0:n;
bj(k+1)=0.0d0;
dj(k+1)=0.0d0;
by(k+1)=-1.0d+300;
dy(k+1)=1.0d+300;
end;  k=n+1;
if(v0 == 0.0);
bj(0+1)=1.0d0;
dj(1+1)=0.5d0;
else;
dj(0+1)=1.0d+300;
end;
vm=v;
return;
end;
if(x <= 12.0);
for  l=0:1;
vl=v0+l;
bjvl=1.0d0;
r=1.0d0;
for  k=1:40;
r=-0.25d0.*r.*x2./(k.*(k+vl));
bjvl=bjvl+r;
if(abs(r)< abs(bjvl).*1.0d-15)break; end;
end;
vg=1.0d0+vl;
[vg,ga]=gamma(vg,ga);
a=(0.5d0.*x).^vl./ga;
if(l == 0)bjv0=bjvl.*a; end;
if(l == 1)bjv1=bjvl.*a; end;
end;
else;
k0=11;
if(x >= 35.0)k0=10; end;
if(x >= 50.0)k0=8; end;
for  j=0:1;
vv=4.0d0.*(j+v0).*(j+v0);
px=1.0d0;
rp=1.0d0;
for  k=1:k0;
rp=-0.78125d-2.*rp.*(vv-(4.0.*k-3.0).^2.0).*(vv-(4.0.*k-1.0).^2.0)./(k.*(2.0.*k-1.0).*x2);
px=px+rp;
end;  k=k0+1;
qx=1.0d0;
rq=1.0d0;
for  k=1:k0;
rq=-0.78125d-2.*rq.*(vv-(4.0.*k-1.0).^2.0).*(vv-(4.0.*k+1.0).^2.0)./(k.*(2.0.*k+1.0).*x2);
qx=qx+rq;
end;  k=k0+1;
qx=0.125d0.*(vv-1.0).*qx./x;
xk=x-(0.5d0.*(j+v0)+0.25d0).*pi;
a0=sqrt(rp2./x);
ck=cos(xk);
sk=sin(xk);
if(j == 0);
bjv0=a0.*(px.*ck-qx.*sk);
byv0=a0.*(px.*sk+qx.*ck);
elseif(j == 1);
bjv1=a0.*(px.*ck-qx.*sk);
byv1=a0.*(px.*sk+qx.*ck);
end;
end;  j=1+1;
end;
bj(0+1)=bjv0;
bj(1+1)=bjv1;
dj(0+1)=v0./x.*bj(0+1)-bj(1+1);
dj(1+1)=-(1.0d0+v0)./x.*bj(1+1)+bj(0+1);
if(n >= 2&n <= fix(0.9.*x));
f0=bjv0;
f1=bjv1;
for  k=2:n;
f=2.0d0.*(k+v0-1.0d0)./x.*f1-f0;
bj(k+1)=f;
f0=f1;
f1=f;
end;  k=n+1;
elseif(n >= 2);
m=msta1(x,200);
if(m < n);
n=m;
else;
m=msta2(x,n,15);
end;
f2=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(v0+k+1.0d0)./x.*f1-f2;
if(k <= n)bj(k+1)=f; end;
f2=f1;
f1=f;
end;  k=0-1;
if(abs(bjv0)> abs(bjv1));
cs=bjv0./f;
else;
cs=bjv1./f2;
end;
for  k=0:n;
bj(k+1)=cs.*bj(k+1);
end;  k=n+1;
end;
for  k=2:n;
dj(k+1)=-(k+v0)./x.*bj(k+1)+bj(k-1+1);
end;  k=n+1;
if(x <= 12.0d0);
if(v0 ~= 0.0);
for  l=0:1;
vl=v0+l;
bjvl=1.0d0;
r=1.0d0;
for  k=1:40;
r=-0.25d0.*r.*x2./(k.*(k-vl));
bjvl=bjvl+r;
if(abs(r)< abs(bjvl).*1.0d-15)break; end;
end;
vg=1.0d0-vl;
[vg,gb]=gamma(vg,gb);
b=(2.0d0./x).^vl./gb;
if(l == 0)bju0=bjvl.*b; end;
if(l == 1)bju1=bjvl.*b; end;
end;
pv0=pi.*v0;
pv1=pi.*(1.0d0+v0);
byv0=(bjv0.*cos(pv0)-bju0)./sin(pv0);
byv1=(bjv1.*cos(pv1)-bju1)./sin(pv1);
else;
ec=log(x./2.0d0)+el;
cs0=0.0d0;
w0=0.0d0;
r0=1.0d0;
for  k=1:30;
w0=w0+1.0d0./k;
r0=-0.25d0.*r0./(k.*k).*x2;
cs0=cs0+r0.*w0;
end;  k=30+1;
byv0=rp2.*(ec.*bjv0-cs0);
cs1=1.0d0;
w1=0.0d0;
r1=1.0d0;
for  k=1:30;
w1=w1+1.0d0./k;
r1=-0.25d0.*r1./(k.*(k+1)).*x2;
cs1=cs1+r1.*(2.0d0.*w1+1.0d0./(k+1.0d0));
end;  k=30+1;
byv1=rp2.*(ec.*bjv1-1.0d0./x-0.25d0.*x.*cs1);
end;
end;
by(0+1)=byv0;
by(1+1)=byv1;
for  k=2:n;
byvk=2.0d0.*(v0+k-1.0d0)./x.*byv1-byv0;
by(k+1)=byvk;
byv0=byv1;
byv1=byvk;
end;  k=n+1;
dy(0+1)=v0./x.*by(0+1)-by(1+1);
for  k=1:n;
dy(k+1)=-(k+v0)./x.*by(k+1)+by(k-1+1);
end;  k=n+1;
vm=n+v0;
return;
end
function [x,ga]=gamma(x,ga,varargin);
%     ==================================================
%     Purpose: Compute gamma function ג(x)
%     Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%     Output:  GA --- ג(x)
%     ==================================================
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

