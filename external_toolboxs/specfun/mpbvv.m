function mpbvv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes the parabolic cylinder
%                functions Vv(x)and Vv'(x)using subroutine
%                PBVV
%       Input:   x --- Argument of Vv(x)
%                v --- Order of Vv(x)
%       Output:  VV(na)--- Vv(x)
%                VP(na)--- Vv'(x)
%(na = |n|, v = n+v0, n = int(v), |v0| < 1
%                  n = 0,ס1,ס2,תתת, |n| ף 100)
%                PVF --- Vv(x)
%                PVD --- Vv'(x)
%       Example: v = 5.5,  x =10.0,  v0 = 0.5,  n = 0,1,2,...,5
%                  n+v0      Vv(x)Vv'(x)
%                ---------------------------------------
%                  0.5   .18522719D+10   .89761157D+10
%                  1.5   .19016268D+09   .90145854D+09
%                  2.5   .19741946D+08   .91452949D+08
%                  3.5   .20733667D+07   .93751130D+07
%                  4.5   .22038231D+06   .97145511D+06
%                  5.5   .23719356D+05   .10178553D+06
%                Vv(x)= .23719356D+05,  Vv'(x)= .10178553D+06
%       =========================================================
v=[];x=[];vv=[];vp=[];pvf=[];pvd=[];
 vv=zeros(1,100+1);
vp=zeros(1,100+1);
fprintf(1,'%s \n','please enter v and  x ');
%        READ(*,*)V,X
v=5.5;
x=10.0;
fprintf(1,[repmat(' ',1,1),'v =','%6.2g',',    ','x =','%6.2g' ' \n'],v,x);
nv=fix(v);
v0=v-nv;
na=abs(nv);
[v,x,vv,vp,pvf,pvd]=pbvv(v,x,vv,vp,pvf,pvd);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   v       vv(x)vv''(x)');
fprintf(1,'%s \n','---------------------------------------');
for  k=0:na;
vk=k.*(abs(1).*sign(nv))+v0;
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,2) ' \n'],vk,vv(k+1),vp(k+1));
end;  k=na+1;
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,1),'v =','%5.1g',',  vv(x)=','%14.8g',',   vv''(x)=','%14.8g' ' \n'],v,pvf,pvd);
%format(1x,',f6.2,',    ',',f6.2);
%format(1x,f5.1,2d16.8);
%format(1x,',f5.1,',d14.8,',d14.8);
end
function [v,x,vv,vp,pvf,pvd]=pbvv(v,x,vv,vp,pvf,pvd,varargin);
%       ===================================================
%       Purpose: Compute parabolic cylinder functions Vv(x)
%                and their derivatives
%       Input:   x --- Argument of Vv(x)
%                v --- Order of Vv(x)
%       Output:  VV(na)--- Vv(x)
%                VP(na)--- Vv'(x)
%(na = |n|, v = n+v0, |v0| < 1
%                  n = 0,ס1,ס2,תתת)
%                PVF --- Vv(x)
%                PVD --- Vv'(x)
%       Routines called:
%(1)VVSA for computing Vv(x)for small |x|
%(2)VVLA for computing Vv(x)for large |x|
%       ===================================================
v0=[];pv0=[];v1=[];f1=[];v2=[];f0=[];
pi=3.141592653589793d0;
xa=abs(x);
vh=v;
v=v+(abs(1.0d0).*sign(v));
nv=fix(v);
v0=v-nv;
na=abs(nv);
qe=exp(0.25d0.*x.*x);
q2p=sqrt(2.0d0./pi);
if(na >= 1)ja=1; end;
if(v <= 0.0);
if(v0 == 0.0);
if(xa <= 7.5)[v0,x,pv0]=vvsa(v0,x,pv0); end;
if(xa > 7.5)[v0,x,pv0]=vvla(v0,x,pv0); end;
f0=q2p.*qe;
f1=x.*f0;
vv(0+1)=pv0;
vv(1+1)=f0;
vv(2+1)=f1;
else;
for  l=0:ja;
v1=v0-l;
if(xa <= 7.5)[v1,x,f1]=vvsa(v1,x,f1); end;
if(xa > 7.5)[v1,x,f1]=vvla(v1,x,f1); end;
if(l == 0);
f0=f1;
end;
end;  l=ja+1;
vv(0+1)=f0;
vv(1+1)=f1;
end;
kv=2;
if(v0 == 0.0)kv=3; end;
for  k=kv:na;
f=x.*f1+(k-v0-2.0d0).*f0;
vv(k+1)=f;
f0=f1;
f1=f;
end;  k=na+1;
else;
if(x >= 0.0&x <= 7.5d0);
v2=v;
if(v2 < 1.0)v2=v2+1.0d0; end;
[v2,x,f1]=vvsa(v2,x,f1);
v1=v2-1.0d0;
kv=fix(v2);
[v1,x,f0]=vvsa(v1,x,f0);
vv(kv+1)=f1;
vv(kv-1+1)=f0;
for  k=kv-2:-1:0;
f=x.*f0-(k+v0+2.0d0).*f1;
if(k <= na)vv(k+1)=f; end;
f1=f0;
f0=f;
end;  k=0-1;
elseif(x > 7.5d0);
[v0,x,pv0]=vvla(v0,x,pv0);
m=100+abs(na);
vv(1+1)=pv0;
f1=0.0d0;
f0=1.0d-40;
for  k=m:-1:0;
f=x.*f0-(k+v0+2.0d0).*f1;
if(k <= na)vv(k+1)=f; end;
f1=f0;
f0=f;
end;  k=0-1;
s0=pv0./f;
for  k=0:na;
vv(k+1)=s0.*vv(k+1);
end;  k=na+1;
else;
if(xa <= 7.5d0);
[v0,x,f0]=vvsa(v0,x,f0);
v1=v0+1.0;
[v1,x,f1]=vvsa(v1,x,f1);
else;
[v0,x,f0]=vvla(v0,x,f0);
v1=v0+1.0d0;
[v1,x,f1]=vvla(v1,x,f1);
end;
vv(0+1)=f0;
vv(1+1)=f1;
for  k=2:na;
f=(x.*f1-f0)./(k+v0);
vv(k+1)=f;
f0=f1;
f1=f;
end;  k=na+1;
end;
end;
for  k=0:na-1;
v1=v0+k;
if(v >= 0.0d0);
vp(k+1)=0.5d0.*x.*vv(k+1)-(v1+1.0d0).*vv(k+1+1);
else;
vp(k+1)=-0.5d0.*x.*vv(k+1)+vv(k+1+1);
end;
end;  k=na-1+1;
pvf=vv(na-1+1);
pvd=vp(na-1+1);
v=vh;
return;
end
function [va,x,pv]=vvsa(va,x,pv,varargin);
%       ===================================================
%       Purpose: Compute parabolic cylinder function Vv(x)
%                for small argument
%       Input:   x  --- Argument
%                va --- Order
%       Output:  PV --- Vv(x)
%       Routine called : GAMMA for computing ג(x)
%       ===================================================
va0=[];ga0=[];v1=[];g1=[];vm=[];gm=[];
eps=1.0d-15;
pi=3.141592653589793d0;
ep=exp(-.25d0.*x.*x);
va0=1.0d0+0.5d0.*va;
if(x == 0.0);
if(va0 <= 0.0&va0 == fix(va0)|va == 0.0);
pv=0.0d0;
else;
vb0=-0.5d0.*va;
sv0=sin(va0.*pi);
[va0,ga0]=gamma(va0,ga0);
pv=2.0d0.^vb0.*sv0./ga0;
end;
else;
sq2=sqrt(2.0d0);
a0=2.0d0.^(-.5d0.*va).*ep./(2.0d0.*pi);
sv=sin(-(va+.5d0).*pi);
v1=-.5d0.*va;
[v1,g1]=gamma(v1,g1);
pv=(sv+1.0d0).*g1;
r=1.0d0;
fac=1.0d0;
for  m=1:250;
vm=.5d0.*(m-va);
[vm,gm]=gamma(vm,gm);
r=r.*sq2.*x./m;
fac=-fac;
gw=fac.*sv+1.0d0;
r1=gw.*r.*gm;
pv=pv+r1;
if(abs(r1./pv)< eps&gw ~= 0.0)break; end;
end;
pv=a0.*pv;
end;
return;
end
function [va,x,pv]=vvla(va,x,pv,varargin);
%       ===================================================
%       Purpose: Compute parabolic cylinder function Vv(x)
%                for large argument
%       Input:   x  --- Argument
%                va --- Order
%       Output:  PV --- Vv(x)
%       Routines called:
%(1)DVLA for computing Dv(x)for large |x|
%(2)GAMMA for computing ג(x)
%       ===================================================
x1=[];pdl=[];gl=[];
pi=3.141592653589793d0;
eps=1.0d-12;
qe=exp(0.25.*x.*x);
a0=abs(x).^(-va-1.0d0).*sqrt(2.0d0./pi).*qe;
r=1.0d0;
pv=1.0d0;
for  k=1:18;
r=0.5d0.*r.*(2.0.*k+va-1.0).*(2.0.*k+va)./(k.*x.*x);
pv=pv+r;
if(abs(r./pv)< eps)break; end;
end;
pv=a0.*pv;
if(x < 0.0d0);
x1=-x;
[va,x1,pdl]=dvla(va,x1,pdl);
[dumvar1,gl]=gamma(-va,gl);
dsl=sin(pi.*va).*sin(pi.*va);
pv=dsl.*gl./pi.*pdl-cos(pi.*va).*pv;
end;
return;
end
function [va,x,pd]=dvla(va,x,pd,varargin);
%       ====================================================
%       Purpose: Compute parabolic cylinder functions Dv(x)
%                for large argument
%       Input:   x  --- Argument
%                va --- Order
%       Output:  PD --- Dv(x)
%       Routines called:
%(1)VVLA for computing Vv(x)for large |x|
%(2)GAMMA for computing ג(x)
%       ====================================================
x1=[];vl=[];gl=[];
pi=3.141592653589793d0;
eps=1.0d-12;
ep=exp(-.25.*x.*x);
a0=x.^va.*ep;
r=1.0d0;
pd=1.0d0;
for  k=1:16;
r=-0.5d0.*r.*(2.0.*k-va-1.0).*(2.0.*k-va-2.0)./(k.*x.*x);
pd=pd+r;
if(abs(r./pd)< eps)break; end;
end;
pd=a0.*pd;
if(x < 0.0d0);
x1=-x;
[va,x1,vl]=vvla(va,x1,vl);
[dumvar1,gl]=gamma(-va,gl);
pd=pi.*vl./gl+cos(pi.*va).*pd;
end;
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

