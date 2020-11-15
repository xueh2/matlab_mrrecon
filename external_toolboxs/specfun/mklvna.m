function mklvna
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =======================================================
%       Purpose: This program computes Kelvin functions ber x,
%                bei x, ker x and kei x, and their derivatives
%                using subroutine KLVNA
%       Input :  x   --- Argument of Kelvin functions
%       Output:  BER --- ber x
%                BEI --- bei x
%                GER --- ker x
%                GEI --- kei x
%                DER --- ber'x
%                DEI --- bei'x
%                HER --- ker'x
%                HEI --- kei'x
%       Example:
%      x       ber x          bei x          ker x          kei x
%    -----------------------------------------------------------------
%      0    .1000000D+01    0              ì            -.7853982D+00
%      5   -.6230082D+01   .1160344D+00  -.1151173D-01   .1118759D-01
%     10    .1388405D+03   .5637046D+02   .1294663D-03  -.3075246D-03
%     15   -.2967255D+04  -.2952708D+04  -.1514347D-07   .7962894D-05
%     20    .4748937D+05   .1147752D+06  -.7715233D-07  -.1858942D-06
%      x       ber'x          bei'x          ker'x          kei'x
%    -----------------------------------------------------------------
%      0     0              0            - ì              0
%      5   -.3845339D+01  -.4354141D+01   .1719340D-01  -.8199865D-03
%     10    .5119526D+02   .1353093D+03  -.3155969D-03   .1409138D-03
%     15    .9105533D+02  -.4087755D+04   .5644678D-05  -.5882223D-05
%     20   -.4880320D+05   .1118550D+06  -.7501859D-07   .1906243D-06
%       =======================================================
x=[];ber=[];bei=[];ger=[];gei=[];der=[];dei=[];her=[];hei=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=20;
[x,ber,bei,ger,gei,der,dei,her,hei]=klvna(x,ber,bei,ger,gei,der,dei,her,hei);
fprintf(1,'%s ','   x        ber x           bei x');fprintf(1,'%s \n','           ker x           kei x');
fprintf(1,'%s ','--------------------------------');fprintf(1,'%s \n','--------------------------------------');
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x,ber,bei,ger,gei);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   x        ber''x           bei''x');fprintf(1,'%s \n','           ker''x           kei''x');
fprintf(1,'%s ','--------------------------------');fprintf(1,'%s \n','--------------------------------------');
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x,der,dei,her,hei);
%format(1x,f5.1,4d16.8);
end
function [x,ber,bei,ger,gei,der,dei,her,hei]=klvna(x,ber,bei,ger,gei,der,dei,her,hei,varargin);
%       ======================================================
%       Purpose: Compute Kelvin functions ber x, bei x, ker x
%                and kei x, and their derivatives(x > 0)
%       Input :  x   --- Argument of Kelvin functions
%       Output:  BER --- ber x
%                BEI --- bei x
%                GER --- ker x
%                GEI --- kei x
%                DER --- ber'x
%                DEI --- bei'x
%                HER --- ker'x
%                HEI --- kei'x
%       ================================================
pi=3.141592653589793d0;
el=.5772156649015329d0;
eps=1.0d-15;
if(x == 0.0d0);
ber=1.0d0;
bei=0.0d0;
ger=1.0d+300;
gei=-0.25d0.*pi;
der=0.0d0;
dei=0.0d0;
her=-1.0d+300;
hei=0.0d0;
return;
end;
x2=0.25d0.*x.*x;
x4=x2.*x2;
if(abs(x)< 10.0d0);
ber=1.0d0;
r=1.0d0;
for  m=1:60;
r=-0.25d0.*r./(m.*m)./(2.0d0.*m-1.0d0).^2.*x4;
ber=ber+r;
if(abs(r)< abs(ber).*eps)break; end;
end;
bei=x2;
r=x2;
for  m=1:60;
r=-0.25d0.*r./(m.*m)./(2.0d0.*m+1.0d0).^2.*x4;
bei=bei+r;
if(abs(r)< abs(bei).*eps)break; end;
end;
ger=-(log(x./2.0d0)+el).*ber+0.25d0.*pi.*bei;
r=1.0d0;
gs=0.0d0;
for  m=1:60;
r=-0.25d0.*r./(m.*m)./(2.0d0.*m-1.0d0).^2.*x4;
gs=gs+1.0d0./(2.0d0.*m-1.0d0)+1.0d0./(2.0d0.*m);
ger=ger+r.*gs;
if(abs(r.*gs)< abs(ger).*eps)break; end;
end;
gei=x2-(log(x./2.0d0)+el).*bei-0.25d0.*pi.*ber;
r=x2;
gs=1.0d0;
for  m=1:60;
r=-0.25d0.*r./(m.*m)./(2.0d0.*m+1.0d0).^2.*x4;
gs=gs+1.0d0./(2.0d0.*m)+1.0d0./(2.0d0.*m+1.0d0);
gei=gei+r.*gs;
if(abs(r.*gs)< abs(gei).*eps)break; end;
end;
der=-0.25d0.*x.*x2;
r=der;
for  m=1:60;
r=-0.25d0.*r./m./(m+1.0d0)./(2.0d0.*m+1.0d0).^2.*x4;
der=der+r;
if(abs(r)< abs(der).*eps)break; end;
end;
dei=0.5d0.*x;
r=dei;
for  m=1:60;
r=-0.25d0.*r./(m.*m)./(2.d0.*m-1.d0)./(2.d0.*m+1.d0).*x4;
dei=dei+r;
if(abs(r)< abs(dei).*eps)break; end;
end;
r=-0.25d0.*x.*x2;
gs=1.5d0;
her=1.5d0.*r-ber./x-(log(x./2.d0)+el).*der+0.25.*pi.*dei;
for  m=1:60;
r=-0.25d0.*r./m./(m+1.0d0)./(2.0d0.*m+1.0d0).^2.*x4;
gs=gs+1.0d0./(2.*m+1.0d0)+1.0d0./(2.*m+2.0d0);
her=her+r.*gs;
if(abs(r.*gs)< abs(her).*eps)break; end;
end;
r=0.5d0.*x;
gs=1.0d0;
hei=0.5d0.*x-bei./x-(log(x./2.d0)+el).*dei-0.25.*pi.*der;
for  m=1:60;
r=-0.25d0.*r./(m.*m)./(2.*m-1.0d0)./(2.*m+1.0d0).*x4;
gs=gs+1.0d0./(2.0d0.*m)+1.0d0./(2.*m+1.0d0);
hei=hei+r.*gs;
if(abs(r.*gs)< abs(hei).*eps)return; end;
end;  m=60+1;
else;
pp0=1.0d0;
pn0=1.0d0;
qp0=0.0d0;
qn0=0.0d0;
r0=1.0d0;
km=18;
if(abs(x)>= 40.0)km=10; end;
fac=1.0d0;
for  k=1:km;
fac=-fac;
xt=0.25d0.*k.*pi-fix(0.125d0.*k).*2.0d0.*pi;
cs=cos(xt);
ss=sin(xt);
r0=0.125d0.*r0.*(2.0d0.*k-1.0d0).^2./k./x;
rc=r0.*cs;
rs=r0.*ss;
pp0=pp0+rc;
pn0=pn0+fac.*rc;
qp0=qp0+rs;
qn0=qn0+fac.*rs;
end;  k=km+1;
xd=x./sqrt(2.0d0);
xe1=exp(xd);
xe2=exp(-xd);
xc1=1.d0./sqrt(2.0d0.*pi.*x);
xc2=sqrt(.5d0.*pi./x);
cp0=cos(xd+0.125d0.*pi);
cn0=cos(xd-0.125d0.*pi);
sp0=sin(xd+0.125d0.*pi);
sn0=sin(xd-0.125d0.*pi);
ger=xc2.*xe2.*(pn0.*cp0-qn0.*sp0);
gei=xc2.*xe2.*(-pn0.*sp0-qn0.*cp0);
ber=xc1.*xe1.*(pp0.*cn0+qp0.*sn0)-gei./pi;
bei=xc1.*xe1.*(pp0.*sn0-qp0.*cn0)+ger./pi;
pp1=1.0d0;
pn1=1.0d0;
qp1=0.0d0;
qn1=0.0d0;
r1=1.0d0;
fac=1.0d0;
for  k=1:km;
fac=-fac;
xt=0.25d0.*k.*pi-fix(0.125d0.*k).*2.0d0.*pi;
cs=cos(xt);
ss=sin(xt);
r1=0.125d0.*r1.*(4.d0-(2.0d0.*k-1.0d0).^2)./k./x;
rc=r1.*cs;
rs=r1.*ss;
pp1=pp1+fac.*rc;
pn1=pn1+rc;
qp1=qp1+fac.*rs;
qn1=qn1+rs;
end;  k=km+1;
her=xc2.*xe2.*(-pn1.*cn0+qn1.*sn0);
hei=xc2.*xe2.*(pn1.*sn0+qn1.*cn0);
der=xc1.*xe1.*(pp1.*cp0+qp1.*sp0)-hei./pi;
dei=xc1.*xe1.*(pp1.*sp0-qp1.*cp0)+her./pi;
end;
return;
end

