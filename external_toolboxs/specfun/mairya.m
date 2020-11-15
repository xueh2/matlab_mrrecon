function mairya
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes Airy functions and their
%                derivatives using subroutine AIRYA
%       Input:   x  --- Argument of Airy function
%       Output:  AI --- Ai(x)
%                BI --- Bi(x)
%                AD --- Ai'(x)
%                BD --- Bi'(x)
%       Example:
%   x       Ai(x)Bi(x)Ai'(x)Bi'(x)
%  ----------------------------------------------------------------
%   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
%  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
%  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
%  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48
%   x       Ai(-x)Bi(-x)Ai'(-x)Bi'(-x)
%  ----------------------------------------------------------------
%   0       .35502805      .61492663     -.25881940      .44828836
%  10       .04024124     -.31467983      .99626504      .11941411
%  20      -.17640613     -.20013931      .89286286     -.79142903
%  30      -.08796819     -.22444694     1.22862060     -.48369473
%       ============================================================
x=[];ai=[];bi=[];ad=[];bd=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=10;
[x,ai,bi,ad,bd]=airya(x,ai,bi,ad,bd);
fprintf(1,[repmat(' ',1,4),'x',repmat(' ',1,8),'ai(x)',repmat(' ',1,11),'bi(x)',repmat(' ',1,11),'ai''(x)',repmat(' ',1,10),'bi''(x)' ' \n']);
fprintf(1,[repmat(' ',1,2),'----------------------------------','-----------------------------------' ' \n']);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x,ai,bi,ad,bd);
fprintf(1,'%0.15g \n');
[dumvar1,ai,bi,ad,bd]=airya(-x,ai,bi,ad,bd);
fprintf(1,[repmat(' ',1,4),'x',repmat(' ',1,8),'ai(-x)',repmat(' ',1,10),'bi(-x)',repmat(' ',1,10),'ai''(-x)',repmat(' ',1,9),'bi''(-x)' ' \n']);
fprintf(1,[repmat(' ',1,2),'----------------------------------','-----------------------------------' ' \n']);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x,ai,bi,ad,bd);
%format(1x,f5.1,4d16.8);
%format(1x,f5.1,4d16.8);
%format(4x,'x',8x,'ai(x)',11x,'bi(x)',11x,'ai''(x)', 10x,'bi''(x)');
%format(2x,'----------------------------------','-----------------------------------');
%format(4x,'x',8x,'ai(-x)',10x,'bi(-x)',10x, 'ai''(-x)',9x,'bi''(-x)');
end
function [x,ai,bi,ad,bd]=airya(x,ai,bi,ad,bd,varargin);
%       ======================================================
%       Purpose: Compute Airy functions and their derivatives
%       Input:   x  --- Argument of Airy function
%       Output:  AI --- Ai(x)
%                BI --- Bi(x)
%                AD --- Ai'(x)
%                BD --- Bi'(x)
%       Routine called:
%                AJYIK for computing Jv(x), Yv(x), Iv(x)and
%                Kv(x)with v=1/3 and 2/3
%       ======================================================
z=[];vj1=[];vj2=[];vy1=[];vy2=[];vi1=[];vi2=[];vk1=[];vk2=[];
xa=abs(x);
pir=0.318309886183891d0;
c1=0.355028053887817d0;
c2=0.258819403792807d0;
sr3=1.732050807568877d0;
z=xa.^1.5./1.5d0;
xq=sqrt(xa);
[z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2]=ajyik(z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2);
if(x == 0.0d0);
ai=c1;
bi=sr3.*c1;
ad=-c2;
bd=sr3.*c2;
elseif(x > 0.0d0);
ai=pir.*xq./sr3.*vk1;
bi=xq.*(pir.*vk1+2.0d0./sr3.*vi1);
ad=-xa./sr3.*pir.*vk2;
bd=xa.*(pir.*vk2+2.0d0./sr3.*vi2);
else;
ai=0.5d0.*xq.*(vj1-vy1./sr3);
bi=-0.5d0.*xq.*(vj1./sr3+vy1);
ad=0.5d0.*xa.*(vj2+vy2./sr3);
bd=0.5d0.*xa.*(vj2./sr3-vy2);
end;
return;
end
function [x,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2]=ajyik(x,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2,varargin);
%       =======================================================
%       Purpose: Compute Bessel functions Jv(x)and Yv(x),
%                and modified Bessel functions Iv(x)and
%                Kv(x), and their derivatives with v=1/3,2/3
%       Input :  x --- Argument of Jv(x),Yv(x),Iv(x)and
%                      Kv(x,x ò 0)
%       Output:  VJ1 --- J1/3(x)
%                VJ2 --- J2/3(x)
%                VY1 --- Y1/3(x)
%                VY2 --- Y2/3(x)
%                VI1 --- I1/3(x)
%                VI2 --- I2/3(x)
%                VK1 --- K1/3(x)
%                VK2 --- K2/3(x)
%       =======================================================
if(x == 0.0d0);
vj1=0.0d0;
vj2=0.0d0;
vy1=-1.0d+300;
vy2=1.0d+300;
vi1=0.0d0;
vi2=0.0d0;
vk1=-1.0d+300;
vk2=-1.0d+300;
return;
end;
pi=3.141592653589793d0;
rp2=.63661977236758d0;
gp1=.892979511569249d0;
gp2=.902745292950934d0;
gn1=1.3541179394264d0;
gn2=2.678938534707747d0;
vv0=0.444444444444444d0;
uu0=1.1547005383793d0;
x2=x.*x;
k0=12;
if(x >= 35.0)k0=10; end;
if(x >= 50.0)k0=8; end;
if(x <= 12.0);
for  l=1:2;
vl=l./3.0d0;
vjl=1.0d0;
r=1.0d0;
for  k=1:40;
r=-0.25d0.*r.*x2./(k.*(k+vl));
vjl=vjl+r;
if(abs(r)< 1.0d-15)break; end;
end;
a0=(0.5d0.*x).^vl;
if(l == 1)vj1=a0./gp1.*vjl; end;
if(l == 2)vj2=a0./gp2.*vjl; end;
end;
else;
for  l=1:2;
vv=vv0.*l.*l;
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
xk=x-(0.5d0.*l./3.0d0+0.25d0).*pi;
a0=sqrt(rp2./x);
ck=cos(xk);
sk=sin(xk);
if(l == 1);
vj1=a0.*(px.*ck-qx.*sk);
vy1=a0.*(px.*sk+qx.*ck);
elseif(l == 2);
vj2=a0.*(px.*ck-qx.*sk);
vy2=a0.*(px.*sk+qx.*ck);
end;
end;  l=2+1;
end;
if(x <= 12.0d0);
for  l=1:2;
vl=l./3.0d0;
vjl=1.0d0;
r=1.0d0;
for  k=1:40;
r=-0.25d0.*r.*x2./(k.*(k-vl));
vjl=vjl+r;
if(abs(r)< 1.0d-15)break; end;
end;
b0=(2.0d0./x).^vl;
if(l == 1)uj1=b0.*vjl./gn1; end;
if(l == 2)uj2=b0.*vjl./gn2; end;
end;
pv1=pi./3.0d0;
pv2=pi./1.5d0;
vy1=uu0.*(vj1.*cos(pv1)-uj1);
vy2=uu0.*(vj2.*cos(pv2)-uj2);
end;
if(x <= 18.0);
for  l=1:2;
vl=l./3.0d0;
vil=1.0d0;
r=1.0d0;
for  k=1:40;
r=0.25d0.*r.*x2./(k.*(k+vl));
vil=vil+r;
if(abs(r)< 1.0d-15)break; end;
end;
a0=(0.5d0.*x).^vl;
if(l == 1)vi1=a0./gp1.*vil; end;
if(l == 2)vi2=a0./gp2.*vil; end;
end;
else;
c0=exp(x)./sqrt(2.0d0.*pi.*x);
for  l=1:2;
vv=vv0.*l.*l;
vsl=1.0d0;
r=1.0d0;
for  k=1:k0;
r=-0.125d0.*r.*(vv-(2.0d0.*k-1.0d0).^2.0)./(k.*x);
vsl=vsl+r;
end;  k=k0+1;
if(l == 1)vi1=c0.*vsl; end;
if(l == 2)vi2=c0.*vsl; end;
end;  l=2+1;
end;
if(x <= 9.0d0);
for  l=1:2;
vl=l./3.0d0;
if(l == 1)gn=gn1; end;
if(l == 2)gn=gn2; end;
a0=(2.0d0./x).^vl./gn;
sum=1.0d0;
r=1.0d0;
for  k=1:60;
r=0.25d0.*r.*x2./(k.*(k-vl));
sum=sum+r;
if(abs(r)< 1.0d-15)break; end;
end;
if(l == 1)vk1=0.5d0.*uu0.*pi.*(sum.*a0-vi1); end;
if(l == 2)vk2=0.5d0.*uu0.*pi.*(sum.*a0-vi2); end;
end;
else;
c0=exp(-x).*sqrt(0.5d0.*pi./x);
for  l=1:2;
vv=vv0.*l.*l;
sum=1.0d0;
r=1.0d0;
for  k=1:k0;
r=0.125d0.*r.*(vv-(2.0.*k-1.0).^2.0)./(k.*x);
sum=sum+r;
end;  k=k0+1;
if(l == 1)vk1=c0.*sum; end;
if(l == 2)vk2=c0.*sum; end;
end;  l=2+1;
end;
return;
end

