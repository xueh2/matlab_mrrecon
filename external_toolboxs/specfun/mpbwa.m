function mpbwa
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes the parabolic cylinder
%                functions W(a,ñx)and their derivatives using
%                subroutine PBWA
%       Input  : a --- Parameter(0 ó |a| ó 5)
%                x --- Argument of W(a,ñx,0 ó |x| ó 5)
%       Output : W1F --- W(a,x)
%                W1D --- W'(a,x)
%                W2F --- W(a,-x)
%                W2D --- W'(a,-x)
%       Example: x = 5.0
%                 a      W(a,x)W'(a,x)W(a,-x)W'(a,-x)
%              ----------------------------------------------------
%                0.5   .1871153    .1915744  -.8556585   4.4682493
%                1.5  -.0215853    .0899870 -8.8586002  -9.3971967
%                0.0   .3009549   -.7148233   .6599634   1.7552224
%               -0.5  -.1934088  -1.3474400   .6448148   -.6781011
%               -1.5  -.5266539    .8219516  -.2822774  -1.4582283
%               -5.0   .0893618  -1.8118641   .5386084    .2698553
%       ============================================================
a=[];x=[];w1f=[];w1d=[];w2f=[];w2d=[];
fprintf(1,'%s \n','please enter a and x ');
%        READ(*,*)A,X
a=-5.0;
x=5.0;
fprintf(1,[repmat(' ',1,1),'a=','%5.1g',repmat(' ',1,3),'x=','%5.1g' ' \n'],a,x);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   a       w(a,x)w''(a,x)');fprintf(1,'%s \n','         w(a,-x)w''(a,-x)');
fprintf(1,'%s ',' -----------------------------------------');fprintf(1,'%s \n','----------------------------');
[a,x,w1f,w1d,w2f,w2d]=pbwa(a,x,w1f,w1d,w2f,w2d);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],a,w1f,w1d,w2f,w2d);
%format(1x,',f5.1,3x,',f5.1);
%format(1x,f5.1,4d16.8);
end
function [a,x,w1f,w1d,w2f,w2d]=pbwa(a,x,w1f,w1d,w2f,w2d,varargin);
%       ======================================================
%       Purpose: Compute parabolic cylinder functions W(a,ñx)
%                and their derivatives
%       Input  : a --- Parameter(0 ó |a| ó 5)
%                x --- Argument of W(a,ñx,0 ó |x| ó 5)
%       Output : W1F --- W(a,x)
%                W1D --- W'(a,x)
%                W2F --- W(a,-x)
%                W2D --- W'(a,-x)
%       Routine called:
%               CGAMA for computing complex gamma function
%       ======================================================
x1=[];y1=[];ugr=[];ugi=[];x2=[];vgr=[];vgi=[];
 h=zeros(1,100);
d=zeros(1,100);
x1=0.0;
eps=1.0d-15;
p0=0.59460355750136d0;
if(a == 0.0d0);
g1=3.625609908222d0;
g2=1.225416702465d0;
else;
x1=0.25d0;
y1=0.5d0.*a;
[x1,y1,dumvar3,ugr,ugi]=cgama(x1,y1,1,ugr,ugi);
g1=sqrt(ugr.*ugr+ugi.*ugi);
x2=0.75d0;
[x2,y1,dumvar3,vgr,vgi]=cgama(x2,y1,1,vgr,vgi);
g2=sqrt(vgr.*vgr+vgi.*vgi);
end;
f1=sqrt(g1./g2);
f2=sqrt(2.0d0.*g2./g1);
h0=1.0d0;
h1=a;
h(1)=a;
for  l1=4:2:200;
m=l1./2;
hl=a.*h1-0.25d0.*(l1-2.0d0).*(l1-3.0d0).*h0;
h(m)=hl;
h0=h1;
h1=hl;
end;  l1=200+1;
y1f=1.0d0;
r=1.0d0;
for  k=1:100;
r=0.5d0.*r.*x.*x./(k.*(2.0d0.*k-1.0d0));
r1=h(k).*r;
y1f=y1f+r1;
if(abs(r1./y1f)<= eps&k > 30)break; end;
end;
y1d=a;
r=1.0d0;
for  k=1:100;
r=0.5d0.*r.*x.*x./(k.*(2.0d0.*k+1.0d0));
r1=h(k+1).*r;
y1d=y1d+r1;
if(abs(r1./y1d)<= eps&k > 30)break; end;
end;
y1d=x.*y1d;
d1=1.0d0;
d2=a;
d(1)=1.0d0;
d(2)=a;
for  l2=5:2:160;
m=(l2+1)./2;
dl=a.*d2-0.25d0.*(l2-2.0d0).*(l2-3.0d0).*d1;
d(m)=dl;
d1=d2;
d2=dl;
end;  l2=160+1;
y2f=1.0d0;
r=1.0d0;
for  k=1:100;
r=0.5d0.*r.*x.*x./(k.*(2.0d0.*k+1.0d0));
r1=d(k+1).*r;
y2f=y2f+r1;
if(abs(r1./y2f)<= eps&k > 30)break; end;
end;
y2f=x.*y2f;
y2d=1.0d0;
r=1.0d0;
for  k=1:100;
r=0.5d0.*r.*x.*x./(k.*(2.0d0.*k-1.0d0));
r1=d(k+1).*r;
y2d=y2d+r1;
if(abs(r1./y2d)<= eps&k > 30)break; end;
end;
w1f=p0.*(f1.*y1f-f2.*y2f);
w2f=p0.*(f1.*y1f+f2.*y2f);
w1d=p0.*(f1.*y1d-f2.*y2d);
w2d=p0.*(f1.*y1d+f2.*y2d);
return;
end
function [x,y,kf,gr,gi]=cgama(x,y,kf,gr,gi,varargin);
%       =========================================================
%       Purpose: Compute complex gamma function â(z)or Ln[â(z)]
%       Input :  x  --- Real part of z
%                y  --- Imaginary part of z
%                KF --- Function code
%                       KF=0 for Ln[â(z)]
%                       KF=1 for â(z)
%       Output:  GR --- Real part of Ln[â(z)]or â(z)
%                GI --- Imaginary part of Ln[â(z)]or â(z)
%       ========================================================
 a=zeros(1,10);
x1=0.0;
pi=3.141592653589793d0;
a(:)=[8.333333333333333d-02,-2.777777777777778d-03,7.936507936507937d-04,-5.952380952380952d-04,8.417508417508418d-04,-1.917526917526918d-03,6.410256410256410d-03,-2.955065359477124d-02,1.796443723688307d-01,-1.39243221690590d+00];
if(y == 0.0d0&x == fix(x)&x <= 0.0d0);
gr=1.0d+300;
gi=0.0d0;
return;
elseif(x < 0.0d0);
x1=x;
y1=y;
x=-x;
y=-y;
end;
x0=x;
if(x <= 7.0);
na=fix(7-x);
x0=x+na;
end;
z1=sqrt(x0.*x0+y.*y);
th=atan(y./x0);
gr=(x0-.5d0).*log(z1)-th.*y-x0+0.5d0.*log(2.0d0.*pi);
gi=th.*(x0-0.5d0)+y.*log(z1)-y;
for  k=1:10;
t=z1.^(1-2.*k);
gr=gr+a(k).*t.*cos((2.0d0.*k-1.0d0).*th);
gi=gi-a(k).*t.*sin((2.0d0.*k-1.0d0).*th);
end;  k=10+1;
if(x <= 7.0);
gr1=0.0d0;
gi1=0.0d0;
for  j=0:na-1;
gr1=gr1+.5d0.*log((x+j).^2+y.*y);
gi1=gi1+atan(y./(x+j));
end;  j=na-1+1;
gr=gr-gr1;
gi=gi-gi1;
end;
if(x1 < 0.0d0);
z1=sqrt(x.*x+y.*y);
th1=atan(y./x);
sr=-sin(pi.*x).*cosh(pi.*y);
si=-cos(pi.*x).*sinh(pi.*y);
z2=sqrt(sr.*sr+si.*si);
th2=atan(si./sr);
if(sr < 0.0d0);
th2=pi+th2;
end;
gr=log(pi./(z1.*z2))-gr;
gi=-th1-th2-gi;
x=x1;
y=y1;
end;
if(kf == 1);
g0=exp(gr);
gr=g0.*cos(gi);
gi=g0.*sin(gi);
end;
return;
end

