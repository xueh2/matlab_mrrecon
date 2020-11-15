function mlpmv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes the associated Legendre
%                function Pmv(x)with an integer order and an
%                arbitrary nonnegative degree using subroutine
%                LPMV
%       Input :  x   --- Argument of Pm(x,-1 ó x ó 1)
%                m   --- Order of Pmv(x)
%                v   --- Degree of Pmv(x)
%       Output:  PMV --- Pmv(x)
%       Example:    m = 4,  x = 0.5
%                    v          Pmv(x)
%                 -----------------------
%                   1.5       .46218726
%                   1.6       .48103143
%                   1.7       .45031429
%                   1.8       .36216902
%                   1.9       .21206446
%                   2.0       .00000000
%                   2.5     -1.51996235
%       =========================================================
v=[];m=[];x=[];pmv=[];
fprintf(1,'%s \n','please enter m,v,x = ?');
%        READ(*,*)M,V,X
m=4;
v=2.5;
x=0.5;
fprintf(1,[repmat(' ',1,3),'m =','%2g',',    ','x =','%6.2g' ' \n'],m,x);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','     v        pmv(x)');
fprintf(1,'%s \n','  -----------------------');
[v,m,x,pmv]=lpmv(v,m,x,pmv);
fprintf(1,[repmat(' ',1,3),'%5.1g','%16.8g' ' \n'],v,pmv);
%format(3x,f5.1,e16.8);
%format(3x,',i2,',    ',',f6.2);
end
function [v,m,x,pmv]=lpmv(v,m,x,pmv,varargin);
%       =======================================================
%       Purpose: Compute the associated Legendre function
%                Pmv(x)with an integer order and an arbitrary
%                nonnegative degree v
%       Input :  x   --- Argument of Pm(x,-1 ó x ó 1)
%                m   --- Order of Pmv(x)
%                v   --- Degree of Pmv(x)
%       Output:  PMV --- Pmv(x)
%       Routine called:  PSI for computing Psi function
%       =======================================================
psv=[];
pi=3.141592653589793d0;
el=.5772156649015329d0;
eps=1.0d-14;
nv=fix(v);
v0=v-nv;
if(x == -1.0d0&v ~= nv);
if(m == 0)pmv=-1.0d+300; end;
if(m ~= 0)pmv=1.0d+300; end;
return;
end;
c0=1.0d0;
if(m ~= 0);
rg=v.*(v+fix(m));
for  j=1:m-1;
rg=rg.*(v.*v-j.*j);
end;  j=fix(m)-1+1;
xq=sqrt(1.0d0-x.*x);
r0=1.0d0;
for  j=1:m;
r0=.5d0.*r0.*xq./j;
end;  j=fix(m)+1;
c0=r0.*rg;
end;
if(v0 == 0.0d0);
pmv=1.0d0;
r=1.0d0;
for  k=1:nv-m;
r=0.5d0.*r.*(-nv+fix(m)+k-1.0d0).*(nv+fix(m)+k)./(k.*(k+fix(m))).*(1.0d0+x);
pmv=pmv+r;
end;  k=nv-fix(m)+1;
pmv=(-1).^nv.*c0.*pmv;
else;
if(x >= -0.35d0);
pmv=1.0d0;
r=1.0d0;
for  k=1:100;
r=0.5d0.*r.*(-v+fix(m)+k-1.0d0).*(v+fix(m)+k)./(k.*(fix(m)+k)).*(1.0d0-x);
pmv=pmv+r;
if(k > 12&abs(r./pmv)< eps)break; end;
end;
pmv=(-1).^fix(m).*c0.*pmv;
else;
vs=sin(v.*pi)./pi;
pv0=0.0d0;
if(m ~= 0);
qr=sqrt((1.0d0-x)./(1.0d0+x));
r2=1.0d0;
for  j=1:m;
r2=r2.*qr.*j;
end;  j=fix(m)+1;
s0=1.0d0;
r1=1.0d0;
for  k=1:m-1;
r1=0.5d0.*r1.*(-v+k-1).*(v+k)./(k.*(k-fix(m))).*(1.0d0+x);
s0=s0+r1;
end;  k=fix(m)-1+1;
pv0=-vs.*r2./fix(m).*s0;
end;
[v,psv]=psi(v,psv);
pa=2.0d0.*(psv+el)+pi./tan(pi.*v)+1.0d0./v;
s1=0.0d0;
for  j=1:m;
s1=s1+(j.*j+v.*v)./(j.*(j.*j-v.*v));
end;  j=fix(m)+1;
pmv=pa+s1-1.0d0./(fix(m)-v)+log(0.5d0.*(1.0d0+x));
r=1.0d0;
for  k=1:100;
r=0.5d0.*r.*(-v+fix(m)+k-1.0d0).*(v+fix(m)+k)./(k.*(k+fix(m))).*(1.0d0+x);
s=0.0d0;
for  j=1:m;
s=s+((k+j).^2+v.*v)./((k+j).*((k+j).^2-v.*v));
end;  j=fix(m)+1;
s2=0.0d0;
for  j=1:k;
s2=s2+1.0d0./(j.*(j.*j-v.*v));
end;  j=k+1;
pss=pa+s+2.0d0.*v.*v.*s2-1.0d0./(fix(m)+k-v)+log(0.5d0.*(1.0d0+x));
r2=pss.*r;
pmv=pmv+r2;
if(abs(r2./pmv)< eps)break; end;
end;
pmv=pv0+pmv.*vs.*c0;
end;
end;
return;
end
function [x,ps]=psi(x,ps,varargin);
%       ======================================
%       Purpose: Compute psi function
%       Input :  x  --- Argument of psi(x)
%       Output:  PS --- psi(x)
%       ======================================
xa=abs(x);
pi=3.141592653589793d0;
el=.5772156649015329d0;
s=0.0d0;
if(x == fix(x)&x <= 0.0);
ps=1.0d+300;
return;
elseif(xa == fix(xa));
n=xa;
for  k=1 :n-1;
s=s+1.0d0./k;
end;  k=n-1+1;
ps=-el+s;
elseif(xa+0.5 == fix(xa+0.5));
n=xa-.5;
for  k=1:n;
s=s+1.0./(2.0d0.*k-1.0d0);
end;  k=n+1;
ps=-el+2.0d0.*s-1.386294361119891d0;
else;
if(xa < 10.0);
n=10-fix(xa);
for  k=0:n-1;
s=s+1.0d0./(xa+k);
end;  k=n-1+1;
xa=xa+n;
end;
x2=1.0d0./(xa.*xa);
a1=-.8333333333333d-01;
a2=.83333333333333333d-02;
a3=-.39682539682539683d-02;
a4=.41666666666666667d-02;
a5=-.75757575757575758d-02;
a6=.21092796092796093d-01;
a7=-.83333333333333333d-01;
a8=.4432598039215686d0;
ps=log(xa)-.5d0./xa+x2.*(((((((a8.*x2+a7).*x2+a6).*x2+a5).*x2+a4).*x2+a3).*x2+a2).*x2+a1);
ps=ps-s;
end;
if(x < 0.0)ps=ps-pi.*cos(pi.*x)./sin(pi.*x)-1.0d0./x; end;
return;
end

