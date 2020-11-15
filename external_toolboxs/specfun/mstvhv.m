function mstvhv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ======================================================
%       Purpose:  This program computes Struve function Hv(x)
%                 for an arbitrary order using subroutine
%                 STVHV
%       Input :   v  --- Order of Hv(x,-8.0 ף v ף 12.5)
%                 x  --- Argument of Hv(x,x ע 0)
%       Output:   HV --- Hv(x)
%       Example:  x = 10.0
%                   v           Hv(x)
%                 -----------------------
%                   .5     .46402212D+00
%                  1.5     .14452322D+01
%                  2.5     .31234632D+01
%                  3.5     .53730255D+01
%                  4.5     .72083122D+01
%                  5.5     .76851132D+01
%       ======================================================
v=[];x=[];hv=[];
fprintf(1,'%s \n','please v and x ');
%        READ(*,*)V,X
v=5.5;
x=10.0;
fprintf(1,[repmat(' ',1,1),'v =','%5.1g',repmat(' ',1,6),'x =','%5.1g' ' \n'],v,x);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   v           hv(x)');
fprintf(1,'%s \n',' -----------------------');
[v,x,hv]=stvhv(v,x,hv);
fprintf(1,[repmat(' ',1,1),'%5.1g','%18.8g' ' \n'],v,hv);
%format(1x,f5.1,d18.8);
%format(1x,',f5.1,6x,',f5.1);
end
function [v,x,hv]=stvhv(v,x,hv,varargin);
%       =====================================================
%       Purpose: Compute Struve function Hv(x)with an
%                arbitrary order v
%       Input :  v  --- Order of Hv(x,-8.0 ף v ף 12.5)
%                x  --- Argument of Hv(x,x ע 0)
%       Output:  HV --- Hv(x)
%       Routine called: GAMMA to compute the gamma function
%       =====================================================
v0=[];ga=[];va=[];vb=[];gb=[];
pi=3.141592653589793d0;
if(x == 0.0d0);
if(v > -1.0|fix(v)-v == 0.5d0);
hv=0.0d0;
elseif(v < -1.0d0);
hv=(-1).^(fix(0.5d0-v)-1).*1.0d+300;
elseif(v == -1.0d0);
hv=2.0d0./pi;
end;
return;
end;
if(x <= 20.0d0);
v0=v+1.5d0;
[v0,ga]=gamma(v0,ga);
s=2.0d0./(sqrt(pi).*ga);
r1=1.0d0;
for  k=1:100;
va=k+1.5d0;
[va,ga]=gamma(va,ga);
vb=v+k+1.5d0;
[vb,gb]=gamma(vb,gb);
r1=-r1.*(0.5d0.*x).^2;
r2=r1./(ga.*gb);
s=s+r2;
if(abs(r2)< abs(s).*1.0d-12)break; end;
end;
hv=(0.5d0.*x).^(v+1.0d0).*s;
else;
sa=(0.5d0.*x).^(v-1.0)./pi;
v0=v+0.5d0;
[v0,ga]=gamma(v0,ga);
s=sqrt(pi)./ga;
r1=1.0d0;
for  k=1:12;
va=k+0.5d0;
[va,ga]=gamma(va,ga);
vb=-k+v+0.5d0;
[vb,gb]=gamma(vb,gb);
r1=r1./(0.5d0.*x).^2;
s=s+r1.*ga./gb;
end;  k=12+1;
s0=sa.*s;
u=abs(v);
n=fix(u);
u0=u-n;
for  l=0:1;
vt=4.0d0.*(u0+l).^2;
r1=1.0d0;
pu1=1.0d0;
for  k=1:12;
r1=-0.0078125d0.*r1.*(vt-(4.0.*k-3.0d0).^2).*(vt-(4.0d0.*k-1.0).^2)./((2.0d0.*k-1.0).*k.*x.*x);
pu1=pu1+r1;
end;  k=12+1;
qu1=1.0d0;
r2=1.0d0;
for  k=1:12;
r2=-0.0078125d0.*r2.*(vt-(4.0d0.*k-1.0).^2).*(vt-(4.0d0.*k+1.0).^2)./((2.0d0.*k+1.0).*k.*x.*x);
qu1=qu1+r2;
end;  k=12+1;
qu1=0.125d0.*(vt-1.0d0)./x.*qu1;
if(l == 0);
pu0=pu1;
qu0=qu1;
end;
end;  l=1+1;
t0=x-(0.5.*u0+0.25d0).*pi;
t1=x-(0.5.*u0+0.75d0).*pi;
sr=sqrt(2.0d0./(pi.*x));
by0=sr.*(pu0.*sin(t0)+qu0.*cos(t0));
by1=sr.*(pu1.*sin(t1)+qu1.*cos(t1));
bf0=by0;
bf1=by1;
for  k=2:n;
bf=2.0d0.*(k-1.0+u0)./x.*bf1-bf0;
bf0=bf1;
bf1=bf;
end;  k=n+1;
if(n == 0)byv=by0; end;
if(n == 1)byv=by1; end;
if(n > 1)byv=bf; end;
hv=byv+s0;
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

