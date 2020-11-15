function mittjya
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program computes the integral of[1-J0(t)]/t
%                with respect to t from 0 to x and Y0(t)/t with
%                respect to t from x to ì using subroutine ITTJYA
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTJ --- Integration of[1-J0(t)]/t from 0 to x
%                TTY --- Integration of Y0(t)/t from x to ì
%       Example:
%                  x[1-J0(t)]/tdt       Y0(t)/tdt
%               -------------------------------------------
%                 5.0     .15403472D+01    -.46322055D-01
%                10.0     .21778664D+01    -.22987934D-01
%                15.0     .25785507D+01     .38573574D-03
%                20.0     .28773106D+01     .85031527D-02
%                25.0     .31082313D+01     .35263393D-02
%       ===========================================================
x=[];ttj=[];tty=[];
  x=0;
 ttj=0;
 tty=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x[1-j0(t)]/tdt       y0(t)/tdt');
fprintf(1,'%s \n','-------------------------------------------');
[x,ttj,tty]=ittjya(x,ttj,tty);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%18.8g',1,2) ' \n'],x,ttj,tty);
%format(1x,f5.1,2d18.8);
end
function [x,ttj,tty]=ittjya(x,ttj,tty,varargin);
%       =========================================================
%       Purpose: Integrate[1-J0(t)]/t with respect to t from 0
%                to x, and Y0(t)/t with respect to t from x to ì
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTJ --- Integration of[1-J0(t)]/t from 0 to x
%                TTY --- Integration of Y0(t)/t from x to ì
%       =========================================================
pi=3.141592653589793d0;
el=.5772156649015329d0;
if(x == 0.0d0);
ttj=0.0d0;
tty=-1.0d+300;
elseif(x <= 20.0d0);
ttj=1.0d0;
r=1.0d0;
for  k=2:100;
r=-.25d0.*r.*(k-1.0d0)./(k.*k.*k).*x.*x;
ttj=ttj+r;
if(abs(r)< abs(ttj).*1.0d-12)break; end;
end;
ttj=ttj.*.125d0.*x.*x;
e0=.5d0.*(pi.*pi./6.0d0-el.*el)-(.5d0.*log(x./2.0d0)+el).*log(x./2.0d0);
b1=el+log(x./2.0d0)-1.5d0;
rs=1.0d0;
r=-1.0d0;
for  k=2:100;
r=-.25d0.*r.*(k-1.0d0)./(k.*k.*k).*x.*x;
rs=rs+1.0d0./k;
r2=r.*(rs+1.0d0./(2.0d0.*k)-(el+log(x./2.0d0)));
b1=b1+r2;
if(abs(r2)< abs(b1).*1.0d-12)break; end;
end;
tty=2.0d0./pi.*(e0+.125d0.*x.*x.*b1);
else;
a0=sqrt(2.0d0./(pi.*x));
for  l=0:1;
vt=4.0d0.*l.*l;
px=1.0d0;
r=1.0d0;
for  k=1:14;
r=-.0078125d0.*r.*(vt-(4.0d0.*k-3.0d0).^2)./(x.*k).*(vt-(4.0d0.*k-1.0d0).^2)./((2.0d0.*k-1.0d0).*x);
px=px+r;
if(abs(r)< abs(px).*1.0d-12)break; end;
end;
qx=1.0d0;
r=1.0d0;
for  k=1:14;
r=-.0078125d0.*r.*(vt-(4.0d0.*k-1.0d0).^2)./(x.*k).*(vt-(4.0d0.*k+1.0d0).^2)./(2.0d0.*k+1.0d0)./x;
qx=qx+r;
if(abs(r)< abs(qx).*1.0d-12)break; end;
end;
qx=.125d0.*(vt-1.0d0)./x.*qx;
xk=x-(.25d0+.5d0.*l).*pi;
bj1=a0.*(px.*cos(xk)-qx.*sin(xk));
by1=a0.*(px.*sin(xk)+qx.*cos(xk));
if(l == 0);
bj0=bj1;
by0=by1;
end;
end;
t=2.0d0./x;
g0=1.0d0;
r0=1.0d0;
for  k=1:10;
r0=-k.*k.*t.*t.*r0;
g0=g0+r0;
end;  k=10+1;
g1=1.0d0;
r1=1.0d0;
for  k=1:10;
r1=-k.*(k+1.0d0).*t.*t.*r1;
g1=g1+r1;
end;  k=10+1;
ttj=2.0d0.*g1.*bj0./(x.*x)-g0.*bj1./x+el+log(x./2.0d0);
tty=2.0d0.*g1.*by0./(x.*x)-g0.*by1./x;
end;
return;
end

