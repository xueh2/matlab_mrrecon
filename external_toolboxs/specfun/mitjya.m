function mitjya
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program evaluates the integral of Bessel
%                functions J0(t)and Y0(t)with respect to t
%                from 0 to x using subroutine ITJYA
%       Input :  x  --- Upper limit of the integral(x ò 0)
%       Output:  TJ --- Integration of J0(t)from 0 to x
%                TY --- Integration of Y0(t)from 0 to x
%       Example:
%                   x         J0(t)dt          Y0(t)dt
%                ---------------------------------------
%                  5.0       .71531192       .19971938
%                 10.0      1.06701130       .24129032
%                 15.0      1.20516194       .00745772
%                 20.0      1.05837882      -.16821598
%                 25.0       .87101492      -.09360793
%                 30.0       .88424909       .08822971
%       ===========================================================
x=[];tj=[];ty=[];
  x=0;
 tj=0;
 ty=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=30.0;
fprintf(1,'%s \n','   x         j0(t)dt          y0(t)dt');
fprintf(1,'%s \n','---------------------------------------');
[x,tj,ty]=itjya(x,tj,ty);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,2) ' \n'],x,tj,ty);
%format(1x,f5.1,2f16.8);
end
function [x,tj,ty]=itjya(x,tj,ty,varargin);
%       ==========================================================
%       Purpose: Integrate Bessel functions J0(t)& Y0(t)with
%                respect to t from 0 to x
%       Input :  x  --- Upper limit of the integral(x ò 0)
%       Output:  TJ --- Integration of J0(t)from 0 to x
%                TY --- Integration of Y0(t)from 0 to x
%       =======================================================
 a=zeros(1,18);
pi=3.141592653589793d0;
el=.5772156649015329d0;
eps=1.0d-12;
if(x == 0.0d0);
tj=0.0d0;
ty=0.0d0;
elseif(x <= 20.0d0);
x2=x.*x;
tj=x;
r=x;
for  k=1:60;
r=-.25d0.*r.*(2.*k-1.0d0)./(2.*k+1.0d0)./(k.*k).*x2;
tj=tj+r;
if(abs(r)< abs(tj).*eps)break; end;
end;
ty1=(el+log(x./2.0d0)).*tj;
rs=0.0d0;
ty2=1.0d0;
r=1.0d0;
for  k=1:60;
r=-.25d0.*r.*(2.*k-1.0d0)./(2.*k+1.0d0)./(k.*k).*x2;
rs=rs+1.0d0./k;
r2=r.*(rs+1.0d0./(2.0d0.*k+1.0d0));
ty2=ty2+r2;
if(abs(r2)< abs(ty2).*eps)break; end;
end;
ty=(ty1-x.*ty2).*2.0d0./pi;
else;
a0=1.0d0;
a1=5.0d0./8.0d0;
a(1)=a1;
for  k=1:16;
af=((1.5d0.*(k+.5d0).*(k+5.0d0./6.0d0).*a1-.5d0.*(k+.5d0).*(k+.5d0).*(k-.5d0).*a0))./(k+1.0d0);
a(k+1)=af;
a0=a1;
a1=af;
end;  k=16+1;
bf=1.0d0;
r=1.0d0;
for  k=1:8;
r=-r./(x.*x);
bf=bf+a(2.*k).*r;
end;  k=8+1;
bg=a(1)./x;
r=1.0d0./x;
for  k=1:8;
r=-r./(x.*x);
bg=bg+a(2.*k+1).*r;
end;  k=8+1;
xp=x+.25d0.*pi;
rc=sqrt(2.0d0./(pi.*x));
tj=1.0d0-rc.*(bf.*cos(xp)+bg.*sin(xp));
ty=rc.*(bg.*cos(xp)-bf.*sin(xp));
end;
return;
end

