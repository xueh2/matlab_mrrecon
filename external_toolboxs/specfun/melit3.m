function melit3
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes the elliptic integral of
%                the third kind using subroutine ELIT3
%       Input :  Phi --- Argument(in degrees)
%                 k  --- Modulus(0 ó k ó 1)
%                 c  --- Parameter(0 ó c ó 1)
%       Output:  EL3 ÄÄÄ Value of the elliptic integral of the
%                        third kind
%       ==========================================================
phi=[];hk=[];c=[];el3=[];
fprintf(1,'%s \n','please enter phi, k and c ');
%        READ(*,*)PHI,HK,C
phi=45;
hk=.5;
c=.5;
[phi,hk,c,el3]=elit3(phi,hk,c,el3);
fprintf(1,[repmat(' ',1,1),'el3=','%12.8g' ' \n'],el3);
%format(1x,',f12.8);
end
function [phi,hk,c,el3]=elit3(phi,hk,c,el3,varargin);
%       =========================================================
%       Purpose: Compute the elliptic integral of the third kind
%                using Gauss-Legendre quadrature
%       Input :  Phi --- Argument(in degrees)
%                 k  --- Modulus(0 ó k ó 1.0)
%                 c  --- Parameter(0 ó c ó 1.0)
%       Output:  EL3 --- Value of the elliptic integral of the
%                        third kind
%       =========================================================
 t=zeros(1,10);
w=zeros(1,10);
  lb1=false;
 lb2=false;
t(:)=[.9931285991850949,.9639719272779138,.9122344282513259,.8391169718222188,.7463319064601508,.6360536807265150,.5108670019508271,.3737060887154195,.2277858511416451,.7652652113349734d-1];
w(:)=[.1761400713915212d-1,.4060142980038694d-1,.6267204833410907d-1,.8327674157670475d-1,.1019301198172404,.1181945319615184,.1316886384491766,.1420961093183820,.1491729864726037,.1527533871307258];
lb1=hk == 1.0d0&abs(phi-90.0)<= 1.0d-8;
lb2=c == 1.0d0&abs(phi-90.0)<= 1.0d-8;
if(lb1|lb2);
el3=1.0d+300;
return;
end;
c1=0.87266462599716d-2.*phi;
c2=c1;
el3=0.0d0;
for  i=1:10;
c0=c2.*t(i);
t1=c1+c0;
t2=c1-c0;
f1=1.0d0./((1.0d0-c.*sin(t1).*sin(t1)).* sqrt(1.0d0-hk.*hk.*sin(t1).*sin(t1)));
f2=1.0d0./((1.0d0-c.*sin(t2).*sin(t2)).* sqrt(1.0d0-hk.*hk.*sin(t2).*sin(t2)));
el3=el3+w(i).*(f1+f2);
end;  i=10+1;
el3=c1.*el3;
return;
end

