function melit
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes complete and incomplete
%                elliptic integrals F(k,phi)and E(k,phi)using
%                subroutine ELIT
%       Input  : HK  --- Modulus k(0 ó k ó 1)
%                Phi --- Argument(in degrees)
%       Output : FE  --- F(k,phi)
%                EE  --- E(k,phi)
%       Example:
%                k = .5
%                 phi     F(k,phi)E(k,phi)
%                -----------------------------------
%                   0      .00000000      .00000000
%                  15      .26254249      .26106005
%                  30      .52942863      .51788193
%                  45      .80436610      .76719599
%                  60     1.08955067     1.00755556
%                  75     1.38457455     1.23988858
%                  90     1.68575035     1.46746221
%       ==========================================================
hk=[];phi=[];fe=[];ee=[];
fprintf(1,'%s \n','please enter k and phi(in degs.)');
%        READ(*,*)HK,PHI
hk=.5;
phi=90;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  phi     f(k,phi)e(k,phi)');
fprintf(1,'%s \n',' -----------------------------------');
[hk,phi,fe,ee]=elit(hk,phi,fe,ee);
fprintf(1,[repmat(' ',1,1),'%6.2g',repmat('%13.8g',1,2) ' \n'],phi,fe,ee);
%format(1x,f6.2,2f13.8);
end
function [hk,phi,fe,ee]=elit(hk,phi,fe,ee,varargin);
%       ==================================================
%       Purpose: Compute complete and incomplete elliptic
%                integrals F(k,phi)and E(k,phi)
%       Input  : HK  --- Modulus k(0 ó k ó 1)
%                Phi --- Argument(in degrees)
%       Output : FE  --- F(k,phi)
%                EE  --- E(k,phi)
%       ==================================================
g=0.0d0;
pi=3.14159265358979d0;
a0=1.0d0;
b0=sqrt(1.0d0-hk.*hk);
d0=(pi./180.0d0).*phi;
r=hk.*hk;
if(hk == 1.0d0&phi == 90.0d0);
fe=1.0d+300;
ee=1.0d0;
elseif(hk == 1.0d0);
fe=log((1.0d0+sin(d0))./cos(d0));
ee=sin(d0);
else;
fac=1.0d0;
for  n=1:40;
a=(a0+b0)./2.0d0;
b=sqrt(a0.*b0);
c=(a0-b0)./2.0d0;
fac=2.0d0.*fac;
r=r+fac.*c.*c;
if(phi ~= 90.0d0);
d=d0+atan((b0./a0).*tan(d0));
g=g+c.*sin(d);
d0=d+pi.*fix(d./pi+.5d0);
end;
a0=a;
b0=b;
if(c < 1.0d-7)break; end;
end;
ck=pi./(2.0d0.*a);
ce=pi.*(2.0d0-r)./(4.0d0.*a);
if(phi == 90.0d0);
fe=ck;
ee=ce;
else;
fe=d./(fac.*a);
ee=fe.*ce./ck+g;
end;
end;
return;
end

