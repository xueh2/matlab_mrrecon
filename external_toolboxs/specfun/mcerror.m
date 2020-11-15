function mcerror
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes the error function erf(z)
%                for a complex argument using subroutine CERROR
%       Input :  x   --- Real part of z
%                y   --- Imaginary part of z(y ó 3.0)
%       Output:  ERR --- Real part of erf(z)
%                ERI --- Imaginary part of erf(z)
%       Example:
%                   x       y       Re[erf(z)]Im[erf(z)]
%                 ---------------------------------------------
%                  1.0     2.0      -.53664357     -5.04914370
%                  2.0     2.0      1.15131087       .12729163
%                  3.0     2.0       .99896328      -.00001155
%                  4.0     2.0      1.00000057      -.00000051
%                  5.0     2.0      1.00000000       .00000000
%       ============================================================
z=[];cer=[];
  x=0;
 y=0;
fprintf(1,'%s \n','x,y=?');
% READ(*,*)X,Y
x=2.0;
y=2.0;
fprintf(1,'%s \n','   x      y      re[erf(z)]im[erf(z)]');
fprintf(1,'%s \n',' ---------------------------------------------');
z=complex(x,y);
[z,cer]=cerror(z,cer);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,2),'%5.1g',repmat(' ',1,1),repmat('%16.8g',1,2) ' \n'],z,cer);
fprintf(1,'%0.15g ',z);fprintf(1,'%0.15g \n',cer);
%format(1x,f5.1,2x,f5.1,1x,2e16.8);
end
function [z,cer]=cerror(z,cer,varargin);
%       ====================================================
%       Purpose: Compute error function erf(z)for a complex
%                argument(z=x+iy)
%       Input :  z   --- Complex argument
%       Output:  CER --- erf(z)
%       ====================================================
  a0=0;
 pi=0;
a0=abs(z);
c0=exp(-z.*z);
pi=3.141592653589793d0;
z1=z;
if(real(z)< 0.0);
z1=-z;
end;
if(a0 <= 5.8d0);
cs=z1;
cr=z1;
for  k=1:120;
cr=cr.*z1.*z1./(k+0.5d0);
cs=cs+cr;
if(abs(cr./cs)< 1.0d-15)break; end;
end;
cer=2.0d0.*c0.*cs./sqrt(pi);
else;
cl=1.0d0./z1;
cr=cl;
for  k=1:13;
cr=-cr.*(k-0.5d0)./(z1.*z1);
cl=cl+cr;
if(abs(cr./cl)< 1.0d-15)break; end;
end;
cer=1.0d0-c0.*cl./sqrt(pi);
end;
if(real(z)< 0.0);
cer=-cer;
end;
return;
end

