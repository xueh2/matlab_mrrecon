function mstvl0
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =====================================================
%       Purpose: This program computes modified Struve
%                function L0(x)using subroutine STVL0
%       Input :  x   --- Argument of L0(x,x ò 0)
%       Output:  SL0 --- L0(x)
%       Example:
%                   x        L0(x)
%               ------------------------
%                  0.0   .00000000D+00
%                  5.0   .27105917D+02
%                 10.0   .28156522D+04
%                 15.0   .33964933D+06
%                 20.0   .43558283D+08
%                 30.0   .78167230D+12
%                 40.0   .14894775D+17
%                 50.0   .29325538D+21
%       =====================================================
x=[];sl0=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=50.0;
fprintf(1,'%s \n','   x        l0(x)');
fprintf(1,'%s \n','-----------------------');
[x,sl0]=stvl0(x,sl0);
fprintf(1,[repmat(' ',1,1),'%5.1g','%16.8g' ' \n'],x,sl0);
%format(1x,f5.1,d16.8);
end
function [x,sl0]=stvl0(x,sl0,varargin);
%       ================================================
%       Purpose: Compute modified Struve function L0(x)
%       Input :  x   --- Argument of L0(x,x ò 0)
%       Output:  SL0 --- L0(x)
%       ================================================
pi=3.141592653589793d0;
s=1.0d0;
r=1.0d0;
if(x <= 20.0d0);
a0=2.0d0.*x./pi;
for  k=1:60;
r=r.*(x./(2.0d0.*k+1.0d0)).^2;
s=s+r;
if(abs(r./s)< 1.0d-12)break; end;
end;
sl0=a0.*s;
else;
km=fix(.5.*(x+1.0));
if(x >= 50.0)km=25; end;
for  k=1:km;
r=r.*((2.0d0.*k-1.0d0)./x).^2;
s=s+r;
if(abs(r./s)< 1.0d-12)break; end;
end;
a1=exp(x)./sqrt(2.0d0.*pi.*x);
r=1.0d0;
bi0=1.0d0;
for  k=1:16;
r=0.125d0.*r.*(2.0d0.*k-1.0d0).^2./(k.*x);
bi0=bi0+r;
if(abs(r./bi0)< 1.0d-12)break; end;
end;
bi0=a1.*bi0;
sl0=-2.0d0./(pi.*x).*s+bi0;
end;
return;
end

