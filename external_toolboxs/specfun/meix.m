function meix
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes the exponential integral
%                Ei(x)using subroutine EIX
%       Example:
%                  x        Ei(x)
%                -----------------------
%                  0    -.10000000+301
%                  1     .18951178E+01
%                  2     .49542344E+01
%                  3     .99338326E+01
%                  4     .19630874E+02
%                  5     .40185275E+02
%       =========================================================
x=[];ei=[];
  ei=0;
 x=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=5;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   x         ei(x)');
fprintf(1,'%s \n','------------------------');
[x,ei]=eix(x,ei);
fprintf(1,[repmat(' ',1,1),'%5.1g','%18.8g' ' \n'],x,ei);
%format(1x,f5.1,e18.8);
end
function [x,ei]=eix(x,ei,varargin);
%       ============================================
%       Purpose: Compute exponential integral Ei(x)
%       Input :  x  --- Argument of Ei(x)
%       Output:  EI --- Ei(x,x > 0)
%       ============================================
if(x == 0.0);
ei=-1.0d+300;
elseif(x <= 40.0);
ei=1.0d0;
r=1.0d0;
for  k=1:100;
r=r.*k.*x./(k+1.0d0).^2;
ei=ei+r;
if(abs(r./ei)<= 1.0d-15)break; end;
end;
ga=0.5772156649015328d0;
ei=ga+log(x)+x.*ei;
else;
ei=1.0d0;
r=1.0d0;
for  k=1:20;
r=r.*k./x;
ei=ei+r;
end;  k=20+1;
ei=exp(x)./x.*ei;
end;
return;
end

