function me1xb
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
%                E1(x)using subroutine E1XB
%       Input :  x  --- Argument of E1(x,x > 0)
%       Output:  E1 --- E1(x)
%       Example:
%                  x          E1(x)
%                -------------------------
%                 0.0     .1000000000+301
%                 1.0     .2193839344E+00
%                 2.0     .4890051071E-01
%                 3.0     .1304838109E-01
%                 4.0     .3779352410E-02
%                 5.0     .1148295591E-02
%       =========================================================
x=[];e1=[];
  e1=0;
 x=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=5.0;
fprintf(1,'%s \n','   x          e1(x)');
fprintf(1,'%s \n',' -------------------------');
[x,e1]=e1xb(x,e1);
fprintf(1,[repmat(' ',1,1),'%5.1g','%20.10g' ' \n'],x,e1);
%format(1x,f5.1,e20.10);
end
function [x,e1]=e1xb(x,e1,varargin);
%       ============================================
%       Purpose: Compute exponential integral E1(x)
%       Input :  x  --- Argument of E1(x)
%       Output:  E1 --- E1(x,x > 0)
%       ============================================
if(x == 0.0);
e1=1.0d+300;
elseif(x <= 1.0);
e1=1.0d0;
r=1.0d0;
for  k=1:25;
r=-r.*k.*x./(k+1.0d0).^2;
e1=e1+r;
if(abs(r)<= abs(e1).*1.0d-15)break; end;
end;
ga=0.5772156649015328d0;
e1=-ga-log(x)+x.*e1;
else;
m=20+fix(80.0./x);
t0=0.0d0;
for  k=m:-1:1;
t0=k./(1.0d0+k./(x+t0));
end;  k=1-1;
t=1.0d0./(x+t0);
e1=exp(-x).*t;
end;
return;
end

