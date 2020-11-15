function mitth0
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program evaluates the integral of H0(t)/t
%                with respect to t from x to infinity using
%                subroutine ITTH0
%       Input :  x   --- Lower limit(x ò 0)
%       Output:  TTH --- Integration of H0(t)/t from x to infinity
%       Example:
%                    x        H0(t)/t dt
%                 -----------------------
%                   0.0      1.57079633
%                   5.0       .07954575
%                  10.0       .04047175
%                  15.0       .04276558
%                  20.0       .04030796
%                  30.0       .01815256
%                  40.0       .01621331
%                  50.0       .01378661
%       =======================================================
x=[];tth=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=50.0;
fprintf(1,'%s \n','   x       h0(t)/t dt');
fprintf(1,'%s \n','-----------------------');
[x,tth]=itth0(x,tth);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,1),'%16.8g' ' \n'],x,tth);
%format(1x,f5.1,1x,e16.8);
end
function [x,tth]=itth0(x,tth,varargin);
%       ===========================================================
%       Purpose: Evaluate the integral H0(t)/t with respect to t
%                from x to infinity
%       Input :  x   --- Lower limit(x ò 0)
%       Output:  TTH --- Integration of H0(t)/t from x to infinity
%       ===========================================================
pi=3.141592653589793d0;
s=1.0d0;
r=1.0d0;
if(x < 24.5d0);
for  k=1:60;
r=-r.*x.*x.*(2.0.*k-1.0d0)./(2.0.*k+1.0d0).^3;
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
tth=pi./2.0d0-2.0d0./pi.*x.*s;
else;
for  k=1:10;
r=-r.*(2.0.*k-1.0d0).^3./((2.0.*k+1.0d0).*x.*x);
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
tth=2.0d0./(pi.*x).*s;
t=8.0d0./x;
xt=x+.25d0.*pi;
f0=(((((.18118d-2.*t-.91909d-2).*t+.017033d0).*t-.9394d-3).*t-.051445d0).*t-.11d-5).*t+.7978846d0;
g0=(((((-.23731d-2.*t+.59842d-2).*t+.24437d-2).*t-.0233178d0).*t+.595d-4).*t+.1620695d0).*t;
tty=(f0.*sin(xt)-g0.*cos(xt))./(sqrt(x).*x);
tth=tth+tty;
end;
return;
end

