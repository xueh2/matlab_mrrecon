function mitsl0
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program evaluates the integral of modified
%                Struve function L0(t)with respect to t from 0
%                to x using subroutine ITSL0
%       Input :  x   --- Upper limit(x ò 0)
%       Output:  TL0 --- Integration of L0(t)from 0 to x
%       Example:
%                      x        L0(t)dt
%                   -----------------------
%                     0.0    .0000000D+00
%                     5.0    .3003079D+02
%                    10.0    .2990773D+04
%                    15.0    .3526179D+06
%                    20.0    .4475860D+08
%                    30.0    .7955389D+12
%                    40.0    .1508972D+17
%                    50.0    .2962966D+21
%       ===========================================================
x=[];tl0=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=50.0;
fprintf(1,'%s \n','   x        l0(t)dt');
fprintf(1,'%s \n','-----------------------');
[x,tl0]=itsl0(x,tl0);
fprintf(1,[repmat(' ',1,1),'%5.1g','%16.7g' ' \n'],x,tl0);
%format(1x,f5.1,d16.7);
end
function [x,tl0]=itsl0(x,tl0,varargin);
%       ===========================================================
%       Purpose: Evaluate the integral of modified Struve function
%                L0(t)with respect to t from 0 to x
%       Input :  x   --- Upper limit(x ò 0)
%       Output:  TL0 --- Integration of L0(t)from 0 to x
%       ===========================================================
 a=zeros(1,18);
pi=3.141592653589793d0;
r=1.0d0;
if(x <= 20.0);
s=0.5d0;
for  k=1:100;
rd=1.0d0;
if(k == 1)rd=0.5d0; end;
r=r.*rd.*k./(k+1.0d0).*(x./(2.0d0.*k+1.0d0)).^2;
s=s+r;
if(abs(r./s)< 1.0d-12)break; end;
end;
tl0=2.0d0./pi.*x.*x.*s;
else;
s=1.0d0;
for  k=1:10;
r=r.*k./(k+1.0d0).*((2.0d0.*k+1.0d0)./x).^2;
s=s+r;
if(abs(r./s)< 1.0d-12)break; end;
end;
el=.57721566490153d0;
s0=-s./(pi.*x.*x)+2.0d0./pi.*(log(2.0d0.*x)+el);
a0=1.0d0;
a1=5.0d0./8.0d0;
a(1)=a1;
for  k=1:10;
af=((1.5d0.*(k+.50d0).*(k+5.0d0./6.0d0).*a1-.5d0.*(k+.5d0).^2.*(k-.5d0).*a0))./(k+1.0d0);
a(k+1)=af;
a0=a1;
a1=af;
end;  k=10+1;
ti=1.0d0;
r=1.0d0;
for  k=1:11;
r=r./x;
ti=ti+a(k).*r;
end;  k=11+1;
tl0=ti./sqrt(2.*pi.*x).*exp(x)+s0;
end;
return;
end

