function mitsh0
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ====================================================
%       Purpose: This program evaluates the integral of
%                Struve function H0(t)with respect to t
%                from 0 and x using subroutine ITSH0
%       Input :  x   --- Upper limit(x ò 0)
%       Output:  TH0 --- Integration of H0(t)from 0 and x
%       Example:
%                    x        H0(t)dt
%                 ----------------------
%                   0.0       .0000000
%                   5.0      2.0442437
%                  10.0      2.5189577
%                  15.0      2.5415824
%                  20.0      2.5484517
%                  30.0      3.0625848
%                  40.0      3.1484123
%                  50.0      3.2445168
%       ====================================================
x=[];th0=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=50.0;
fprintf(1,'%s \n','   x         h0(t)dt');
fprintf(1,'%s \n','----------------------ä');
[x,th0]=itsh0(x,th0);
fprintf(1,[repmat(' ',1,1),'%5.1g','%16.7g' ' \n'],x,th0);
%format(1x,f5.1,e16.7);
end
function [x,th0]=itsh0(x,th0,varargin);
%       ===================================================
%       Purpose: Evaluate the integral of Struve function
%                H0(t)with respect to t from 0 and x
%       Input :  x   --- Upper limit(x ò 0)
%       Output:  TH0 --- Integration of H0(t)from 0 and x
%       ===================================================
 a=zeros(1,25);
pi=3.141592653589793d0;
r=1.0d0;
if(x <= 30.0);
s=0.5d0;
for  k=1:100;
rd=1.0d0;
if(k == 1)rd=0.5d0; end;
r=-r.*rd.*k./(k+1.0d0).*(x./(2.0d0.*k+1.0d0)).^2;
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
th0=2.0d0./pi.*x.*x.*s;
else;
s=1.0d0;
for  k=1:12;
r=-r.*k./(k+1.0d0).*((2.0d0.*k+1.0d0)./x).^2;
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
el=.57721566490153d0;
s0=s./(pi.*x.*x)+2.0d0./pi.*(log(2.0d0.*x)+el);
a0=1.0d0;
a1=5.0d0./8.0d0;
a(1)=a1;
for  k=1:20;
af=((1.5d0.*(k+.5d0).*(k+5.0d0./6.0d0).*a1-.5d0.*(k+.5d0).*(k+.5d0).*(k-.5d0).*a0))./(k+1.0d0);
a(k+1)=af;
a0=a1;
a1=af;
end;  k=20+1;
bf=1.0d0;
r=1.0d0;
for  k=1:10;
r=-r./(x.*x);
bf=bf+a(2.*k).*r;
end;  k=10+1;
bg=a(1)./x;
r=1.0d0./x;
for  k=1:10;
r=-r./(x.*x);
bg=bg+a(2.*k+1).*r;
end;  k=10+1;
xp=x+.25d0.*pi;
ty=sqrt(2.0d0./(pi.*x)).*(bg.*cos(xp)-bf.*sin(xp));
th0=ty+s0;
end;
return;
end

