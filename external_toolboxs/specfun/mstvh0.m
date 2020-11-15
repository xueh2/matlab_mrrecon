function mstvh0
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ====================================================
%       Purpose: This program computes Struve function
%                H0(x)using subroutine STVH0
%       Input :  x   --- Argument of H0(x,x ò 0)
%       Output:  SH0 --- H0(x)
%       Example:
%                   x          H0(x)
%                ----------------------
%                  0.0       .00000000
%                  5.0      -.18521682
%                 10.0       .11874368
%                 15.0       .24772383
%                 20.0       .09439370
%                 25.0      -.10182519
%       ====================================================
x=[];sh0=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x          h0(x)');
fprintf(1,'%s \n','----------------------');
[x,sh0]=stvh0(x,sh0);
fprintf(1,[repmat(' ',1,1),'%5.1g','%16.8g' ' \n'],x,sh0);
%format(1x,f5.1,e16.8);
end
function [x,sh0]=stvh0(x,sh0,varargin);
%       =============================================
%       Purpose: Compute Struve function H0(x)
%       Input :  x   --- Argument of H0(x,x ò 0)
%       Output:  SH0 --- H0(x)
%       =============================================
pi=3.141592653589793d0;
s=1.0d0;
r=1.0d0;
if(x <= 20.0d0);
a0=2.0.*x./pi;
for  k=1:60;
r=-r.*x./(2.0d0.*k+1.0d0).*x./(2.0d0.*k+1.0d0);
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
sh0=a0.*s;
else;
km=fix(.5.*(x+1.0));
if(x >= 50.0)km=25; end;
for  k=1:km;
r=-r.*((2.0d0.*k-1.0d0)./x).^2;
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
t=4.0d0./x;
t2=t.*t;
p0=((((-.37043d-5.*t2+.173565d-4).*t2-.487613d-4).*t2+.17343d-3).*t2-.1753062d-2).*t2+.3989422793d0;
q0=t.*(((((.32312d-5.*t2-.142078d-4).*t2+.342468d-4).*t2-.869791d-4).*t2+.4564324d-3).*t2-.0124669441d0);
ta0=x-.25d0.*pi;
by0=2.0d0./sqrt(x).*(p0.*sin(ta0)+q0.*cos(ta0));
sh0=2.0d0./(pi.*x).*s+by0;
end;
return;
end

