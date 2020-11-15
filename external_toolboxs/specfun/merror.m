function merror
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===================================================
%       Purpose: This program computes the error function
%                erf(x)using subroutine ERROR
%       Input:   x   --- Argument of erf(x)
%       Output:  ERR --- erf(x)
%       Example:
%                  x         erf(x)
%                ---------------------
%                 1.0       .84270079
%                 2.0       .99532227
%                 3.0       .99997791
%                 4.0       .99999998
%                 5.0      1.00000000
%       ===================================================
x=[];er=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=1.0;
fprintf(1,'%s \n','   x         erf(x)');
fprintf(1,'%s \n','---------------------');
[x,er]=errorf(x,er);
fprintf(1,[repmat(' ',1,1),'%5.2g','%15.8g' ' \n'],x,er);
%format(1x,f5.2,f15.8);
end
function [x,err]=errorf(x,err,varargin);
%       =========================================
%       Purpose: Compute error function erf(x)
%       Input:   x   --- Argument of erf(x)
%       Output:  ERR --- erf(x)
%       =========================================
eps=1.0d-15;
pi=3.141592653589793d0;
x2=x.*x;
if(abs(x)< 3.5d0);
er=1.0d0;
r=1.0d0;
for  k=1:50;
r=r.*x2./(k+0.5d0);
er=er+r;
if(abs(r)<= abs(er).*eps)break; end;
end;
c0=2.0d0./sqrt(pi).*x.*exp(-x2);
err=c0.*er;
else;
er=1.0d0;
r=1.0d0;
for  k=1:12;
r=-r.*(k-0.5d0)./x2;
er=er+r;
end;  k=12+1;
c0=exp(-x2)./(abs(x).*sqrt(pi));
err=1.0d0-c0.*er;
if(x < 0.0)err=-err; end;
end;
return;
end

