function mstvl1
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =====================================================
%       Purpose: This program computes the modified Struve
%                function L1(x)using subroutine STVL1
%       Input :  x   --- Argument of L1(x,x ò 0)
%       Output:  SL1 --- L1(x)
%       Example:
%                     x        L1(x)
%                 -----------------------
%                   0.0   .00000000D+00
%                   5.0   .23728216D+02
%                  10.0   .26703583D+04
%                  15.0   .32812429D+06
%                  20.0   .42454973D+08
%                  30.0   .76853204D+12
%                  40.0   .14707396D+17
%                  50.0   .29030786D+21
%       =====================================================
x=[];sl1=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=50.0;
fprintf(1,'%s \n','   x        l1(x)');
fprintf(1,'%s \n','-----------------------');
[x,sl1]=stvl1(x,sl1);
fprintf(1,[repmat(' ',1,1),'%5.1g','%16.8g' ' \n'],x,sl1);
%format(1x,f5.1,d16.8);
end
function [x,sl1]=stvl1(x,sl1,varargin);
%       ================================================
%       Purpose: Compute modified Struve function L1(x)
%       Input :  x   --- Argument of L1(x,x ò 0)
%       Output:  SL1 --- L1(x)
%       ================================================
pi=3.141592653589793d0;
r=1.0d0;
if(x <= 20.0d0);
s=0.0d0;
for  k=1:60;
r=r.*x.*x./(4.0d0.*k.*k-1.0d0);
s=s+r;
if(abs(r)< abs(s).*1.0d-12)break; end;
end;
sl1=2.0d0./pi.*s;
else;
s=1.0d0;
km=fix(.50.*x);
if(x > 50)km=25; end;
for  k=1:km;
r=r.*(2.0d0.*k+3.0d0).*(2.0d0.*k+1.0d0)./(x.*x);
s=s+r;
if(abs(r./s)< 1.0d-12)break; end;
end;
sl1=2.0d0./pi.*(-1.0d0+1.0d0./(x.*x)+3.0d0.*s./x.^4);
a1=exp(x)./sqrt(2.0d0.*pi.*x);
r=1.0d0;
bi1=1.0d0;
for  k=1:16;
r=-0.125d0.*r.*(4.0d0-(2.0d0.*k-1.0d0).^2)./(k.*x);
bi1=bi1+r;
if(abs(r./bi1)< 1.0d-12)break; end;
end;
sl1=sl1+a1.*bi1;
end;
return;
end

