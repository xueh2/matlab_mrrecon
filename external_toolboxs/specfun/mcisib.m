function mcisib
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes the cosine and sine
%                integrals using subroutine CISIB
%       Input :  x  --- Argument of Ci(x)and Si(x)
%       Output:  CI --- Ci(x)
%                SI --- Si(x)
%       Example:
%                   x        Ci(x)Si(x)
%                ------------------------------------
%                  0.0    - ì                 0
%                  5.0    -.190030D+00      1.549931
%                 10.0    -.454563D-01      1.658348
%                 20.0     .444201D-01      1.548241
%                 30.0    -.330326D-01      1.566757
%                 40.0     .190201D-01      1.586985
%       ========================================================
x=[];ci=[];si=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=10.0;
fprintf(1,'%s \n','   x        ci(x)si(x)');
fprintf(1,'%s \n','------------------------------------');
[x,ci,si]=cisib(x,ci,si);
if(x ~= 0.0d0)fprintf(1,[repmat(' ',1,1),'%5.1g','%16.6g','%14.6g' ' \n'],x,ci,si); end;
if(x == 0.0d0)fprintf(1,[repmat(' ',1,3),' .0',repmat(' ',1,3),' - ì',repmat(' ',1,17),'0' ' \n']); end;
%format(1x,f5.1,d16.6,f14.6);
%format(3x,' .0',3x,' - ì',17x,'0');
end
function [x,ci,si]=cisib(x,ci,si,varargin);
%       =============================================
%       Purpose: Compute cosine and sine integrals
%                Si(x)and Ci(x,x ò 0)
%       Input :  x  --- Argument of Ci(x)and Si(x)
%       Output:  CI --- Ci(x)
%                SI --- Si(x)
%       =============================================
x2=x.*x;
if(x == 0.0);
ci=-1.0d+300;
si=0.0d0;
elseif(x <= 1.0d0);
ci=((((-3.0d-8.*x2+3.10d-6).*x2-2.3148d-4).*x2+1.041667d-2).*x2-0.25).*x2+0.577215665d0+log(x);
si=((((3.1d-7.*x2-2.834d-5).*x2+1.66667d-003).*x2-5.555556d-002).*x2+1.0).*x;
else;
fx=((((x2+38.027264d0).*x2+265.187033d0).*x2+335.67732d0).*x2+38.102495d0)./((((x2 +40.021433d0).*x2+322.624911d0).*x2+570.23628d0).*x2+157.105423d0);
gx=((((x2+42.242855d0).*x2+302.757865d0).*x2+352.018498d0).*x2+21.821899d0)./((((x2 +48.196927d0).*x2+482.485984d0).*x2+1114.978885d0).*x2+449.690326d0)./x;
ci=fx.*sin(x)./x-gx.*cos(x)./x;
si=1.570796327d0-fx.*cos(x)./x-gx.*sin(x)./x;
end;
return;
end

