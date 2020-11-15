function me1xa
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
%                E1(x)using subroutine E1XA
%       Input :  x  --- Argument of E1(x,x > 0)
%       Output:  E1 --- E1(x)
%       Example:
%                  x        E1(x)
%                ----------------------
%                 0.0     .1000000+301
%                 1.0     .2193839E+00
%                 2.0     .4890051E-01
%                 3.0     .1304838E-01
%                 4.0     .3779352E-02
%                 5.0     .1148296E-02
%       =========================================================
x=[];e1=[];
  e1=0;
 x=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=5.0;
fprintf(1,'%s \n','   x        e1(x)');
fprintf(1,'%s \n',' ----------------------');
[x,e1]=e1xa(x,e1);
fprintf(1,[repmat(' ',1,1),'%5.1g','%17.7g' ' \n'],x,e1);
%format(1x,f5.1,e17.7);
end
function [x,e1]=e1xa(x,e1,varargin);
%       ============================================
%       Purpose: Compute exponential integral E1(x)
%       Input :  x  --- Argument of E1(x)
%       Output:  E1 --- E1(x,x > 0)
%       ============================================
if(x == 0.0);
e1=1.0d+300;
elseif(x <= 1.0);
e1=-log(x)+((((1.07857d-3.*x-9.76004d-3).*x+5.519968d-2).*x-0.24991055d0).*x+0.99999193d0).*x-0.57721566d0;
else;
es1=(((x+8.5733287401d0).*x+18.059016973d0).*x+8.6347608925d0).*x+0.2677737343d0;
es2=(((x+9.5733223454d0).*x+25.6329561486d0).*x+21.0996530827d0).*x+3.9584969228d0;
e1=exp(-x)./x.*es1./es2;
end;
return;
end

