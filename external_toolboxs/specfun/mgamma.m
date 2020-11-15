function mgamma
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ====================================================
%       Purpose: This program computes the gamma function
%                ג(x)using subroutine GAMMA
%       Examples:
%                   x            ג(x)
%                ----------------------------
%                  1/3       2.678938534708
%                  0.5       1.772453850906
%                 -0.5      -3.544907701811
%                 -1.5       2.363271801207
%                  5.0      24.000000000000
%       ====================================================
x=[];ga=[];
  a=0;
 x=0;
 ga=0;
 a=zeros(1,5);
a(:)=[.333333333333333333d0,0.5d0,-0.5d0,-1.5,5.0d0];
fprintf(1,'%s \n','     x            ג(x)');
fprintf(1,'%s \n',' ----------------------------');
for  k=1:5;
x=a(k);
[x,ga]=gamma(x,ga);
fprintf(1,[repmat(' ',1,1),'%8.4g','%20.12g' ' \n'],x,ga);
end;  k=5+1;
fprintf(1,'%s \n','please enter x:');
%        READ(*,*)X
x=.5;
[x,ga]=gamma(x,ga);
fprintf(1,[repmat(' ',1,1),'%8.4g','%20.12g' ' \n'],x,ga);
%format(1x,f8.4,e20.12);
end
function [x,ga]=gamma(x,ga,varargin);
%       ==================================================
%       Purpose: Compute the gamma function ג(x)
%       Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%       Output:  GA --- ג(x)
%       ==================================================
 g=zeros(1,26);
pi=3.141592653589793d0;
if(x == fix(x));
if(x > 0.0d0);
ga=1.0d0;
m1=x-1;
for  k=2:m1;
ga=ga.*k;
end;  k=m1+1;
else;
ga=1.0d+300;
end;
else;
if(abs(x)> 1.0d0);
z=abs(x);
m=fix(z);
r=1.0d0;
for  k=1:m;
r=r.*(z-k);
end;  k=m+1;
z=z-m;
else;
z=x;
end;
g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
gr=g(26);
for  k=25:-1:1;
gr=gr.*z+g(k);
end;  k=1-1;
ga=1.0d0./(gr.*z);
if(abs(x)> 1.0d0);
ga=ga.*r;
if(x < 0.0d0)ga=-pi./(x.*ga.*sin(pi.*x)); end;
end;
end;
return;
end

