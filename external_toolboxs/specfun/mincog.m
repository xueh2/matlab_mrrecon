function mincog
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes the incomplete gamma
%                function r(a,x), ג(a,x)and P(a,x)using
%                subroutine INCOG
%       Input :  a   --- Parameter
%                x   --- Argument
%       Output:  GIN --- r(a,x)
%                GIM --- ג(a,x)
%                GIP --- P(a,x)
%       Example:
%            a     x      r(a,x)ג(a,x)P(a,x)
%           -------------------------------------------------------
%           3.0   5.0  .17506960D+01  .24930404D+00  .87534798D+00
%       ===========================================================
a=[];x=[];gin=[];gim=[];gip=[];
  a=0;
 x=0;
 gin=0;
 gim=0;
 gip=0;
fprintf(1,'%s \n','plese enter a and x');
%        READ(*,*)A,X
a=3.0;
x=5.0;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   a     x      r(a,x)g(a,x)p(a,x)');
fprintf(1,'%s ',' --------------------------------------------');fprintf(1,'%s \n', '------------');
[a,x,gin,gim,gip]=incog(a,x,gin,gim,gip);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,1),'%5.1g',repmat('%15.8g',1,3) ' \n'],a,x,gin,gim,gip);
%format(1x,f5.1,1x,f5.1,3d15.8);
end
function [a,x,gin,gim,gip]=incog(a,x,gin,gim,gip,varargin);
%       ===================================================
%       Purpose: Compute the incomplete gamma function
%                r(a,x), ג(a,x)and P(a,x)
%       Input :  a   --- Parameter(a ף 170)
%                x   --- Argument
%       Output:  GIN --- r(a,x)
%                GIM --- ג(a,x)
%                GIP --- P(a,x)
%       Routine called: GAMMA for computing ג(x)
%       ===================================================
ga=[];
xam=-x+a.*log(x);
if(xam > 700.0|a > 170.0);
fprintf(1,'%s \n','a and/or x too large');
error('stop encountered in original fortran code');
end;
if(x == 0.0);
gin=0.0;
[a,ga]=gamma(a,ga);
gim=ga;
gip=0.0;
elseif(x <= 1.0+a);
s=1.0d0./a;
r=s;
for  k=1:60;
r=r.*x./(a+k);
s=s+r;
if(abs(r./s)< 1.0d-15)break; end;
end;
gin=exp(xam).*s;
[a,ga]=gamma(a,ga);
gip=gin./ga;
gim=ga-gin;
elseif(x > 1.0+a);
t0=0.0d0;
for  k=60:-1:1;
t0=(k-a)./(1.0d0+k./(x+t0));
end;  k=1-1;
gim=exp(xam)./(x+t0);
[a,ga]=gamma(a,ga);
gin=ga-gim;
gip=1.0d0-gim./ga;
end;
end
function [x,ga]=gamma(x,ga,varargin);
%       ==================================================
%       Purpose: Compute gamma function ג(x)
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

