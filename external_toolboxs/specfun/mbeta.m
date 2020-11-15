function mbeta
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ====================================================
%       Purpose: This program computes the beta function
%                B(p,q)for p > 0 and q > 0 using
%                subroutine BETA
%       Input :  p  --- Parameter(p > 0)
%                q  --- Parameter(q > 0)
%       Output:  BT --- B(p,q)
%       Examples:
%                 p       q           B(p,q)
%               ---------------------------------
%                1.5     2.0     .2666666667D+00
%                2.5     2.0     .1142857143D+00
%                1.5     3.0     .1523809524D+00
%       ====================================================
p=[];q=[];bt=[];
fprintf(1,'%s \n','please enter p and q');
%        READ(*,*)P,Q
p=2.5;
q=2.0;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','    p       q           b(p,q)');
fprintf(1,'%s \n','  ---------------------------------');
[p,q,bt]=beta(p,q,bt);
fprintf(1,[repmat(' ',1,2),'%5.1g',repmat(' ',1,3),'%5.1g','%20.10g' ' \n'],p,q,bt);
%format(2x,f5.1,3x,f5.1,d20.10);
end
function [p,q,bt]=beta(p,q,bt,varargin);
%       ==========================================
%       Purpose: Compute the beta function B(p,q)
%       Input :  p  --- Parameter(p > 0)
%                q  --- Parameter(q > 0)
%       Output:  BT --- B(p,q)
%       Routine called: GAMMA for computing ג(x)
%       ==========================================
gp=[];gq=[];ppq=[];gpq=[];
[p,gp]=gamma(p,gp);
[q,gq]=gamma(q,gq);
ppq=p+q;
[ppq,gpq]=gamma(ppq,gpq);
bt=gp.*gq./gpq;
return;
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

