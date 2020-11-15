function mincob
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes the incomplete beta
%                function Ix(a,b)using subroutine INCOB
%       Input :  a --- Parameter
%                b --- Parameter
%                x --- Argument(0 ף x ף 1)
%       Output:  BIX --- Ix(a,b)
%       Example:
%                  a       b       x       Ix(a,b)
%                -----------------------------------
%                 1.0     3.0     .25     .57812500
%       =========================================================
a=[];b=[];x=[];bix=[];
fprintf(1,'%s \n','please enter a, b and x(0 ף x ף 1)');
%        READ(*,*)A,B,X
a=1.0;
b=3.0;
x=.25;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   a       b       x       ix(a,b)');
fprintf(1,'%s \n',' -----------------------------------');
[a,b,x,bix]=incob(a,b,x,bix);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,3),'%5.1g',repmat(' ',1,3),'%5.2g','%14.8g' ' \n'],a,b,x,bix);
%format(1x,f5.1,3x,f5.1,3x,f5.2,f14.8);
end
function [a,b,x,bix]=incob(a,b,x,bix,varargin);
%       ========================================================
%       Purpose: Compute the incomplete beta function Ix(a,b)
%       Input :  a --- Parameter
%                b --- Parameter
%                x --- Argument(0 ף x ף 1)
%       Output:  BIX --- Ix(a,b)
%       Routine called: BETA for computing beta function B(p,q)
%       ========================================================
bt=[];
 dk=zeros(1,51);
fk=zeros(1,51);
s0=(a+1.0d0)./(a+b+2.0d0);
[a,b,bt]=beta(a,b,bt);
if(x <= s0);
for  k=1:20;
dk(2.*k)=k.*(b-k).*x./(a+2.0d0.*k-1.0d0)./(a+2.0d0.*k);
end;  k=20+1;
for  k=0:20;
dk(2.*k+1)=-(a+k).*(a+b+k).*x./(a+2.d0.*k)./(a+2.0.*k+1.0);
end;  k=20+1;
t1=0.0d0;
for  k=20:-1:1;
t1=dk(k)./(1.0d0+t1);
end;  k=1-1;
ta=1.0d0./(1.0d0+t1);
bix=x.^a.*(1.0d0-x).^b./(a.*bt).*ta;
else;
for  k=1:20;
fk(2.*k)=k.*(a-k).*(1.0d0-x)./(b+2..*k-1.0)./(b+2.0.*k);
end;  k=20+1;
for  k=0:20;
fk(2.*k+1)=-(b+k).*(a+b+k).*(1.d0-x)./(b+2.d0.*k)./(b+2.d0.*k+1.d0);
end;  k=20+1;
t2=0.0d0;
for  k=20:-1:1;
t2=fk(k)./(1.0d0+t2);
end;  k=1-1;
tb=1.0d0./(1.0d0+t2);
bix=1.0d0-x.^a.*(1.0d0-x).^b./(b.*bt).*tb;
end;
return;
end
function [p,q,bt]=beta(p,q,bt,varargin);
%       ==========================================
%       Purpose: Compute the beta function B(p,q)
%       Input :  p --- Parameter(p > 0)
%                q --- Parameter(q > 0)
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

