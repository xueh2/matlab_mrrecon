function menxa
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
%                En(x)using subroutine ENXA
%       Example: x = 10.0
%                   n         En(x)
%                 ----------------------
%                   0     .45399930D-05
%                   1     .41569689D-05
%                   2     .38302405D-05
%                   3     .35487626D-05
%                   4     .33041014D-05
%                   5     .30897289D-05
%       =========================================================
n=[];x=[];en=[];
 en=zeros(1,100+1);
fprintf(1,'%s \n','please enter n and x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,5),'%3g',',   ','x=','%5.1g' ' \n'],n,x);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   n         en(x)');
fprintf(1,'%s \n',' ----------------------');
[n,x,en]=enxa(n,x,en);
for  k=0:n;
fprintf(1,[repmat(' ',1,2),'%3g','%18.8g' ' \n'],k,en(k+1));
end;  k=n+1;
%format(5x,i3,',   ',',f5.1);
%format(2x,i3,d18.8);
end
function [n,x,en]=enxa(n,x,en,varargin);
%       ============================================
%       Purpose: Compute exponential integral En(x)
%       Input :  x --- Argument of En(x,x ó 20)
%                n --- Order of En(x)
%       Output:  EN(n)--- En(x)
%       Routine called: E1XB for computing E1(x)
%       ============================================
e1=[];
en(0+1)=exp(-x)./x;
[x,e1]=e1xb(x,e1);
en(1+1)=e1;
for  k=2:n;
ek=(exp(-x)-x.*e1)./(k-1.0d0);
en(k+1)=ek;
e1=ek;
end;  k=fix(n)+1;
return;
end
function [x,e1]=e1xb(x,e1,varargin);
%       ============================================
%       Purpose: Compute exponential integral E1(x)
%       Input :  x  --- Argument of E1(x)
%       Output:  E1 --- E1(x)
%       ============================================
if(x == 0.0);
e1=1.0d+300;
elseif(x <= 1.0);
e1=1.0d0;
r=1.0d0;
for  k=1:25;
r=-r.*k.*x./(k+1.0d0).^2;
e1=e1+r;
if(abs(r)<= abs(e1).*1.0d-15)break; end;
end;
ga=0.5772156649015328d0;
e1=-ga-log(x)+x.*e1;
else;
m=20+fix(80.0./x);
t0=0.0d0;
for  k=m:-1:1;
t0=k./(1.0d0+k./(x+t0));
end;  k=1-1;
t=1.0d0./(x+t0);
e1=exp(-x).*t;
end;
return;
end

