function menxb
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
%                En(x)using subroutine ENXB
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
[n,x,en]=enxb(n,x,en);
for  k=0:n;
fprintf(1,[repmat(' ',1,2),'%3g','%18.8g' ' \n'],k,en(k+1));
end;  k=n+1;
%format(5x,i3,',   ',',f5.1);
%format(2x,i3,d18.8);
end
function [n,x,en]=enxb(n,x,en,varargin);
%       ===============================================
%       Purpose: Compute exponential integral En(x)
%       Input :  x --- Argument of En(x)
%                n --- Order of En(x,n = 0,1,2,...)
%       Output:  EN(n)--- En(x)
%       ===============================================
if(x == 0.0);
en(0+1)=1.0d+300;
en(1+1)=1.0d+300;
for  k=2:n;
en(k+1)=1.0d0./(k-1.0);
end;  k=fix(n)+1;
return;
elseif(x <= 1.0);
en(0+1)=exp(-x)./x;
for  l=1:n;
rp=1.0d0;
for  j=1:l-1;
rp=-rp.*x./j;
end;  j=l-1+1;
ps=-0.5772156649015328d0;
for  m=1:l-1;
ps=ps+1.0d0./m;
end;  m=l-1+1;
ens=rp.*(-log(x)+ps);
s=0.0d0;
for  m=0:20;
if(~(m == l-1));
r=1.0d0;
for  j=1:m;
r=-r.*x./j;
end;  j=m+1;
s=s+r./(m-l+1.0d0);
if(abs(s-s0)< abs(s).*1.0d-15)break; end;
s0=s;
end;
end;
en(l+1)=ens-s;
end;
else;
en(0+1)=exp(-x)./x;
m=15+fix(100.0./x);
for  l=1:n;
t0=0.0d0;
for  k=m:-1:1;
t0=(l+k-1.0d0)./(1.0d0+k./(x+t0));
end;  k=1-1;
t=1.0d0./(x+t0);
en(l+1)=exp(-x).*t;
end;  l=fix(n)+1;
end;
end

