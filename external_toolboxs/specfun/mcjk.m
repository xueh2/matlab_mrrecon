function mcjk
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes the expansion coefficients
%                for the asymptotic expansion of Bessel functions
%                with large orders using subroutine CJK
%       Input :  Km   --- Maximum k
%       Output:  A(L)--- Cj(k)where j and k are related to L by
%                         L=j+1+[k*(k+1)]/2; j,k=0,1,2,...,Km
%       ============================================================
km=[];a=[];
 a=zeros(1,231);
fprintf(1,'%s \n','please enter km(ó 20)');
%        READ(*,*)KM
km=2;
lm=km+1+(km.*(km+1))./2;
[km,a]=cjk(km,a);
for  k=1:lm;
fprintf(1,[repmat(' ',1,1),'%3g','%25.14g' ' \n'],k,a(k));
end;  k=lm+1;
%format(1x,i3,d25.14);
end
function [km,a]=cjk(km,a,varargin);
%       ========================================================
%       Purpose: Compute the expansion coefficients for the
%                asymptotic expansion of Bessel functions
%                with large orders
%       Input :  Km   --- Maximum k
%       Output:  A(L)--- Cj(k)where j and k are related to L
%                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
%       ========================================================
a(1)=1.0d0;
f0=1.0d0;
g0=1.0d0;
for  k=0:km-1;
l1=(k+1).*(k+2)./2+1;
l2=(k+1).*(k+2)./2+k+2;
f=(0.5d0.*k+0.125d0./(k+1)).*f0;
g=-(1.5d0.*k+0.625d0./(3.0.*(k+1.0d0))).*g0;
a(l1)=f;
a(l2)=g;
f0=f;
g0=g;
end;  k=fix(km)-1+1;
for  k=1:km-1;
for  j=1:k;
l3=k.*(k+1)./2+j+1;
l4=(k+1).*(k+2)./2+j+1;
a(l4)=(j+0.5d0.*k+0.125d0./(2.0.*j+k+1.0)).*a(l3)-(j+0.5d0.*k-1.0+0.625d0./(2.0.*j+k+1.0)).*a(l3-1);
end;  j=k+1;
end;  k=fix(km)-1+1;
return;
end

