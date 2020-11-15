function msphk
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ======================================================
%       Purpose: This program computes the modified spherical
%                Bessel functions kn(x)and kn'(x)using
%                subroutine SPHK
%       Input :  x --- Argument of kn(x,x ò 0)
%                n --- Order of kn(x,n ó 250)
%       Output:  SK(n)--- kn(x)
%                DK(n)--- kn'(x)
%       Example: x= 10.0
%                  n          kn(x)kn'(x)
%                --------------------------------------------
%                  0     .7131404291D-05    -.7844544720D-05
%                  1     .7844544720D-05    -.8700313235D-05
%                  2     .9484767707D-05    -.1068997503D-04
%                  3     .1258692857D-04    -.1451953914D-04
%                  4     .1829561771D-04    -.2173473743D-04
%                  5     .2905298451D-04    -.3572740841D-04
%       ======================================================
n=[];x=[];nm=[];sk=[];dk=[];
 sk=zeros(1,250+1);
dk=zeros(1,250+1);
fprintf(1,'%s \n','please enter n and x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',     ','x =','%6.1g' ' \n'],n,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,x,nm,sk,dk]=sphk(n,x,nm,sk,dk);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n          kn(x)kn''(x)');
fprintf(1,'%s \n','--------------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%20.10g',1,2) ' \n'],k,sk(k+1),dk(k+1));
end;  k=nm+1;
%format(1x,i3,2d20.10);
%format(3x,',i3,',     ',',f6.1);
end
function [n,x,nm,sk,dk]=sphk(n,x,nm,sk,dk,varargin);
%       =====================================================
%       Purpose: Compute modified spherical Bessel functions
%                of the second kind, kn(x)and kn'(x)
%       Input :  x --- Argument of kn(x,x ò 0)
%                n --- Order of kn(x,n = 0,1,2,...)
%       Output:  SK(n)--- kn(x)
%                DK(n)--- kn'(x)
%                NM --- Highest order computed
%       =====================================================
pi=3.141592653589793d0;
nm=fix(fix(n));
if(x < 1.0d-60);
for  k=0:n;
sk(k+1)=1.0d+300;
dk(k+1)=-1.0d+300;
end;  k=fix(n)+1;
return;
end;
sk(0+1)=0.5d0.*pi./x.*exp(-x);
sk(1+1)=sk(0+1).*(1.0d0+1.0d0./x);
f0=sk(0+1);
f1=sk(1+1);
for  k=2:n;
f=(2.0d0.*k-1.0d0).*f1./x+f0;
sk(k+1)=f;
if(abs(f)> 1.0d+300)break; end;
f0=f1;
f1=f;
end;
nm=fix(k-1);
dk(0+1)=-sk(1+1);
for  k=1:nm;
dk(k+1)=-sk(k-1+1)-(k+1.0d0)./x.*sk(k+1);
end;  k=fix(nm)+1;
return;
end

