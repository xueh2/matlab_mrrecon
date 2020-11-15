function mothpl
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes orthogonal polynomials:
%                Tn(x)or Un(x)or Ln(x)or Hn(x), and their
%                derivatives using subroutine OTHPL
%       Input :  KF --- Function code
%                       KF=1 for Chebyshev polynomial Tn(x)
%                       KF=2 for Chebyshev polynomial Un(x)
%                       KF=3 for Laguerre polynomial Ln(x)
%                       KF=4 for Hermite polynomial Hn(x)
%                n ---  Order of orthogonal polynomials
%                x ---  Argument
%       Output:  PL(n)--- Tn(x)or Un(x)or Ln(x)or Hn(x)
%                DPL(n)--- Tn'(x)or Un'(x)or Ln'(x)or Hn'(x)
%                          n = 0,1,2,...,N(N ó 100)
%       =========================================================
kf=[];n=[];x=[];pl=[];dpl=[];
 pl=zeros(1,100+1);
dpl=zeros(1,100+1);
fprintf(1,'%s \n','kf,n,x = ?');
%        READ(*,*)KF,N,X
kf=4;
n=5;
x=1.2;
fprintf(1,[repmat(' ',1,1),'kf=','%3g',repmat(' ',1,5),'nmax=','%3g',repmat(' ',1,5),'x=','%6.3g' ' \n'],kf,n,x);
fprintf(1,'%0.15g \n');
[kf,n,x,pl,dpl]=othpl(kf,n,x,pl,dpl);
if(kf == 1)fprintf(1,'%s \n','  n          tn(x)tn''(x)'); end;
if(kf == 2)fprintf(1,'%s \n','  n          un(x)un''(x)'); end;
if(kf == 3)fprintf(1,'%s \n','  n          ln(x)ln''(x)'); end;
if(kf == 4)fprintf(1,'%s \n','  n          hn(x)hn''(x)'); end;
fprintf(1,'%s \n','-----------------------------------------');
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%18.8g',1,2) ' \n'],k,pl(k+1),dpl(k+1));
end;  k=n+1;
%format(1x,,i3,5x,,i3,5x,,f6.3);
%format(1x,i3,2d18.8);
end
function [kf,n,x,pl,dpl]=othpl(kf,n,x,pl,dpl,varargin);
%       ==========================================================
%       Purpose: Compute orthogonal polynomials: Tn(x)or Un(x),
%                or Ln(x)or Hn(x), and their derivatives
%       Input :  KF --- Function code
%                       KF=1 for Chebyshev polynomial Tn(x)
%                       KF=2 for Chebyshev polynomial Un(x)
%                       KF=3 for Laguerre polynomial Ln(x)
%                       KF=4 for Hermite polynomial Hn(x)
%                n ---  Order of orthogonal polynomials
%                x ---  Argument of orthogonal polynomials
%       Output:  PL(n)--- Tn(x)or Un(x)or Ln(x)or Hn(x)
%                DPL(n)--- Tn'(x)or Un'(x)or Ln'(x)or Hn'(x)
%       =========================================================
a=2.0d0;
b=0.0d0;
c=1.0d0;
y0=1.0d0;
y1=2.0d0.*x;
dy0=0.0d0;
dy1=2.0d0;
pl(0+1)=1.0d0;
pl(1+1)=2.0d0.*x;
dpl(0+1)=0.0d0;
dpl(1+1)=2.0d0;
if(kf == 1);
y1=x;
dy1=1.0d0;
pl(1+1)=x;
dpl(1+1)=1.0d0;
elseif(kf == 3);
y1=1.0d0-x;
dy1=-1.0d0;
pl(1+1)=1.0d0-x;
dpl(1+1)=-1.0d0;
end;
for  k=2:n;
if(kf == 3);
a=-1.0d0./k;
b=2.0d0+a;
c=1.0d0+a;
elseif(kf == 4);
c=2.0d0.*(k-1.0d0);
end;
yn=(a.*x+b).*y1-c.*y0;
dyn=a.*y1+(a.*x+b).*dy1-c.*dy0;
pl(k+1)=yn;
dpl(k+1)=dyn;
y0=y1;
y1=yn;
dy0=dy1;
dy1=dyn;
end;  k=fix(n)+1;
return;
end

