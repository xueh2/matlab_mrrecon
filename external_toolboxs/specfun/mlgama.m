function mlgama
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===================================================
%       Purpose: This program computes the gamma function
%                â(x)for x > 0 using subroutine LGAMA
%       Examples:
%                  x           â(x)
%                -------------------------
%                 0.5     .1772453851D+01
%                 2.5     .1329340388D+01
%                 5.0     .2400000000D+02
%                 7.5     .1871254306D+04
%                10.0     .3628800000D+06
%       ===================================================
x=[];gl=[];
fprintf(1,'%s \n','   x           â(x)');
fprintf(1,'%s \n',' -------------------------');
for  l=0:5:20;
x=0.5d0.*l;
if(l == 0)x=0.5; end;
[dumvar1,x,gl]=lgama(1,x,gl);
fprintf(1,[repmat(' ',1,1),'%5.1g','%20.10g' ' \n'],x,gl);
end;  l=20+1;
fprintf(1,'%s \n','please enter x:');
%        READ(*,*)X
x=10.0;
[dumvar1,x,gl]=lgama(1,x,gl);
fprintf(1,[repmat(' ',1,1),'%5.1g','%20.10g' ' \n'],x,gl);
%format(1x,f5.1,d20.10);
end
function [kf,x,gl]=lgama(kf,x,gl,varargin);
%       ==================================================
%       Purpose: Compute gamma function â(x)or ln[â(x)]
%       Input:   x  --- Argument of â(x,x > 0)
%                KF --- Function code
%                       KF=1 for â(x); KF=0 for ln[â(x)]
%       Output:  GL --- â(x)or ln[â(x)]
%       ==================================================
 a=zeros(1,10);
a(:)=[8.333333333333333d-02,-2.777777777777778d-03,7.936507936507937d-04,-5.952380952380952d-04,8.417508417508418d-04,-1.917526917526918d-03,6.410256410256410d-03,-2.955065359477124d-02,1.796443723688307d-01,-1.39243221690590d+00];
x0=x;
ifoo=1;
if(x == 1.0|x == 2.0);
gl=0.0d0;
ifoo=0;
elseif(x <= 7.0);
n=fix(7-x);
x0=x+n;
end;
if(ifoo == 1);
x2=1.0d0./(x0.*x0);
xp=6.283185307179586477d0;
gl0=a(10);
for  k=9:-1:1;
gl0=gl0.*x2+a(k);
end;  k=1-1;
gl=gl0./x0+0.5d0.*log(xp)+(x0-.5d0).*log(x0)-x0;
if(x <= 7.0);
for  k=1:n;
gl=gl-log(x0-1.0d0);
x0=x0-1.0d0;
end;  k=n+1;
end;
end;
if(kf == 1)gl=exp(gl); end;
return;
end

