function mpsi
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==================================================
%       Purpose: This program computes the psi function
%                using subroutine PSI
%       Input :  x  --- Argument of psi(x)
%       Output:  PS --- psi(x)
%       Examples:
%                   x          Psi(x)
%                ------------------------
%                  .25      -4.227453533
%                  .50      -1.963510026
%                  .75      -1.085860880
%                 1.00       -.577215665
%                 1.25       -.227453533
%                 1.50        .036489974
%                 1.75        .247472454
%                 2.00        .422784335
%       ==================================================
x=[];ps=[];
  x=0;
 ps=0;
fprintf(1,'%s \n','please enter x');
%        READ(*,*)X
x=2.0;
fprintf(1,'%s \n','    x          psi(x)');
fprintf(1,'%s \n',' ------------------------');
[x,ps]=psi(x,ps);
fprintf(1,[repmat(' ',1,1),'%6.2g','%18.9g' ' \n'],x,ps);
%format(1x,f6.2,f18.9);
end
function [x,ps]=psi(x,ps,varargin);
%       ======================================
%       Purpose: Compute the psi function
%       Input :  x  --- Argument of psi(x)
%       Output:  PS --- psi(x)
%       ======================================
xa=abs(x);
pi=3.141592653589793d0;
el=.5772156649015329d0;
s=0.0d0;
if(x == fix(x)&x <= 0.0);
ps=1.0d+300;
return;
elseif(xa == fix(xa));
n=xa;
for  k=1 :n-1;
s=s+1.0d0./k;
end;  k=n-1+1;
ps=-el+s;
elseif(xa+.5 == fix(xa+.5));
n=xa-.5;
for  k=1:n;
s=s+1.0./(2.0d0.*k-1.0d0);
end;  k=n+1;
ps=-el+2.0d0.*s-1.386294361119891d0;
else;
if(xa < 10.0);
n=10-fix(xa);
for  k=0:n-1;
s=s+1.0d0./(xa+k);
end;  k=n-1+1;
xa=xa+n;
end;
x2=1.0d0./(xa.*xa);
a1=-.8333333333333d-01;
a2=.83333333333333333d-02;
a3=-.39682539682539683d-02;
a4=.41666666666666667d-02;
a5=-.75757575757575758d-02;
a6=.21092796092796093d-01;
a7=-.83333333333333333d-01;
a8=.4432598039215686d0;
ps=log(xa)-.5d0./xa+x2.*(((((((a8.*x2+a7).*x2+a6).*x2+a5).*x2+a4).*x2+a3).*x2+a2).*x2+a1);
ps=ps-s;
end;
if(x < 0.0)ps=ps-pi.*cos(pi.*x)./sin(pi.*x)-1.0d0./x; end;
return;
end

