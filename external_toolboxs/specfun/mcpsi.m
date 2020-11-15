function mcpsi
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     =========================================================
%     Purpose: This program computes the psi function psi(z)
%     for a complex argument using subroutine CPSI
%     Input :  x   --- Real part of z
%     y   --- Imaginary part of z
%     Output:  PSR --- Real part of psi(z)
%     PSI --- Imaginary part of psi(z)
%     Examples:
%     x       y      Re[psi(z)]Im[psi(z)]
%     -------------------------------------------
%     3.0     2.0     1.16459152      .67080728
%     3.0    -2.0     1.16459152     -.67080728
%     -3.0     2.0     1.39536075     2.62465344
%     -3.0    -2.0     1.39536075    -2.62465344
%     =========================================================
x=[];y=[];psr=[];psi=[];
  x=0;
 y=0;
 psr=0;
 psi=0;
fprintf(1,'%s \n','please enter x and y(z=x+iy)');
%     READ(*,*)X,Y
x=3.0;
y=2.0;
fprintf(1,[repmat(' ',1,1),'x=','%6.2g',repmat(' ',1,6),'y=','%6.2g' ' \n'],x,y);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   x       y      re[psi(z)]im[psi(z)]');
fprintf(1,'%s \n',' ----------------------------------------------');
[x,y,psr,psi]=cpsi(x,y,psr,psi);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,3),'%5.1g',repmat(' ',1,2),repmat('%16.8g',1,2) ' \n'],x,y,psr,psi);
%format(1x,f5.1,3x,f5.1,2x,2e16.8);
%format(1x,,f6.2,6x,,f6.2);
end
function [x,y,psr,psi]=cpsi(x,y,psr,psi,varargin);
%     =============================================
%     Purpose: Compute the psi function for a
%     complex argument
%     Input :  x   --- Real part of z
%     y   --- Imaginary part of z
%     Output:  PSR --- Real part of psi(z)
%     PSI --- Imaginary part of psi(z)
%     =============================================
 a=zeros(1,8);
a(:)=[-.8333333333333d-01,.83333333333333333d-02,-.39682539682539683d-02,.41666666666666667d-02,-.75757575757575758d-02,.21092796092796093d-01,-.83333333333333333d-01,.4432598039215686d0];
x1=0.0;
pi=3.141592653589793d0;
if(y == 0.0d0&x == fix(x)&x <= 0.0d0);
psr=1.0d+300;
psi=0.0d0;
else;
if(x < 0.0d0);
x1=x;
y1=y;
x=-x;
y=-y;
end;
x0=x;
if(x < 8.0d0);
n=fix(8-fix(x));
x0=x+n;
end;
if(x0 == 0.0d0&y ~= 0.0d0)th=0.5d0.*pi; end;
if(x0 ~= 0.0d0)th=atan(y./x0); end;
z2=x0.*x0+y.*y;
z0=sqrt(z2);
psr=log(z0)-0.5d0.*x0./z2;
psi=th+0.5d0.*y./z2;
for  k=1:8;
psr=psr+a(k).*z2.^(-k).*cos(2.0d0.*k.*th);
psi=psi-a(k).*z2.^(-k).*sin(2.0d0.*k.*th);
end;  k=8+1;
if(x < 8.0d0);
rr=0.0d0;
ri=0.0d0;
for  k=1:n;
rr=rr+(x0-k)./((x0-k).^2.0d0+y.*y);
ri=ri+y./((x0-k).^2.0d0+y.*y);
end;  k=n+1;
psr=psr-rr;
psi=psi+ri;
end;
if(x1 < 0.0d0);
tn=tan(pi.*x);
tm=tanh(pi.*y);
ct2=tn.*tn+tm.*tm;
psr=psr+x./(x.*x+y.*y)+pi.*(tn-tn.*tm.*tm)./ct2;
psi=psi-y./(x.*x+y.*y)-pi.*tm.*(1.0d0+tn.*tn)./ct2;
x=x1;
y=y1;
end;
end;
return;
end

