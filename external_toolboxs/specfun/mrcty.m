function mrcty
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =======================================================
%       Purpose: This program computes the Riccati-Bessel
%                functions of the second kind and their
%                derivatives using subroutine RCTY
%       Input:   x --- Argument of Riccati-Bessel function
%                n --- Order of yn(x)
%       Output:  RY(n)--- xúyn(x)
%                DY(n)---[xúyn(x)]'
%       Example: x = 10.0
%                  n        xúyn(x)[xúyn(x)]'
%                --------------------------------------------
%                  0     .8390715291D+00    -.5440211109D+00
%                  1     .6279282638D+00     .7762787027D+00
%                  2    -.6506930499D+00     .7580668738D+00
%                  3    -.9532747888D+00    -.3647106133D+00
%                  4    -.1659930220D-01    -.9466350679D+00
%                  5     .9383354168D+00    -.4857670106D+00
%       =======================================================
n=[];x=[];nm=[];ry=[];dy=[];
 ry=zeros(1,250+1);
dy=zeros(1,250+1);
fprintf(1,'%s \n','  please enter n and x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',    ','x =','%6.2g' ' \n'],n,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
fprintf(1,'%0.15g \n');
[n,x,nm,ry,dy]=rcty(n,x,nm,ry,dy);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n        xúyn(x)[xúyn(x)]''');
fprintf(1,'%s \n','--------------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%20.10g',1,2) ' \n'],k,ry(k+1),dy(k+1));
end;  k=nm+1;
%format(1x,i3,2d20.10);
%format(3x,,i3,',    ',,f6.2);
end
function [n,x,nm,ry,dy]=rcty(n,x,nm,ry,dy,varargin);
%       ========================================================
%       Purpose: Compute Riccati-Bessel functions of the second
%                kind and their derivatives
%       Input:   x --- Argument of Riccati-Bessel function
%                n --- Order of yn(x)
%       Output:  RY(n)--- xúyn(x)
%                DY(n)---[xúyn(x)]'
%                NM --- Highest order computed
%       ========================================================
nm=fix(fix(n));
if(x < 1.0d-60);
for  k=0:n;
ry(k+1)=-1.0d+300;
dy(k+1)=1.0d+300;
end;  k=fix(n)+1;
ry(0+1)=-1.0d0;
dy(0+1)=0.0d0;
return;
end;
ry(0+1)=-cos(x);
ry(1+1)=ry(0+1)./x-sin(x);
rf0=ry(0+1);
rf1=ry(1+1);
for  k=2:n;
rf2=(2.0d0.*k-1.0d0).*rf1./x-rf0;
if(abs(rf2)> 1.0d+300)break; end;
ry(k+1)=rf2;
rf0=rf1;
rf1=rf2;
end;
nm=fix(k-1);
dy(0+1)=sin(x);
for  k=1:nm;
dy(k+1)=-k.*ry(k+1)./x+ry(k-1+1);
end;  k=fix(nm)+1;
return;
end

