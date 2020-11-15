function msphy
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes the spherical Bessel
%                functions yn(x)and yn'(x)using subroutine
%                SPHY
%       Input :  x --- Argument of yn(x,x ע 0)
%                n --- Order of yn(x,n = 0,1,תתת, ף 250)
%       Output:  SY(n)--- yn(x)
%                DY(n)--- yn'(x)
%       Example:   x = 10.0
%                  n          yn(x)yn'(x)
%                --------------------------------------------
%                  0     .8390715291D-01    -.6279282638D-01
%                  1     .6279282638D-01     .7134858763D-01
%                  2    -.6506930499D-01     .8231361788D-01
%                  3    -.9532747888D-01    -.2693831344D-01
%                  4    -.1659930220D-02    -.9449751377D-01
%                  5     .9383354168D-01    -.5796005523D-01
%       ========================================================
n=[];x=[];nm=[];sy=[];dy=[];
 sy=zeros(1,250+1);
dy=zeros(1,250+1);
fprintf(1,'%s \n','please enter n and x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',     ','x=','%6.1g' ' \n'],n,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,x,nm,sy,dy]=sphy(n,x,nm,sy,dy);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n          yn(x)yn''(x)');
fprintf(1,'%s \n','--------------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%20.10g',1,2) ' \n'],k,sy(k+1),dy(k+1));
end;  k=nm+1;
%format(1x,i3,2d20.10);
%format(3x,,i3,',     ',,f6.1);
end
function [n,x,nm,sy,dy]=sphy(n,x,nm,sy,dy,varargin);
%       ======================================================
%       Purpose: Compute spherical Bessel functions yn(x)and
%                their derivatives
%       Input :  x --- Argument of yn(x,x ע 0)
%                n --- Order of yn(x,n = 0,1,תתת)
%       Output:  SY(n)--- yn(x)
%                DY(n)--- yn'(x)
%                NM --- Highest order computed
%       ======================================================
nm=fix(fix(n));
if(x < 1.0d-60);
for  k=0:n;
sy(k+1)=-1.0d+300;
dy(k+1)=1.0d+300;
end;  k=fix(n)+1;
return;
end;
sy(0+1)=-cos(x)./x;
sy(1+1)=(sy(0+1)-sin(x))./x;
f0=sy(0+1);
f1=sy(1+1);
for  k=2:n;
f=(2.0d0.*k-1.0d0).*f1./x-f0;
sy(k+1)=f;
if(abs(f)>= 1.0d+300)break; end;
f0=f1;
f1=f;
end;
nm=fix(k-1);
dy(0+1)=(sin(x)+cos(x)./x)./x;
for  k=1:nm;
dy(k+1)=sy(k-1+1)-(k+1.0d0).*sy(k+1)./x;
end;  k=fix(nm)+1;
return;
end

