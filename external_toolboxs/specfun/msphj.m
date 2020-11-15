function msphj
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =======================================================
%       Purpose: This program computes the spherical Bessel
%                functions jn(x)and jn'(x)using subroutine
%                SPHJ
%       Input :  x --- Argument of jn(x)
%                n --- Order of jn(x,n = 0,1,תתת,ף 250)
%       Output:  SJ(n)--- jn(x)
%                DJ(n)--- jn'(x)
%       Example:   x =10.0
%                  n          jn(x)jn(x)
%                --------------------------------------------
%                  0    -.5440211109D-01    -.7846694180D-01
%                  1     .7846694180D-01    -.7009549945D-01
%                  2     .7794219363D-01     .5508428371D-01
%                  3    -.3949584498D-01     .9374053162D-01
%                  4    -.1055892851D+00     .1329879757D-01
%                  5    -.5553451162D-01    -.7226857814D-01
%       =======================================================
n=[];x=[];nm=[];sj=[];dj=[];
 sj=zeros(1,250+1);
dj=zeros(1,250+1);
fprintf(1,'%s \n','please enter n and x ');
%        READ(*,*)N,X
n=5;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',     ','x=','%5.1g' ' \n'],n,x);
if(n <= 10);
ns=1;
else;
fprintf(1,'%s \n','please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
[n,x,nm,sj,dj]=sphj(n,x,nm,sj,dj);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n          jn(x)jn''(x)');
fprintf(1,'%s \n','--------------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%20.10g',1,2) ' \n'],k,sj(k+1),dj(k+1));
end;  k=nm+1;
%format(1x,i3,2d20.10);
%format(3x,,i3,',     ',,f5.1);
end
function [n,x,nm,sj,dj]=sphj(n,x,nm,sj,dj,varargin);
%       =======================================================
%       Purpose: Compute spherical Bessel functions jn(x)and
%                their derivatives
%       Input :  x --- Argument of jn(x)
%                n --- Order of jn(x,n = 0,1,תתת)
%       Output:  SJ(n)--- jn(x)
%                DJ(n)--- jn'(x)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%       =======================================================
nm=fix(fix(n));
if(abs(x)== 1.0d-100);
for  k=0:n;
sj(k+1)=0.0d0;
dj(k+1)=0.0d0;
end;  k=fix(n)+1;
sj(0+1)=1.0d0;
dj(1+1)=.3333333333333333d0;
return;
end;
sj(0+1)=sin(x)./x;
sj(1+1)=(sj(0+1)-cos(x))./x;
if(n >= 2);
sa=sj(0+1);
sb=sj(1+1);
m=msta1(x,200);
if(m < n);
nm=fix(m);
else;
m=msta2(x,fix(n),15);
end;
f0=0.0d0;
f1=1.0d0-100;
for  k=m:-1:0;
f=(2.0d0.*k+3.0d0).*f1./x-f0;
if(k <= nm)sj(k+1)=f; end;
f0=f1;
f1=f;
end;  k=0-1;
if(abs(sa)> abs(sb))cs=sa./f; end;
if(abs(sa)<= abs(sb))cs=sb./f0; end;
for  k=0:nm;
sj(k+1)=cs.*sj(k+1);
end;  k=fix(nm)+1;
end;
dj(0+1)=(cos(x)-sin(x)./x)./x;
for  k=1:nm;
dj(k+1)=sj(k-1+1)-(k+1.0d0).*sj(k+1)./x;
end;  k=fix(nm)+1;
return;
end
function [msta1Result]=msta1(x,mp,varargin);
%       ===================================================
%       Purpose: Determine the starting point for backward
%                recurrence such that the magnitude of
%                Jn(x)at that point is about 10^(-MP)
%       Input :  x     --- Argument of Jn(x)
%                MP    --- Value of magnitude
%       Output:  MSTA1 --- Starting point
%       ===================================================
a0=abs(x);
n0=fix(1.1.*a0)+1;
f0=envj(n0,a0)-fix(mp);
n1=n0+5;
f1=envj(n1,a0)-fix(mp);
for  it=1:20;
nn=n1-(n1-n0)./(1.0d0-f0./f1);
f=envj(nn,a0)-fix(mp);
if(abs(nn-n1)< 1)break; end;
n0=n1;
f0=f1;
n1=nn;
f1=f;
end;
msta1Result=fix(nn);
return;
end
function [msta2Result]=msta2(x,n,mp,varargin);
%       ===================================================
%       Purpose: Determine the starting point for backward
%                recurrence such that all Jn(x)has MP
%                significant digits
%       Input :  x  --- Argument of Jn(x)
%                n  --- Order of Jn(x)
%                MP --- Significant digit
%       Output:  MSTA2 --- Starting point
%       ===================================================
a0=abs(x);
hmp=0.5d0.*fix(mp);
ejn=envj(fix(n),a0);
if(ejn <= hmp);
obj=fix(mp);
n0=fix(1.1.*a0);
else;
obj=hmp+ejn;
n0=fix(n);
end;
f0=envj(n0,a0)-obj;
n1=n0+5;
f1=envj(n1,a0)-obj;
for  it=1:20;
nn=n1-(n1-n0)./(1.0d0-f0./f1);
f=envj(nn,a0)-obj;
if(abs(nn-n1)< 1)break; end;
n0=n1;
f0=f1;
n1=nn;
f1=f;
end;
msta2Result=fix(nn+10);
return;
end
function [envjResult]=envj(n,x,varargin);
envjResult=0.5d0.*log10(6.28d0.*fix(n))-fix(n).*log10(1.36d0.*x./fix(n));
return;
end

