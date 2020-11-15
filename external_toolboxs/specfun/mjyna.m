function mjyna
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ====================================================
%       Purpose: This program computes Bessel functions
%                Jn(x)and Yn(x), and their derivatives
%                using subroutine JYNA
%       Input :  x --- Argument of Jn(x)& Yn(x,x ע 0)
%                n --- Order of Jn(x)& Yn(x)
%(n = 0,1,2,תתת, n ף 250)
%       Output:  BJ(n)--- Jn(x)
%                DJ(n)--- Jn'(x)
%                BY(n)--- Yn(x)
%                DY(n)--- Yn'(x)
%       Example:
%                x = 10.0
%                n        Jn(x)Jn'(x)
%              -------------------------------------
%                0    -.2459358D+00   -.4347275D-01
%               10     .2074861D+00    .8436958D-01
%               20     .1151337D-04    .2011954D-04
%               30     .1551096D-11    .4396479D-11
%                n        Yn(x)Yn'(x)
%              -------------------------------------
%                0     .5567117D-01   -.2490154D+00
%               10    -.3598142D+00    .1605149D+00
%               20    -.1597484D+04    .2737803D+04
%               30    -.7256142D+10    .2047617D+11
%       ====================================================
n=[];x=[];nm=[];bj=[];dj=[];by=[];dy=[];
 bj=zeros(1,250+1);
by=zeros(1,250+1);
dj=zeros(1,250+1);
dy=zeros(1,250+1);
fprintf(1,'%s \n','  please enter n, x');
%        READ(*,*)N,X
n=0;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',    ','x =','%6.1g' ' \n'],n,x);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%           READ(*,*)NS
ns=1;
end;
fprintf(1,'%0.15g \n');
[n,x,nm,bj,dj,by,dy]=jyna(n,x,nm,bj,dj,by,dy);
fprintf(1,'%s \n','  n        jn(x)jn''(x)');
fprintf(1,'%s \n','--------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,1),repmat('%16.7g',1,2) ' \n'],k,bj(k+1),dj(k+1));
end;  k=nm+1;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n        yn(x)yn''(x)');
fprintf(1,'%s \n','--------------------------------------');
for  k=0:ns:nm;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,1),repmat('%16.7g',1,2) ' \n'],k,by(k+1),dy(k+1));
end;  k=nm+1;
%format(3x,,i3,',    ',,f6.1);
%format(1x,i3,1x,2d16.7);
end
function [n,x,nm,bj,dj,by,dy]=jyna(n,x,nm,bj,dj,by,dy,varargin);
%       ==========================================================
%       Purpose: Compute Bessel functions Jn(x)& Yn(x)and
%                their derivatives
%       Input :  x --- Argument of Jn(x)& Yn(x,x ע 0)
%                n --- Order of Jn(x)& Yn(x)
%       Output:  BJ(n)--- Jn(x)
%                DJ(n)--- Jn'(x)
%                BY(n)--- Yn(x)
%                DY(n)--- Yn'(x)
%                NM --- Highest order computed
%       Routines called:
%(1)JY01B to calculate J0(x), J1(x), Y0(x)& Y1(x)
%(2)MSTA1 and MSTA2 to calculate the starting
%                point for backward recurrence
%       =========================================================
bj0=[];dj0=[];bj1=[];dj1=[];by0=[];dy0=[];by1=[];dy1=[];
nm=fix(fix(n));
if(x < 1.0d-100);
for  k=0:n;
bj(k+1)=0.0d0;
dj(k+1)=0.0d0;
by(k+1)=-1.0d+300;
dy(k+1)=1.0d+300;
end;  k=fix(n)+1;
bj(0+1)=1.0d0;
dj(1+1)=0.5d0;
return;
end;
[x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1]=jy01b(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1);
bj(0+1)=bj0;
bj(1+1)=bj1;
by(0+1)=by0;
by(1+1)=by1;
dj(0+1)=dj0;
dj(1+1)=dj1;
dy(0+1)=dy0;
dy(1+1)=dy1;
if(n <= 1)return; end;
if(n < fix(0.9.*x));
for  k=2:n;
bjk=2.0d0.*(k-1.0d0)./x.*bj1-bj0;
bj(k+1)=bjk;
bj0=bj1;
bj1=bjk;
end;  k=fix(n)+1;
else;
m=msta1(x,200);
if(m < n);
nm=fix(m);
else;
m=msta2(x,fix(n),15);
end;
f2=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0)./x.*f1-f2;
if(k <= nm)bj(k+1)=f; end;
f2=f1;
f1=f;
end;  k=0-1;
if(abs(bj0)> abs(bj1));
cs=bj0./f;
else;
cs=bj1./f2;
end;
for  k=0:nm;
bj(k+1)=cs.*bj(k+1);
end;  k=fix(nm)+1;
end;
for  k=2:nm;
dj(k+1)=bj(k-1+1)-k./x.*bj(k+1);
end;  k=fix(nm)+1;
f0=by(0+1);
f1=by(1+1);
for  k=2:nm;
f=2.0d0.*(k-1.0d0)./x.*f1-f0;
by(k+1)=f;
f0=f1;
f1=f;
end;  k=fix(nm)+1;
for  k=2:nm;
dy(k+1)=by(k-1+1)-k.*by(k+1)./x;
end;  k=fix(nm)+1;
return;
end
function [x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1]=jy01b(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1,varargin);
%       =======================================================
%       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
%                Y1(x), and their derivatives
%       Input :  x   --- Argument of Jn(x)& Yn(x,x ע 0)
%       Output:  BJ0 --- J0(x)
%                DJ0 --- J0'(x)
%                BJ1 --- J1(x)
%                DJ1 --- J1'(x)
%                BY0 --- Y0(x)
%                DY0 --- Y0'(x)
%                BY1 --- Y1(x)
%                DY1 --- Y1'(x)
%       =======================================================
pi=3.141592653589793d0;
if(x == 0.0d0);
bj0=1.0d0;
bj1=0.0d0;
dj0=0.0d0;
dj1=0.5d0;
by0=-1.0d+300;
by1=-1.0d+300;
dy0=1.0d+300;
dy1=1.0d+300;
return;
elseif(x <= 4.0d0);
t=x./4.0d0;
t2=t.*t;
bj0=((((((-.5014415d-3.*t2+.76771853d-2).*t2-.0709253492d0).*t2+.4443584263d0).*t2 -1.7777560599d0).*t2+3.9999973021d0).*t2-3.9999998721d0).*t2+1.0d0;
bj1=t.*(((((((-.1289769d-3.*t2+.22069155d-2).*t2-.0236616773d0).*t2+.1777582922d0).*t2-.8888839649d0).*t2+2.6666660544d0).*t2 -3.9999999710d0).*t2+1.9999999998d0);
by0=(((((((-.567433d-4.*t2+.859977d-3).*t2-.94855882d-2).*t2+.0772975809d0).*t2 -.4261737419d0).*t2+1.4216421221d0).*t2-2.3498519931d0).*t2+1.0766115157).*t2 +.3674669052d0;
by0=2.0d0./pi.*log(x./2.0d0).*bj0+by0;
by1=((((((((.6535773d-3.*t2-.0108175626d0).*t2+.107657606d0).*t2-.7268945577d0).*t2+3.1261399273d0).*t2-7.3980241381d0).*t2+6.8529236342d0).*t2+.3932562018d0).*t2 -.6366197726d0)./x;
by1=2.0d0./pi.*log(x./2.0d0).*bj1+by1;
else;
t=4.0d0./x;
t2=t.*t;
a0=sqrt(2.0d0./(pi.*x));
p0=((((-.9285d-5.*t2+.43506d-4).*t2-.122226d-3).*t2+.434725d-3).*t2-.4394275d-2).*t2+.999999997d0;
q0=t.*(((((.8099d-5.*t2-.35614d-4).*t2+.85844d-4).*t2-.218024d-3).*t2+.1144106d-2).*t2-.031249995d0);
ta0=x-.25d0.*pi;
bj0=a0.*(p0.*cos(ta0)-q0.*sin(ta0));
by0=a0.*(p0.*sin(ta0)+q0.*cos(ta0));
p1=((((.10632d-4.*t2-.50363d-4).*t2+.145575d-3).*t2-.559487d-3).*t2+.7323931d-2).*t2+1.000000004d0;
q1=t.*(((((-.9173d-5.*t2+.40658d-4).*t2-.99941d-4).*t2+.266891d-3).*t2-.1601836d-2).*t2+.093749994d0);
ta1=x-.75d0.*pi;
bj1=a0.*(p1.*cos(ta1)-q1.*sin(ta1));
by1=a0.*(p1.*sin(ta1)+q1.*cos(ta1));
end;
dj0=-bj1;
dj1=bj0-bj1./x;
dy0=-by1;
dy1=by0-by1./x;
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

