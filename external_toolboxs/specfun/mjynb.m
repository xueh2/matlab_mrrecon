function mjynb
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
%                using subroutine JYNB
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
% READ(*,*)N,X
n=0;
x=10.0;
fprintf(1,[repmat(' ',1,3),'nmax =','%3g',',    ','x =','%6.1g' ' \n'],n,x);
if(n <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%    READ(*,*)NS
ns=1;
end;
fprintf(1,'%0.15g \n');
[n,x,nm,bj,dj,by,dy]=jynb(n,x,nm,bj,dj,by,dy);
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
function [n,x,nm,bj,dj,by,dy]=jynb(n,x,nm,bj,dj,by,dy,varargin);
%       =====================================================
%       Purpose: Compute Bessel functions Jn(x), Yn(x)and
%                their derivatives
%       Input :  x --- Argument of Jn(x)and Yn(x,x ע 0)
%                n --- Order of Jn(x)and Yn(x)
%       Output:  BJ(n)--- Jn(x)
%                DJ(n)--- Jn'(x)
%                BY(n)--- Yn(x)
%                DY(n)--- Yn'(x)
%                NM --- Highest order computed
%       Routines called:
%                MSTA1 and MSTA2 to calculate the starting
%                point for backward recurrence
%       =====================================================
  a=zeros(1,4);
b=zeros(1,4);
a1=zeros(1,4);
b1=zeros(1,4);
pi=3.141592653589793d0;
r2p=.63661977236758d0;
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
if(x <= 300.0|n > fix(0.9.*x));
if(n == 0)nm=1; end;
m=msta1(x,200);
if(m < nm);
nm=fix(m);
else;
m=msta2(x,fix(nm),15);
end;
bs=0.0d0;
su=0.0d0;
sv=0.0d0;
f2=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0)./x.*f1-f2;
if(k <= nm)bj(k+1)=f; end;
if(k == 2.*fix(k./2)&k ~= 0);
bs=bs+2.0d0.*f;
su=su+(-1).^(k./2).*f./k;
elseif(k > 1);
sv=sv+(-1).^(k./2).*k./(k.*k-1.0).*f;
end;
f2=f1;
f1=f;
end;  k=0-1;
s0=bs+f;
for  k=0:nm;
bj(k+1)=bj(k+1)./s0;
end;  k=fix(nm)+1;
ec=log(x./2.0d0)+0.5772156649015329d0;
by0=r2p.*(ec.*bj(0+1)-4.0d0.*su./s0);
by(0+1)=by0;
by1=r2p.*((ec-1.0d0).*bj(1+1)-bj(0+1)./x-4.0d0.*sv./s0);
by(1+1)=by1;
else;
a(:)=[-.7031250000000000d-01,.1121520996093750d+00,-.5725014209747314d+00,.6074042001273483d+01];
b(:)=[.7324218750000000d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02];
a1(:)=[.1171875000000000d+00,-.1441955566406250d+00,.6765925884246826d+00,-.6883914268109947d+01];
b1(:)=[-.1025390625000000d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02];
t1=x-0.25d0.*pi;
p0=1.0d0;
q0=-0.125d0./x;
for  k=1:4;
p0=p0+a(k).*x.^(-2.*k);
q0=q0+b(k).*x.^(-2.*k-1);
end;  k=4+1;
cu=sqrt(r2p./x);
bj0=cu.*(p0.*cos(t1)-q0.*sin(t1));
by0=cu.*(p0.*sin(t1)+q0.*cos(t1));
bj(0+1)=bj0;
by(0+1)=by0;
t2=x-0.75d0.*pi;
p1=1.0d0;
q1=0.375d0./x;
for  k=1:4;
p1=p1+a1(k).*x.^(-2.*k);
q1=q1+b1(k).*x.^(-2.*k-1);
end;  k=4+1;
bj1=cu.*(p1.*cos(t2)-q1.*sin(t2));
by1=cu.*(p1.*sin(t2)+q1.*cos(t2));
bj(1+1)=bj1;
by(1+1)=by1;
for  k=2:nm;
bjk=2.0d0.*(k-1.0d0)./x.*bj1-bj0;
bj(k+1)=bjk;
bj0=bj1;
bj1=bjk;
end;  k=fix(nm)+1;
end;
dj(0+1)=-bj(1+1);
for  k=1:nm;
dj(k+1)=bj(k-1+1)-k./x.*bj(k+1);
end;  k=fix(nm)+1;
for  k=2:nm;
byk=2.0d0.*(k-1.0d0).*by1./x-by0;
by(k+1)=byk;
by0=by1;
by1=byk;
end;  k=fix(nm)+1;
dy(0+1)=-by(1+1);
for  k=1:nm;
dy(k+1)=by(k-1+1)-k.*by(k+1)./x;
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

