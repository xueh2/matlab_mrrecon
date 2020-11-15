function mlamv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =======================================================
%       Purpose: This program computes the lambda functions
%                for an arbitrary order, and their derivative
%                using subroutine LAMV
%       Input :  x --- Argument of lambda function
%                v --- Order of lambda function
%(v = n+v0, 0 ó n ó 250, 0 ó v0 < 1)
%       Output:  VL(n)--- Lambda function of order n+v0
%                DL(n)--- Derivative of lambda function
%       Example: x = 10.0
%                   v         Lambda(x)Lambda'(x)
%                 ------------------------------------------
%                  0.25    -.12510515D+00    -.78558916D-01
%                  0.50    -.54402111D-01    -.78466942D-01
%                  0.75    -.13657787D-01    -.66234027D-01
%                  1.00     .86945492D-02    -.50926063D-01
%                  1.25     .19639729D-01    -.36186221D-01
%                  1.50     .23540083D-01    -.23382658D-01
%                  1.75     .23181910D-01    -.12893894D-01
%                  2.00     .20370425D-01    -.46703503D-02
%                  2.25     .16283799D-01     .15101684D-02
%                  2.50     .11691329D-01     .59243767D-02
%       =======================================================
v=[];x=[];vm=[];vl=[];dl=[];
 vl=zeros(1,250+1);
dl=zeros(1,250+1);
fprintf(1,'%s \n','  please enter v and x ');
% READ(*,*)V,X
v=2.5;
x=10.0;
fprintf(1,[repmat(' ',1,1),'v =','%6.2g','    ','x =','%8.2g' ' \n'],v,x);
if(v <= 8);
ns=1;
else;
fprintf(1,'%s \n','  please enter order step ns');
%    READ(*,*)NS
ns=1;
end;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','   v         lambda(x)lambda''(x)');
fprintf(1,'%s \n','-------------------------------------------');
[v,x,vm,vl,dl]=lamv(v,x,vm,vl,dl);
nm=fix(vm);
v0=vm-nm;
for  k=0:ns:nm;
vk=k+v0;
fprintf(1,[repmat(' ',1,1),'%6.2g',repmat('%18.8g',1,2) ' \n'],vk,vl(k+1),dl(k+1));
end;  k=nm+1;
%format(1x,f6.2,2d18.8);
%format(1x,',f6.2,'    ',',f8.2);
end
function [v,x,vm,vl,dl]=lamv(v,x,vm,vl,dl,varargin);
%       =========================================================
%       Purpose: Compute lambda function with arbitrary order v,
%                and their derivative
%       Input :  x --- Argument of lambda function
%                v --- Order of lambda function
%       Output:  VL(n)--- Lambda function of order n+v0
%                DL(n)--- Derivative of lambda function
%                VM --- Highest order computed
%       Routines called:
%(1)MSTA1 and MSTA2 for computing the starting
%                point for backward recurrence
%(2)GAM0 for computing gamma function(|x| ó 1)
%       =========================================================
v0=[];ga=[];
pi=3.141592653589793d0;
rp2=0.63661977236758d0;
x=abs(x);
x2=x.*x;
n=fix(v);
v0=v-n;
vm=v;
if(x <= 12.0d0);
for  k=0:n;
vk=v0+k;
bk=1.0d0;
r=1.0d0;
for  i=1:50;
r=-0.25d0.*r.*x2./(i.*(i+vk));
bk=bk+r;
if(abs(r)< abs(bk).*1.0d-15)break; end;
end;
vl(k+1)=bk;
uk=1.0d0;
r=1.0d0;
for  i=1:50;
r=-0.25d0.*r.*x2./(i.*(i+vk+1.0d0));
uk=uk+r;
if(abs(r)< abs(uk).*1.0d-15)break; end;
end;
dl(k+1)=-0.5d0.*x./(vk+1.0d0).*uk;
end;
return;
end;
k0=11;
if(x >= 35.0d0)k0=10; end;
if(x >= 50.0d0)k0=8; end;
for  j=0:1;
vv=4.0d0.*(j+v0).*(j+v0);
px=1.0d0;
rp=1.0d0;
for  k=1:k0;
rp=-0.78125d-2.*rp.*(vv-(4.0.*k-3.0).^2.0).*(vv-(4.0.*k-1.0).^2.0)./(k.*(2.0.*k-1.0).*x2);
px=px+rp;
end;  k=k0+1;
qx=1.0d0;
rq=1.0d0;
for  k=1:k0;
rq=-0.78125d-2.*rq.*(vv-(4.0.*k-1.0).^2.0).*(vv-(4.0.*k+1.0).^2.0)./(k.*(2.0.*k+1.0).*x2);
qx=qx+rq;
end;  k=k0+1;
qx=0.125d0.*(vv-1.0d0).*qx./x;
xk=x-(0.5d0.*(j+v0)+0.25d0).*pi;
a0=sqrt(rp2./x);
ck=cos(xk);
sk=sin(xk);
if(j == 0)bjv0=a0.*(px.*ck-qx.*sk); end;
if(j == 1)bjv1=a0.*(px.*ck-qx.*sk); end;
end;  j=1+1;
if(v0 == 0.0d0);
ga=1.0d0;
else;
[v0,ga]=gam0(v0,ga);
ga=v0.*ga;
end;
fac=(2.0d0./x).^v0.*ga;
vl(0+1)=bjv0;
dl(0+1)=-bjv1+v0./x.*bjv0;
vl(1+1)=bjv1;
dl(1+1)=bjv0-(1.0d0+v0)./x.*bjv1;
r0=2.0d0.*(1.0d0+v0)./x;
if(n <= 1);
vl(0+1)=fac.*vl(0+1);
dl(0+1)=fac.*dl(0+1)-v0./x.*vl(0+1);
vl(1+1)=fac.*r0.*vl(1+1);
dl(1+1)=fac.*r0.*dl(1+1)-(1.0d0+v0)./x.*vl(1+1);
return;
end;
if(n >= 2&n <= fix(0.9.*x));
f0=bjv0;
f1=bjv1;
for  k=2:n;
f=2.0d0.*(k+v0-1.0d0)./x.*f1-f0;
f0=f1;
f1=f;
vl(k+1)=f;
end;  k=n+1;
elseif(n >= 2);
m=msta1(x,200);
if(m < n);
n=m;
else;
m=msta2(x,n,15);
end;
f2=0.0d0;
f1=1.0d-100;
for  k=m:-1:0;
f=2.0d0.*(v0+k+1.0d0)./x.*f1-f2;
if(k <= n)vl(k+1)=f; end;
f2=f1;
f1=f;
end;  k=0-1;
if(abs(bjv0)> abs(bjv1))cs=bjv0./f; end;
if(abs(bjv0)<= abs(bjv1))cs=bjv1./f2; end;
for  k=0:n;
vl(k+1)=cs.*vl(k+1);
end;  k=n+1;
end;
vl(0+1)=fac.*vl(0+1);
for  j=1:n;
rc=fac.*r0;
vl(j+1)=rc.*vl(j+1);
dl(j-1+1)=-0.5d0.*x./(j+v0).*vl(j+1);
r0=2.0d0.*(j+v0+1)./x.*r0;
end;  j=n+1;
dl(n+1)=2.0d0.*(v0+n).*(vl(n-1+1)-vl(n+1))./x;
vm=n+v0;
return;
end
function [x,ga]=gam0(x,ga,varargin);
%       ================================================
%       Purpose: Compute gamma function â(x)
%       Input :  x  --- Argument of â(x,|x| ó 1)
%       Output:  GA --- â(x)
%       ================================================
 g=zeros(1,25);
g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14];
gr=(25);
for  k=24:-1:1;
gr=gr.*x+g(k);
end;  k=1-1;
ga=1.0d0./(gr.*x);
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

