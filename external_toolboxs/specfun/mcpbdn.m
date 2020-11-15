function mcpbdn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes parabolic cylinder functions
%       Dn(z)for an integer order and a complex argument
%       using subroutine CPBDN
%       Input :  x --- Real part of z
%       y --- Imaginary part of z
%       n --- Order of Dn(z)
%       Output:  CPB(|n|)--- Dn(z)
%       CPD(|n|)--- Dn'(z)
%       Example:
%       z = 5.0+ 5.0 i
%       n     Re[Dn(z)]Im[Dn(z)]Re[Dn'(z)]Im[Dn'(z)]
%       ----------------------------------------------------------------
%       -
%       0   .99779828D+00  .66321897D-01 -.23286910D+01 -.26603004D+01
%       1   .46573819D+01  .53206009D+01  .26558457D+01 -.24878635D+02
%       2  -.43138931D+01  .49823592D+02  .14465848D+03 -.10313305D+03
%       3  -.28000219D+03  .21690729D+03  .12293320D+04  .30720802D+03
%       4  -.24716057D+04 -.46494526D+03  .38966424D+04  .82090067D+04
%       -1   .10813809D+00 -.90921592D-01 -.50014908D+00 -.23280660D-01
%       -2   .24998820D-02 -.19760577D-01 -.52486940D-01  .47769856D-01
%       -3  -.15821033D-02 -.23090595D-02 -.68249161D-03  .10032670D-01
%       -4  -.37829961D-03 -.10158757D-03  .89032322D-03  .11093416D-02
%       =============================================================
n=[];z=[];cpb=[];cpd=[];
 cpb=zeros(1,100+1);
cpd=zeros(1,100+1);
fprintf(1,'%s \n','please enter n, x and y ');
%       READ(*,*)N,X,Y
n=5;
x=5.0;
y=5.0;
fprintf(1,[repmat(' ',1,1),'n =','%3g',',   z =x+iy :','%6.2g','+','%6.2g',' i' ' \n'],n,x,y);
z=complex(x,y);
n0=abs(n);
[n,z,cpb,cpd]=cpbdn(n,z,cpb,cpd);
fprintf(1,'%0.15g \n');
if(n >= 0);
fprintf(1,'%s ','  n     re[dn(z)]im[dn(z)]');fprintf(1,'%s \n','re[dn''(z)]im[dn''(z)]');
else;
fprintf(1,'%s ',' -n     re[dn(z)]im[dn(z)]');fprintf(1,'%s \n','re[dn''(z)]im[dn''(z)]');
end;
fprintf(1,'%s ','-------------------------------------------');fprintf(1,'%s \n','-------------------------');
for  i=0:n0;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%16.8g',1,4) ' \n'],i,cpb(i+1+1),cpd(i+1+1));
end;  i=n0+1;
%format(1x,i3,4d16.8);
%format(1x,',i3,',   z =x+iy :',f6.2,'+',f6.2,' i');
end
function [n,z,cpb,cpd]=cpbdn(n,z,cpb,cpd,varargin);
%       ==================================================
%       Purpose: Compute the parabolic cylinder functions
%       Dn(z)and Dn'(z)for a complex argument
%       Input:   z --- Complex argument of Dn(z)
%       n --- Order of Dn(z,n=0,�1,�2,���)
%       Output:  CPB(|n|)--- Dn(z)
%       CPD(|n|)--- Dn'(z)
%       Routines called:
%(1)CPDSA for computing Dn(z)for a small |z|
%(2)CPDLA for computing Dn(z)for a large |z|
%       ==================================================
z1=[];cf1=[];n0=[];cfa=[];n1=[];cfb=[];
pi=3.141592653589793d0;
x=real(z);
a0=abs(z);
c0=complex(0.0d0,0.0d0);
ca0=exp(-0.25d0.*z.*z);
if(n >= 0);
cf0=ca0;
cf1=z.*ca0;
cpb(0+1+1)=cf0;
cpb(1+1+1)=cf1;
for  k=2:n;
cf=z.*cf1-(k-1.0d0).*cf0;
cpb(k+1+1)=cf;
cf0=cf1;
cf1=cf;
end;  k=fix(n)+1;
else;
n0=-fix(n);
if(x <= 0.0|abs(z)== 0.0);
cf0=ca0;
cpb(0+1+1)=cf0;
z1=-z;
if(a0 <= 7.0);
[dumvar1,z1,cf1]=cpdsa(-1,z1,cf1);
else;
[dumvar1,z1,cf1]=cpdla(-1,z1,cf1);
end;
cf1=sqrt(2.0d0.*pi)./ca0-cf1;
cpb(1+1+1)=cf1;
for  k=2:n0;
cf=(-z.*cf1+cf0)./(k-1.0d0);
cpb(k+1+1)=cf;
cf0=cf1;
cf1=cf;
end;  k=n0+1;
else;
if(a0 < 1.0);
[dumvar1,z,cfa]=cpdsa(-n0,z,cfa);
cpb(n0+1+1)=cfa;
n1=n0+1;
[dumvar1,z,cfb]=cpdsa(-n1,z,cfb);
cpb(n1+1+1)=cfb;
nm1=n0-1;
cpb(n1+1+1)=cfb;
cpb(n0+1+1)=cfa;
for  k=nm1:-1:0;
cf=z.*cfa+(k+1.0d0).*cfb;
cpb(k+1+1)=cf;
cfb=cfa;
cfa=cf;
end;  k=0-1;
else;
m=120+abs(fix(n));
cfa=c0;
cfb=complex(1.0d-30,0.0d0);
for  k=m:-1:0;
cf=z.*cfb+(k+1.0d0).*cfa;
if(k <= n0)cpb(k+1+1)=cf; end;
cfa=cfb;
cfb=cf;
end;  k=0-1;
cs0=ca0./cf;
for  k=0:n0;
cpb(k+1+1)=cs0.*cpb(k+1+1);
end;  k=n0+1;
end;
end;
end;
cpd(0+1+1)=-0.5d0.*z.*cpb(0+1+1);
if(n >= 0);
for  k=1:n;
cpd(k+1+1)=-0.5d0.*z.*cpb(k+1+1)+k.*cpb(k-1+1+1);
end;  k=fix(n)+1;
else;
for  k=1:n0;
cpd(k+1+1)=0.5d0.*z.*cpb(k+1+1)-cpb(k-1+1+1);
end;  k=n0+1;
end;
return;
end
function [n,z,cdn]=cpdsa(n,z,cdn,varargin);
%       ===========================================================
va0=[];ga0=[];xn=[];g1=[];vt=[];g0=[];vm=[];gm=[];
eps=1.0d-15;
pi=3.141592653589793d0;
sq2=sqrt(2.0d0);
ca0=exp(-.25d0.*z.*z);
va0=0.5d0.*(1.0d0-fix(n));
if(n == 0.0);
cdn=ca0;
else;
if(abs(z)== 0.0);
if(va0 <= 0.0&va0 == fix(va0));
cdn=0.0d0;
else;
[va0,ga0]=gaih(va0,ga0);
pd=sqrt(pi)./(2.0d0.^(-.5d0.*fix(n)).*ga0);
cdn=complex(pd,0.0d0);
end;
else;
xn=-fix(n);
[xn,g1]=gaih(xn,g1);
cb0=2.0d0.^(-0.5d0.*fix(n)-1.0d0).*ca0./g1;
vt=-.5d0.*fix(n);
[vt,g0]=gaih(vt,g0);
cdn0=complex(g0,0.0d0);
cdn1=complex(0.0d0,0.0d0);
cdn2=complex(0.0d0,0.0d0);
zq=sq2.*z;
fm=1.0d0;
for  m=1:150;
vm=.5d0.*(m-fix(n));
[vm,gm]=gaih(vm,gm);
fm=1.0d0;
for  j=1:m;
fm=fm.*j;
end;  j=m+1;
cdw=(-1).^m.*gm.*zq.^m./fm;
if(m == 2.*fix(m./2))cdn1=cdn1+cdw; end;
if(m ~= 2.*fix(m./2))cdn2=cdn2+cdw; end;
if(abs(cdw)< abs(cdn1).*eps)break; end;
if(abs(cdw)< abs(cdn2).*eps)break; end;
end;
cdn=cb0.*(cdn0+cdn1+cdn2);
end;
end;
return;
end
function [n,z,cdn]=cpdla(n,z,cdn,varargin);
%       ===========================================================
%       Purpose: Compute complex parabolic cylinder function Dn(z)
%       for large argument
%       Input:   z   --- Complex argument of Dn(z)
%       n   --- Order of Dn(z,n = 0,�1,�2,���)
%       Output:  CDN --- Dn(z)
%       ===========================================================
cb0=z.^fix(n).*exp(-.25d0.*z.*z);
cr=complex(1.0d0,0.0d0);
cdn=complex(1.0d0,0.0d0);
for  k=1:16;
cr=-0.5d0.*cr.*(2.0.*k-fix(n)-1.0).*(2.0.*k-fix(n)-2.0)./(k.*z.*z);
cdn=cdn+cr;
if(abs(cr)< abs(cdn).*1.0d-12)break; end;
end;
cdn=cb0.*cdn;
return;
end
function [x,ga]=gaih(x,ga,varargin);
%       =====================================================
%       Purpose: Compute gamma function �(x)
%       Input :  x  --- Argument of �(x), x = n/2, n=1,2,���
%       Output:  GA --- �(x)
%       =====================================================
pi=3.141592653589793d0;
if(x == fix(x)&x > 0.0);
ga=1.0d0;
m1=fix(x-1.0);
for  k=2:m1;
ga=ga.*k;
end;  k=m1+1;
elseif(x+.5d0 == fix(x+.5d0)&x > 0.0);
m=fix(x);
ga=sqrt(pi);
for  k=1:m;
ga=0.5d0.*ga.*(2.0d0.*k-1.0d0);
end;  k=m+1;
end;
return;
end

