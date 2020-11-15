function mclqmn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes the associated Legendre
%                functions Qmn(z)and their derivatives Qmn'(z)for
%                a complex argument using subroutine CLQMN
%       Definition: Qmn(z)=(-1)**m*(1-z*z)**(m/2)*dm/dzm[Qn(z)]
%                   Q0(z)=1/2*LOG[(1+z)/(1-z)](for |z|<1)
%                   Qmn(z)=(z*z-1)**(m/2)*dm/dzm[Qn(z)]
%                   Q0(z)=1/2*LOG[(z+1)/(z-1)](for |z|>1)
%       Input :  x --- Real part of z
%                y --- Imaginary part of z
%                m --- Order of Qmn(z,m = 0,1,2,תתת)
%                n --- Degree of Qmn(z,n = 0,1,2,תתת)
%       Output:  CQM(m,n)--- Qmn(z)
%                CQD(m,n)--- Qmn'(z)
%       Examples:
%                n = 5, x = 0.5, y = 0.2
%       m     Re[Qmn(z)]Im[Qmn(z)]Re[Qmn'(z)]Im[Qmn'(z)]
%      -------------------------------------------------------------
%       0    .987156D+00   .354345D+00    .324023D+01  -.447297D+01
%       1   -.240328D+01   .436861D+01    .281158D+02   .171437D+02
%       2   -.245853D+02  -.138072D+02   -.106283D+03   .913792D+02
%       3    .102723D+03  -.651233D+02   -.362578D+03  -.429802D+03
%       4    .155510D+03   .357712D+03    .196975D+04  -.287414D+02
%       5   -.167357D+04  -.680954D+03   -.193093D+04  -.925757D+03
%                n = 5, x = 2.5, y = 1.0
%       m     Re[Qmn(z)]Im[Qmn(z)]Re[Qmn'(z)]Im[Qmn'(z)]
%      -------------------------------------------------------------
%       0   -.274023D-04  -.227141D-04    .809834D-04   .210884D-04
%       1    .165620D-03   .136108D-03   -.489095D-03  -.124400D-03
%       2   -.118481D-02  -.948832D-03    .349090D-02   .825057D-03
%       3    .982179D-02   .753264D-02   -.288271D-01  -.596384D-02
%       4   -.927915D-01  -.669521D-01    .270840D+00   .451376D-01
%       5    .985601D+00   .656737D+00   -.285567D+01  -.332533D+00
%       ============================================================
m=[];n=[];x=[];y=[];cqm=[];cqd=[];
 cqm=zeros(40+1,40+1);
cqd=zeros(40+1,40+1);
fprintf(1,'%s \n','  please enter m, n, x and y ');
%        READ(*,*)M,N,X,Y
m=5;
n=5;
x=.5;
y=.2;
fprintf(1,[repmat(' ',1,1),'m =','%2g',', ','n =','%2g',', ','x =','%4.1g',', ','y =','%4.1g' ' \n'],m,n,x,y);
[dumvar1,m,n,x,y,cqm,cqd]=clqmn(40,m,n,x,y,cqm,cqd);
fprintf(1,'%s ','   m   n   re[qmn(z)]im[qmn(z)]');fprintf(1,'%s \n','re[qmn''(z)]im[qmn''(z)]');
fprintf(1,'%s ',' -----------------------------------');fprintf(1,'%s \n','------------------------------');
for  j=0:n;
fprintf(1,[repmat(' ',1,1),repmat('%4g',1,2),repmat('%14.6g',1,2),repmat(' ',1,1),repmat('%14.6g',1,2) ' \n'],m,j,cqm(m+1,j+1),cqd(m+1,j+1));
end;  j=n+1;
%format(1x,2i4,2d14.6,1x,2d14.6);
%format(1x,',i2,', ',',i2,', ',',f4.1, ', ',',f4.1);
end
function [mm,m,n,x,y,cqm,cqd]=clqmn(mm,m,n,x,y,cqm,cqd,varargin);
%       =======================================================
%       Purpose: Compute the associated Legendre functions of
%                the second kind, Qmn(z)and Qmn'(z), for a
%                complex argument
%       Input :  x  --- Real part of z
%                y  --- Imaginary part of z
%                m  --- Order of Qmn(z,m = 0,1,2,תתת)
%                n  --- Degree of Qmn(z,n = 0,1,2,תתת)
%                mm --- Physical dimension of CQM and CQD
%       Output:  CQM(m,n)--- Qmn(z)
%                CQD(m,n)--- Qmn'(z)
%       =======================================================
z=complex(x,y);
if(abs(x)== 1.0d0&y == 0.0d0);
for  i=0:m;
for  j=0:n;
cqm(i+1,j+1)=complex(1.0d+300,0.0d0);
cqd(i+1,j+1)=complex(1.0d+300,0.0d0);
end;  j=fix(n)+1;
end;  i=fix(m)+1;
return;
end;
xc=abs(z);
if(imag(z)== 0.0d0|xc < 1.0d0)ls=1; end;
if(xc > 1.0d0)ls=-1; end;
zq=sqrt(ls.*(1.0d0-z.*z));
zs=ls.*(1.0d0-z.*z);
cq0=0.5d0.*log(ls.*(1.0d0+z)./(1.0d0-z));
if(xc < 1.0001d0);
cqm(0+1,0+1)=cq0;
cqm(0+1,1+1)=z.*cq0-1.0d0;
cqm(1+1,0+1)=-1.0d0./zq;
cqm(1+1,1+1)=-zq.*(cq0+z./(1.0d0-z.*z));
for  i=0:1;
for  j=2:n;
cqm(i+1,j+1)=((2.0d0.*j-1.0d0).*z.*cqm(i+1,j-1+1)-(j+i-1.0d0).*cqm(i+1,j-2+1))./(j-i);
end;  j=fix(n)+1;
end;  i=1+1;
for  j=0:n;
for  i=2:m;
cqm(i+1,j+1)=-2.0d0.*(i-1.0d0).*z./zq.*cqm(i-1+1,j+1)-ls.*(j+i-1.0d0).*(j-i+2.0d0).*cqm(i-2+1,j+1);
end;  i=fix(m)+1;
end;  j=fix(n)+1;
else;
if(xc > 1.1);
km=40+fix(m)+fix(n);
else;
km=(40+fix(m)+fix(n)).*fix(-1.0-1.8.*log(xc-1.0));
end;
cqf2=complex(0.0d0,0.0d0);
cqf1=complex(1.0d0,0.0d0);
for  k=km:-1:0;
cqf0=((2.*k+3.0d0).*z.*cqf1-(k+2.0d0).*cqf2)./(k+1.0d0);
if(k <= n)cqm(0+1,k+1)=cqf0; end;
cqf2=cqf1;
cqf1=cqf0;
end;  k=0-1;
for  k=0:n;
cqm(0+1,k+1)=cq0.*cqm(0+1,k+1)./cqf0;
end;  k=fix(n)+1;
cqf2=0.0d0;
cqf1=1.0d0;
for  k=km:-1:0;
cqf0=((2.*k+3.0d0).*z.*cqf1-(k+1.0d0).*cqf2)./(k+2.0d0);
if(k <= n)cqm(1+1,k+1)=cqf0; end;
cqf2=cqf1;
cqf1=cqf0;
end;  k=0-1;
cq10=-1.0d0./zq;
for  k=0:n;
cqm(1+1,k+1)=cq10.*cqm(1+1,k+1)./cqf0;
end;  k=fix(n)+1;
for  j=0:n;
cq0=cqm(0+1,j+1);
cq1=cqm(1+1,j+1);
for  i=0:m-2;
cqf=-2.0d0.*(i+1).*z./zq.*cq1+(j-i).*(j+i+1.0d0).*cq0;
cqm(i+2+1,j+1)=cqf;
cq0=cq1;
cq1=cqf;
end;  i=fix(m)-2+1;
end;  j=fix(n)+1;
end;
cqd(0+1,0+1)=ls./zs;
for  j=1:n;
cqd(0+1,j+1)=ls.*j.*(cqm(0+1,j-1+1)-z.*cqm(0+1,j+1))./zs;
end;  j=fix(n)+1;
for  j=0:n;
for  i=1:m;
cqd(i+1,j+1)=ls.*i.*z./zs.*cqm(i+1,j+1)+(i+j).*(j-i+1.0d0)./zq.*cqm(i-1+1,j+1);
end;  i=fix(m)+1;
end;  j=fix(n)+1;
return;
end

