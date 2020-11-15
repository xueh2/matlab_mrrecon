function mlqmns
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes the associated Legendre
%                functions Qmn(x)and their derivatives Qmn'(x)
%                for a given order using subroutine LQMNS
%       Input :  x --- Argument of Qmn(x)
%                m --- Order of Qmn(x),  m = 0,1,2,...
%                n --- Degree of Qmn(x), n = 0,1,2,...
%       Output:  QM(n)--- Qmn(x)
%                QD(n)--- Qmn'(x)
%       Examples:
%                m = 1,  N = 5,  x = .5
%                n        Qmn(x)Qmn'(x)
%               -------------------------------------
%                0    .11547005D+01    .76980036D+00
%                1    .10530633D+01    .23771592D+01
%                2   -.72980606D+00    .51853281D+01
%                3   -.24918526D+01    .10914062D+01
%                4   -.19340866D+01   -.11454786D+02
%                5    .93896830D+00   -.18602587D+02
%                m = 2,  N = 5,  x = 2.5
%                n        Qmn(x)Qmn'(x)
%               -------------------------------------
%                0    .95238095D+00   -.52607710D+00
%                1    .38095238D+00   -.36281179D+00
%                2    .12485160D+00   -.17134314D+00
%                3    .36835513D-01   -.66284127D-01
%                4    .10181730D-01   -.22703958D-01
%                5    .26919481D-02   -.71662396D-02
%       =========================================================
m=[];n=[];x=[];qm=[];qd=[];
 qm=zeros(1,200+1);
qd=zeros(1,200+1);
fprintf(1,'%s \n','please enter m, n, and x ');
%        READ(*,*)M,N,X
m=1;
n=5;
x=.5;
fprintf(1,[repmat(' ',1,1),'m =','%2g',',  ','n =','%2g',',  ','x =','%5.1g' ' \n'],m,n,x);
[m,n,x,qm,qd]=lqmns(m,n,x,qm,qd);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n        qmn(x)qmn''(x)');
fprintf(1,'%s \n',' -------------------------------------');
for  j=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%17.8g',1,2) ' \n'],j,qm(j+1),qd(j+1));
end;  j=n+1;
%format(1x,i3,2d17.8);
%format(1x,',i2,',  ',',i2,',  ',',f5.1);
end
function [m,n,x,qm,qd]=lqmns(m,n,x,qm,qd,varargin);
%       ========================================================
%       Purpose: Compute associated Legendre functions Qmn(x)
%                and Qmn'(x)for a given order
%       Input :  x --- Argument of Qmn(x)
%                m --- Order of Qmn(x),  m = 0,1,2,...
%                n --- Degree of Qmn(x), n = 0,1,2,...
%       Output:  QM(n)--- Qmn(x)
%                QD(n)--- Qmn'(x)
%       ========================================================
for  k=0:n;
qm(k+1)=0.0d0;
qd(k+1)=0.0d0;
end;  k=fix(n)+1;
if(abs(x)== 1.0d0);
for  k=0:n;
qm(k+1)=1.0d+300;
qd(k+1)=1.0d+300;
end;  k=fix(n)+1;
return;
end;
ls=1;
if(abs(x)> 1.0d0)ls=-1; end;
xq=sqrt(ls.*(1.0d0-x.*x));
q0=0.5d0.*log(abs((x+1.0)./(x-1.0)));
q00=q0;
q10=-1.0d0./xq;
q01=x.*q0-1.0d0;
q11=-ls.*xq.*(q0+x./(1.0d0-x.*x));
qf0=q00;
qf1=q10;
for  k=2:m;
qm0=-2.0d0.*(k-1.0)./xq.*x.*qf1-ls.*(k-1.0).*(2.0-k).*qf0;
qf0=qf1;
qf1=qm0;
end;  k=fix(m)+1;
if(m == 0)qm0=q00; end;
if(m == 1)qm0=q10; end;
qm(0+1)=qm0;
if(abs(x)< 1.0001d0);
if(m == 0&n > 0);
qf0=q00;
qf1=q01;
for  k=2:n;
qf2=((2.0.*k-1.0d0).*x.*qf1-(k-1.0).*qf0)./k;
qm(k+1)=qf2;
qf0=qf1;
qf1=qf2;
end;  k=fix(n)+1;
end;
qg0=q01;
qg1=q11;
for  k=2:m;
qm1=-2.0d0.*(k-1.0)./xq.*x.*qg1-ls.*k.*(3.0-k).*qg0;
qg0=qg1;
qg1=qm1;
end;  k=fix(m)+1;
if(m == 0)qm1=q01; end;
if(m == 1)qm1=q11; end;
qm(1+1)=qm1;
if(m == 1&n > 1);
qh0=q10;
qh1=q11;
for  k=2:n;
qh2=((2.0.*k-1.0d0).*x.*qh1-k.*qh0)./(k-1.0);
qm(k+1)=qh2;
qh0=qh1;
qh1=qh2;
end;  k=fix(n)+1;
elseif(m >= 2);
qg0=q00;
qg1=q01;
qh0=q10;
qh1=q11;
for  l=2:n;
q0l=((2.0d0.*l-1.0d0).*x.*qg1-(l-1.0d0).*qg0)./l;
q1l=((2.0.*l-1.0d0).*x.*qh1-l.*qh0)./(l-1.0d0);
qf0=q0l;
qf1=q1l;
for  k=2:m;
qmk=-2.0d0.*(k-1.0)./xq.*x.*qf1-ls.*(k+l-1.0).*(l+2.0-k).*qf0;
qf0=qf1;
qf1=qmk;
end;  k=fix(m)+1;
qm(l+1)=qmk;
qg0=qg1;
qg1=q0l;
qh0=qh1;
qh1=q1l;
end;  l=fix(n)+1;
end;
else;
if(abs(x)> 1.1);
km=40+fix(m)+fix(n);
else;
km=(40+fix(m)+fix(n)).*fix(-1.0-1.8.*log(x-1.0));
end;
qf2=0.0d0;
qf1=1.0d0;
for  k=km:-1:0;
qf0=((2.0.*k+3.0d0).*x.*qf1-(k+2.0-fix(m)).*qf2)./(k+fix(m)+1.0);
if(k <= n)qm(k+1)=qf0; end;
qf2=qf1;
qf1=qf0;
end;  k=0-1;
for  k=0:n;
qm(k+1)=qm(k+1).*qm0./qf0;
end;  k=fix(n)+1;
end;
if(abs(x)< 1.0d0);
for  k=0:n;
qm(k+1)=(-1).^fix(m).*qm(k+1);
end;  k=fix(n)+1;
end;
qd(0+1)=((1.0d0-fix(m)).*qm(1+1)-x.*qm(0+1))./(x.*x-1.0);
for  k=1:n;
qd(k+1)=(k.*x.*qm(k+1)-(k+fix(m)).*qm(k-1+1))./(x.*x-1.0);
end;  k=fix(n)+1;
return;
end

