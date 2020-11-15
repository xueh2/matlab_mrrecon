function mlqmn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===============================================================
%       Purpose: This program computes the associated Legendre
%                functions Qmn(x)and their derivatives Qmn'(x)using
%                subroutine LQMN
%       Input :  x --- Argument of Qmn(x)
%                m --- Order of Qmn(x,m = 0,1,2,תתת)
%                n --- Degree of Qmn(x,n = 0,1,2,תתת)
%       Output:  QM(m,n)--- Qmn(x)
%                QD(m,n)--- Qmn'(x)
%       Examples:
%       Qmn(x):  x = 0.5
%       n\m      0           1           2           3           4
%       ---------------------------------------------------------------
%        0     .549306   -1.154701    1.333333   -5.388603   26.666667
%        1    -.725347   -1.053063    2.666667   -6.158403   32.000000
%        2    -.818663     .729806    4.069272  -12.316806   42.666667
%        3    -.198655    2.491853    -.493486  -23.778868   85.333333
%        4     .440175    1.934087  -11.036781   -9.325204  186.818394
%       Qmn'(x): x = 0.5
%       n\m      0           1           2           3           4
%       ---------------------------------------------------------------
%        0    1.333333    -.769800    4.444444  -20.014809  145.777778
%        1    1.215973   -2.377159    3.555556  -24.633611  156.444444
%        2    -.842707   -5.185328    8.796526  -24.633611  199.111111
%        3   -2.877344   -1.091406   28.115454  -50.976710  227.555556
%        4   -2.233291   11.454786   25.483527 -197.068892  412.039838
%       Qmn(x): x = 2.0
%       n\m      0           1           2           3           4
%       ---------------------------------------------------------------
%        0     .549306    -.577350    1.333333   -5.003702   26.666667
%        1     .098612    -.203274     .666667   -3.079201   18.666667
%        2     .021184    -.064946     .277089   -1.539601   10.666667
%        3     .004871    -.019817     .104220    -.679543    5.333333
%        4     .001161    -.005887     .036816    -.276005    2.427640
%       Qmn'(x): x = 2.0
%       n\m      0           1           2           3           4
%       ---------------------------------------------------------------
%        0    -.333333     .384900   -1.111111    5.388603  -36.444444
%        1    -.117361     .249384    -.888889    4.618802  -32.000000
%        2    -.037496     .116680    -.519437    3.079201  -23.111111
%        3    -.011442     .046960    -.253375    1.720114  -14.222222
%        4    -.003399     .017331    -.110263     .849589   -7.748516
%       ===============================================================
m=[];n=[];x=[];qm=[];qd=[];
 qm=zeros(100+1,100+1);
qd=zeros(100+1,100+1);
fprintf(1,'%s \n','  please enter m, n and x');
%        READ(*,*)M,N,X
m=0;
n=4;
x=.5;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  m     n      x          qmn(x)qmn''(x)');
fprintf(1,'%s \n',' ---------------------------------------------------');
[dumvar1,m,n,x,qm,qd]=lqmn(100,m,n,x,qm,qd);
for  j=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,3),'%3g',repmat(' ',1,3),'%5.1g',repmat('%17.8g',1,2) ' \n'],m,j,x,qm(m+1,j+1),qd(m+1,j+1));
end;  j=n+1;
%format(1x,i3,3x,i3,3x,f5.1,2d17.8);
end
function [mm,m,n,x,qm,qd]=lqmn(mm,m,n,x,qm,qd,varargin);
%       ==========================================================
%       Purpose: Compute the associated Legendre functions of the
%                second kind, Qmn(x)and Qmn'(x)
%       Input :  x  --- Argument of Qmn(x)
%                m  --- Order of Qmn(x,m = 0,1,2,תתת)
%                n  --- Degree of Qmn(x,n = 0,1,2,תתת)
%                mm --- Physical dimension of QM and QD
%       Output:  QM(m,n)--- Qmn(x)
%                QD(m,n)--- Qmn'(x)
%       ==========================================================
if(abs(x)== 1.0d0);
for  i=0:m;
for  j=0:n;
qm(i+1,j+1)=1.0d+300;
qd(i+1,j+1)=1.0d+300;
end;  j=fix(n)+1;
end;  i=fix(m)+1;
return;
end;
ls=1;
if(abs(x)> 1.0d0)ls=-1; end;
xs=ls.*(1.0d0-x.*x);
xq=sqrt(xs);
q0=0.5d0.*log(abs((x+1.0d0)./(x-1.0d0)));
if(abs(x)< 1.0001d0);
qm(0+1,0+1)=q0;
qm(0+1,1+1)=x.*q0-1.0d0;
qm(1+1,0+1)=-1.0d0./xq;
qm(1+1,1+1)=-xq.*(q0+x./(1.0d0-x.*x));
for  i=0:1;
for  j=2:n;
qm(i+1,j+1)=((2.0d0.*j-1.0d0).*x.*qm(i+1,j-1+1)-(j+i-1.0d0).*qm(i+1,j-2+1))./(j-i);
end;  j=fix(n)+1;
end;  i=1+1;
for  j=0:n;
for  i=2:m;
qm(i+1,j+1)=-2.0d0.*(i-1.0d0).*x./xq.*qm(i-1+1,j+1)-ls.*(j+i-1.0d0).*(j-i+2.0d0).*qm(i-2+1,j+1);
end;  i=fix(m)+1;
end;  j=fix(n)+1;
else;
if(abs(x)> 1.1d0);
km=40+fix(m)+fix(n);
else;
km=(40+fix(m)+fix(n)).*fix(-1.0-1.8.*log(x-1.0));
end;
qf2=0.0d0;
qf1=1.0d0;
for  k=km:-1:0;
qf0=((2.*k+3.0d0).*x.*qf1-(k+2.0d0).*qf2)./(k+1.0d0);
if(k <= n)qm(0+1,k+1)=qf0; end;
qf2=qf1;
qf1=qf0;
end;  k=0-1;
for  k=0:n;
qm(0+1,k+1)=q0.*qm(0+1,k+1)./qf0;
end;  k=fix(n)+1;
qf2=0.0d0;
qf1=1.0d0;
for  k=km:-1:0;
qf0=((2.*k+3.0d0).*x.*qf1-(k+1.0d0).*qf2)./(k+2.0d0);
if(k <= n)qm(1+1,k+1)=qf0; end;
qf2=qf1;
qf1=qf0;
end;  k=0-1;
q10=-1.0d0./xq;
for  k=0:n;
qm(1+1,k+1)=q10.*qm(1+1,k+1)./qf0;
end;  k=fix(n)+1;
for  j=0:n;
q0=qm(0+1,j+1);
q1=qm(1+1,j+1);
for  i=0:m-2;
qf=-2.0d0.*(i+1).*x./xq.*q1+(j-i).*(j+i+1.0d0).*q0;
qm(i+2+1,j+1)=qf;
q0=q1;
q1=qf;
end;  i=fix(m)-2+1;
end;  j=fix(n)+1;
end;
qd(0+1,0+1)=ls./xs;
for  j=1:n;
qd(0+1,j+1)=ls.*j.*(qm(0+1,j-1+1)-x.*qm(0+1,j+1))./xs;
end;  j=fix(n)+1;
for  j=0:n;
for  i=1:m;
qd(i+1,j+1)=ls.*i.*x./xs.*qm(i+1,j+1)+(i+j).*(j-i+1.0d0)./xq.*qm(i-1+1,j+1);
end;  i=fix(m)+1;
end;  j=fix(n)+1;
return;
end

