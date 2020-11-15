function mclqn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes the Legendre polynomials
%                Qn(z)and Qn'(z)for a complex argument using
%                subroutine CLQN
%       Input :  x --- Real part of z
%                y --- Imaginary part of z
%                n --- Degree of Qn(z), n = 0,1,...
%       Output:  CQN(n)--- Qn(z)
%                CQD(n)--- Qn'(z)
%       Examples:
%       z = 0.5 + 0.5 i
%       n    Re[Qn(z)]Im[Qn(z)]Re[Qn'(z)]Im[Qn'(z)]
%      -----------------------------------------------------------
%       0   .402359D+00   .553574D+00   .800000D+00   .400000D+00
%       1  -.107561D+01   .477967D+00   .602359D+00   .115357D+01
%       2  -.136636D+01  -.725018D+00  -.242682D+01   .183390D+01
%       3   .182619D+00  -.206146D+01  -.622944D+01  -.247151D+01
%       4   .298834D+01  -.110022D+01  -.114849D+01  -.125963D+02
%       5   .353361D+01   .334847D+01   .206656D+02  -.123735D+02
%       z = 3.0 + 2.0 i
%       n    Re[Qn(z)]Im[Qn(z)]Re[Qn'(z)]Im[Qn'(z)]
%      -----------------------------------------------------------
%       0   .229073D+00  -.160875D+00  -.250000D-01   .750000D-01
%       1   .896860D-02  -.244805D-01   .407268D-02   .141247D-01
%       2  -.736230D-03  -.281865D-02   .190581D-02   .155860D-02
%       3  -.264727D-03  -.227023D-03   .391535D-03   .314880D-04
%       4  -.430648D-04  -.443187D-05   .527190D-04  -.305592D-04
%       5  -.481362D-05   .265297D-05   .395108D-05  -.839883D-05
%       ==========================================================
n=[];x=[];y=[];cqn=[];cqd=[];
 cqn=zeros(1,100+1);
cqd=zeros(1,100+1);
fprintf(1,'%s \n','  please enter nmax, x and y(z=x+iy)');
%        READ(*,*)N,X,Y
n=5;
x=.5;
y=.5;
fprintf(1,[repmat(' ',1,3),'x =','%5.1g',',  ','y =','%5.1g' ' \n'],x,y);
fprintf(1,'%0.15g \n');
[n,x,y,cqn,cqd]=clqn(n,x,y,cqn,cqd);
fprintf(1,'%s ','  n    re[qn(z)]im[qn(z)]re[qn''(z)]');fprintf(1,'%s \n', '    im[qn''(z)]');
fprintf(1,'%s ',' ---------------------------------------------');fprintf(1,'%s \n', '--------------');
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%14.6g',1,4) ' \n'],k,cqn(k+1),cqd(k+1));
end;  k=n+1;
%format(1x,i3,4d14.6);
%format(3x,',f5.1,',  ',',f5.1);
end
function [n,x,y,cqn,cqd]=clqn(n,x,y,cqn,cqd,varargin);
%       ==================================================
%       Purpose: Compute the Legendre functions Qn(z)and
%                their derivatives Qn'(z)for a complex
%                argument
%       Input :  x --- Real part of z
%                y --- Imaginary part of z
%                n --- Degree of Qn(z), n = 0,1,2,...
%       Output:  CQN(n)--- Qn(z)
%                CQD(n)--- Qn'(z)
%       ==================================================
z=complex(x,y);
if(z == 1.0d0);
for  k=0:n;
cqn(k+1)=complex(1.0d+300,0.0d0);
cqd(k+1)=complex(1.0d+300,0.0d0);
end;  k=fix(n)+1;
return;
end;
ls=1;
if(abs(z)> 1.0d0)ls=-1; end;
cq0=0.5d0.*log(ls.*(1.0d0+z)./(1.0d0-z));
cq1=z.*cq0-1.0d0;
cqn(0+1)=cq0;
cqn(1+1)=cq1;
if(abs(z)< 1.0001d0);
cqf0=cq0;
cqf1=cq1;
for  k=2:n;
cqf2=((2.0d0.*k-1.0d0).*z.*cqf1-(k-1.0d0).*cqf0)./k;
cqn(k+1)=cqf2;
cqf0=cqf1;
cqf1=cqf2;
end;  k=fix(n)+1;
else;
if(abs(z)> 1.1d0);
km=40+fix(n);
else;
km=(40+fix(n)).*fix(-1.0-1.8.*log(abs(z-1.0)));
end;
cqf2=0.0d0;
cqf1=1.0d0;
for  k=km:-1:0;
cqf0=((2.*k+3.0d0).*z.*cqf1-(k+2.0d0).*cqf2)./(k+1.0d0);
if(k <= n)cqn(k+1)=cqf0; end;
cqf2=cqf1;
cqf1=cqf0;
end;  k=0-1;
for  k=0:n;
cqn(k+1)=cqn(k+1).*cq0./cqf0;
end;  k=fix(n)+1;
end;
cqd(0+1)=(cqn(1+1)-z.*cqn(0+1))./(z.*z-1.0d0);
for  k=1:n;
cqd(k+1)=(k.*z.*cqn(k+1)-k.*cqn(k-1+1))./(z.*z-1.0d0);
end;  k=fix(n)+1;
return;
end

