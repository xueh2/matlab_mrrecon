function mclpn
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
%                Pn(z)and Pn'(z)for a complex argument using
%                subroutine CLPN
%       Input :  x --- Real part of z
%                y --- Imaginary part of z
%                n --- Degree of Pn(z), n = 0,1,...,N
%       Output:  CPN(n)--- Pn(z)
%                CPD(n)--- Pn'(z)
%       Example: z = 3.0 +2.0 i
%       n    Re[Pn(z)]Im[Pn(z)]Re[Pn'(z)]Im[Pn'(z)]
%      -----------------------------------------------------------
%       0   .100000D+01   .000000D+00   .000000D+00   .000000D+00
%       1   .300000D+01   .200000D+01   .100000D+01   .000000D+00
%       2   .700000D+01   .180000D+02   .900000D+01   .600000D+01
%       3  -.270000D+02   .112000D+03   .360000D+02   .900000D+02
%       4  -.539000D+03   .480000D+03  -.180000D+03   .790000D+03
%       5  -.461700D+04   .562000D+03  -.481500D+04   .441000D+04
%       ==========================================================
n=[];x=[];y=[];cpn=[];cpd=[];
 cpn=zeros(1,100+1);
cpd=zeros(1,100+1);
fprintf(1,'%s \n','  please enter nmax, x and y(z=x+iy)');
%        READ(*,*)N,X,Y
n=5;
x=3.0;
y=2.0;
fprintf(1,[repmat(' ',1,3),'x =','%5.1g',',  ','y =','%5.1g' ' \n'],x,y);
fprintf(1,'%0.15g \n');
[n,x,y,cpn,cpd]=clpn(n,x,y,cpn,cpd);
fprintf(1,'%s ','  n    re[pn(z)]im[pn(z)]re[pn''(z)]');fprintf(1,'%s \n', '   im[pn''(z)]');
fprintf(1,'%s ',' ---------------------------------------------');fprintf(1,'%s \n', '--------------');
for  k=0:n;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%14.6g',1,4) ' \n'],k,cpn(k+1),cpd(k+1));
end;  k=n+1;
%format(1x,i3,4d14.6);
%format(3x,',f5.1,',  ',',f5.1);
end
function [n,x,y,cpn,cpd]=clpn(n,x,y,cpn,cpd,varargin);
%       ==================================================
%       Purpose: Compute Legendre polynomials Pn(z)and
%                their derivatives Pn'(z)for a complex
%                argument
%       Input :  x --- Real part of z
%                y --- Imaginary part of z
%                n --- Degree of Pn(z), n = 0,1,2,...
%       Output:  CPN(n)--- Pn(z)
%                CPD(n)--- Pn'(z)
%       ==================================================
z=complex(x,y);
cpn(0+1)=complex(1.0d0,0.0d0);
cpn(1+1)=z;
cpd(0+1)=complex(0.0d0,0.0d0);
cpd(1+1)=complex(1.0d0,0.0d0);
cp0=complex(1.0d0,0.0d0);
cp1=z;
for  k=2:n;
cpf=(2.0d0.*k-1.0d0)./k.*z.*cp1-(k-1.0d0)./k.*cp0;
cpn(k+1)=cpf;
if(abs(x)== 1.0d0&y == 0.0d0);
cpd(k+1)=0.5d0.*x.^(k).*k.*(k+1.0d0);
else;
cpd(k+1)=k.*(cp1-z.*cpf)./(1.0d0-z.*z);
end;
cp0=cp1;
cp1=cpf;
end;  k=fix(n)+1;
return;
end

