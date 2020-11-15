function mcik01
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =============================================================
%       Purpose: This program computes the modified Bessel functions
%                I0(z), I1(z), K0(z), K1(z), and their derivatives
%                for a complex argument using subroutine CIK01
%       Input :  z --- Complex argument
%       Output:  CBI0 --- I0(z)
%                CDI0 --- I0'(z)
%                CBI1 --- I1(z)
%                CDI1 --- I1'(z)
%                CBK0 --- K0(z)
%                CDK0 --- K0'(z)
%                CBK1 --- K1(z)
%                CDK1 --- K1'(z)
%       Example: z = 20.0 + i 10.0
%     n      Re[In(z)]Im[In(z)]Re[In'(z)]Im[In'(z)]
%    -----------------------------------------------------------------
%     0   -.38773811D+08 -.13750292D+08 -.37852037D+08 -.13869150D+08
%     1   -.37852037D+08 -.13869150D+08 -.36982347D+08 -.13952566D+08
%     n      Re[Kn(z)]Im[Kn(z)]Re[Kn'(z)]Im[Kn'(z)]
%    -----------------------------------------------------------------
%     0   -.37692389D-09  .39171613D-09  .38056380D-09 -.40319029D-09
%     1   -.38056380D-09  .40319029D-09  .38408264D-09 -.41545502D-09
%       =============================================================
z=[];cbi0=[];cdi0=[];cbi1=[];cdi1=[];cbk0=[];cdk0=[];cbk1=[];cdk1=[];
fprintf(1,'%s \n','  please enter x and y(z=x+iy)');
%        READ(*,*)X,Y
x=20.0;
y=10.0;
z=complex(x,y);
fprintf(1,[repmat(' ',1,3),'z =','%7.2g',' + i','%7.2g' ' \n'],x,y);
[z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1]=cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[in(z)]im[in(z)]');fprintf(1,'%s \n','      re[in''(z)]im[in''(z)]');
fprintf(1,'%s ',' -------------------------------');fprintf(1,'%s \n','----------------------------------');
fprintf(1,[repmat(' ',1,3),'0',repmat(' ',1,2),repmat('%15.7g',1,4) ' \n'],cbi0,cdi0);
fprintf(1,[repmat(' ',1,3),'1',repmat(' ',1,2),repmat('%15.7g',1,4) ' \n'],cbi1,cdi1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[kn(z)]im[kn(z)]');fprintf(1,'%s \n','      re[kn''(z)]im[kn''(z)]');
fprintf(1,'%s ',' -------------------------------');fprintf(1,'%s \n','----------------------------------');
fprintf(1,[repmat(' ',1,3),'0',repmat(' ',1,2),repmat('%15.7g',1,4) ' \n'],cbk0,cdk0);
fprintf(1,[repmat(' ',1,3),'1',repmat(' ',1,2),repmat('%15.7g',1,4) ' \n'],cbk1,cdk1);
%format(3x,'0',2x,4d15.7);
%format(3x,'1',2x,4d15.7);
%format(3x,,f7.2,' + i',f7.2);
end
function [z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1]=cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1,varargin);
%       ==========================================================
%       Purpose: Compute modified Bessel functions I0(z), I1(z),
%                K0(z), K1(z), and their derivatives for a
%                complex argument
%       Input :  z --- Complex argument
%       Output:  CBI0 --- I0(z)
%                CDI0 --- I0'(z)
%                CBI1 --- I1(z)
%                CDI1 --- I1'(z)
%                CBK0 --- K0(z)
%                CDK0 --- K0'(z)
%                CBK1 --- K1(z)
%                CDK1 --- K1'(z)
%       ==========================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,10);
pi=3.141592653589793d0;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z2=z.*z;
z1=z;
if(a0 == 0.0d0);
cbi0=complex(1.0d0,0.0d0);
cbi1=complex(0.0d0,0.0d0);
cdi0=complex(0.0d0,0.0d0);
cdi1=complex(0.5d0,0.0d0);
cbk0=complex(1.0d+300,0.0d0);
cbk1=complex(1.0d+300,0.0d0);
cdk0=-complex(1.0d+300,0.0d0);
cdk1=-complex(1.0d+300,0.0d0);
return;
end;
if(real(z)< 0.0)z1=-z; end;
if(a0 <= 18.0);
cbi0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:50;
cr=0.25d0.*cr.*z2./(k.*k);
cbi0=cbi0+cr;
if(abs(cr./cbi0)< 1.0d-15)break; end;
end;
cbi1=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:50;
cr=0.25d0.*cr.*z2./(k.*(k+1));
cbi1=cbi1+cr;
if(abs(cr./cbi1)< 1.0d-15)break; end;
end;
cbi1=0.5d0.*z1.*cbi1;
else;
a(:)=[0.125d0,7.03125d-2,7.32421875d-2,1.1215209960938d-1,2.2710800170898d-1,5.7250142097473d-1,1.7277275025845d0,6.0740420012735d0,2.4380529699556d01,1.1001714026925d02,5.5133589612202d02,3.0380905109224d03];
b(:)=[-0.375d0,-1.171875d-1,-1.025390625d-1,-1.4419555664063d-1,-2.7757644653320d-1,-6.7659258842468d-1,-1.9935317337513d0,-6.8839142681099d0,-2.7248827311269d01,-1.2159789187654d02,-6.0384407670507d02,-3.3022722944809d03];
k0=12;
if(a0 >= 35.0)k0=9; end;
if(a0 >= 50.0)k0=7; end;
ca=exp(z1)./sqrt(2.0d0.*pi.*z1);
cbi0=complex(1.0d0,0.0d0);
zr=1.0d0./z1;
for  k=1:k0;
cbi0=cbi0+a(k).*zr.^k;
end;  k=k0+1;
cbi0=ca.*cbi0;
cbi1=complex(1.0d0,0.0d0);
for  k=1:k0;
cbi1=cbi1+b(k).*zr.^k;
end;  k=k0+1;
cbi1=ca.*cbi1;
end;
if(a0 <= 9.0);
cs=complex(0.0d0,0.0d0);
ct=-log(0.5d0.*z1)-0.5772156649015329d0;
w0=0.0d0;
cr=complex(1.0d0,0.0d0);
for  k=1:50;
w0=w0+1.0d0./k;
cr=0.25d0.*cr./(k.*k).*z2;
cs=cs+cr.*(w0+ct);
if(abs((cs-cw)./cs)< 1.0d-15)go to 45; end;
cw=cs;
end;  k=50+1;
cbk0=ct+cs;
else;
a1(:)=[0.125d0,0.2109375d0,1.0986328125d0,1.1775970458984d01,2.1461706161499d02,5.9511522710323d03,2.3347645606175d05,1.2312234987631d07,8.401390346421d08,7.2031420482627d10];
cb=0.5d0./z1;
zr2=1.0d0./z2;
cbk0=complex(1.0d0,0.0d0);
for  k=1:10;
cbk0=cbk0+a1(k).*zr2.^k;
end;  k=10+1;
cbk0=cb.*cbk0./cbi0;
end;
cbk1=(1.0d0./z1-cbi1.*cbk0)./cbi0;
if(real(z)< 0.0);
if(imag(z)< 0.0)cbk0=cbk0+ci.*pi.*cbi0; end;
if(imag(z)> 0.0)cbk0=cbk0-ci.*pi.*cbi0; end;
if(imag(z)< 0.0)cbk1=-cbk1+ci.*pi.*cbi1; end;
if(imag(z)> 0.0)cbk1=-cbk1-ci.*pi.*cbi1; end;
cbi1=-cbi1;
end;
cdi0=cbi1;
cdi1=cbi0-1.0d0./z.*cbi1;
cdk0=-cbk1;
cdk1=-cbk0-1.0d0./z.*cbk1;
return;
end

