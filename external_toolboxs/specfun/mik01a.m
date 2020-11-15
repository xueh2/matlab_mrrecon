function mik01a
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
%                I0(x), I1(x), K0(x), K1(x), and their derivatives
%                using subroutine IK01A
%       Input :  x   --- Argument(x ò 0)
%       Output:  BI0 --- I0(x)
%                DI0 --- I0'(x)
%                BI1 --- I1(x)
%                DI1 --- I1'(x)
%                BK0 --- K0(x)
%                DK0 --- K0'(x)
%                BK1 --- K1(x)
%                DK1 --- K1'(x)
%       Example:
%         x      I0(x)I0'(x)I1(x)I1'(x)
%       -------------------------------------------------------------
%        1.0  .1266066D+01  .5651591D+00  .5651591D+00  .7009068D+00
%       10.0  .2815717D+04  .2670988D+04  .2670988D+04  .2548618D+04
%       20.0  .4355828D+08  .4245497D+08  .4245497D+08  .4143553D+08
%       30.0  .7816723D+12  .7685320D+12  .7685320D+12  .7560546D+12
%       40.0  .1489477D+17  .1470740D+17  .1470740D+17  .1452709D+17
%       50.0  .2932554D+21  .2903079D+21  .2903079D+21  .2874492D+21
%         x      K0(x)K0'(x)K1(x)K1'(x)
%       -------------------------------------------------------------
%        1.0  .4210244D+00 -.6019072D+00  .6019072D+00 -.1022932D+01
%       10.0  .1778006D-04 -.1864877D-04  .1864877D-04 -.1964494D-04
%       20.0  .5741238D-09 -.5883058D-09  .5883058D-09 -.6035391D-09
%       30.0  .2132477D-13 -.2167732D-13  .2167732D-13 -.2204735D-13
%       40.0  .8392861D-18 -.8497132D-18  .8497132D-18 -.8605289D-18
%       50.0  .3410168D-22 -.3444102D-22  .3444102D-22 -.3479050D-22
%       =============================================================
x=[];bi0=[];di0=[];bi1=[];di1=[];bk0=[];dk0=[];bk1=[];dk1=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=1.0;
fprintf(1,[repmat(' ',1,3),'x =','%5.1g' ' \n'],x);
fprintf(1,'%s ','  x       i0(x)i0''(x)i1(x)');fprintf(1,'%s \n','          i1''(x)');
fprintf(1,'%s ','-------------------------------------------');fprintf(1,'%s \n','----------------------');
[x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1]=ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1);
fprintf(1,[repmat(' ',1,1),'%4.1g',repmat('%15.7g',1,4) ' \n'],x,bi0,di0,bi1,di1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  x       k0(x)k0''(x)k1(x)');fprintf(1,'%s \n','          k1''(x)');
fprintf(1,'%s ','-------------------------------------------');fprintf(1,'%s \n','----------------------');
fprintf(1,[repmat(' ',1,1),'%4.1g',repmat('%15.7g',1,4) ' \n'],x,bk0,dk0,bk1,dk1);
%format(',f5.1);
%format(1x,f4.1,4d15.7);
end
function [x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1]=ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1,varargin);
%       =========================================================
%       Purpose: Compute modified Bessel functions I0(x), I1(1),
%                K0(x)and K1(x), and their derivatives
%       Input :  x   --- Argument(x ò 0)
%       Output:  BI0 --- I0(x)
%                DI0 --- I0'(x)
%                BI1 --- I1(x)
%                DI1 --- I1'(x)
%                BK0 --- K0(x)
%                DK0 --- K0'(x)
%                BK1 --- K1(x)
%                DK1 --- K1'(x)
%       =========================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,8);
ww=0.0;
pi=3.141592653589793d0;
el=0.5772156649015329d0;
x2=x.*x;
if(x == 0.0d0);
bi0=1.0d0;
bi1=0.0d0;
bk0=1.0d+300;
bk1=1.0d+300;
di0=0.0d0;
di1=0.5d0;
dk0=-1.0d+300;
dk1=-1.0d+300;
return;
elseif(x <= 18.0d0);
bi0=1.0d0;
r=1.0d0;
for  k=1:50;
r=0.25d0.*r.*x2./(k.*k);
bi0=bi0+r;
if(abs(r./bi0)< 1.0d-15)break; end;
end;
bi1=1.0d0;
r=1.0d0;
for  k=1:50;
r=0.25d0.*r.*x2./(k.*(k+1));
bi1=bi1+r;
if(abs(r./bi1)< 1.0d-15)break; end;
end;
bi1=0.5d0.*x.*bi1;
else;
a(:)=[0.125d0,7.03125d-2,7.32421875d-2,1.1215209960938d-1,2.2710800170898d-1,5.7250142097473d-1,1.7277275025845d0,6.0740420012735d0,2.4380529699556d01,1.1001714026925d02,5.5133589612202d02,3.0380905109224d03];
b(:)=[-0.375d0,-1.171875d-1,-1.025390625d-1,-1.4419555664063d-1,-2.7757644653320d-1,-6.7659258842468d-1,-1.9935317337513d0,-6.8839142681099d0,-2.7248827311269d01,-1.2159789187654d02,-6.0384407670507d02,-3.3022722944809d03];
k0=12;
if(x >= 35.0)k0=9; end;
if(x >= 50.0)k0=7; end;
ca=exp(x)./sqrt(2.0d0.*pi.*x);
bi0=1.0d0;
xr=1.0d0./x;
for  k=1:k0;
bi0=bi0+a(k).*xr.^k;
end;  k=k0+1;
bi0=ca.*bi0;
bi1=1.0d0;
for  k=1:k0;
bi1=bi1+b(k).*xr.^k;
end;  k=k0+1;
bi1=ca.*bi1;
end;
if(x <= 9.0d0);
ct=-(log(x./2.0d0)+el);
bk0=0.0d0;
w0=0.0d0;
r=1.0d0;
for  k=1:50;
w0=w0+1.0d0./k;
r=0.25d0.*r./(k.*k).*x2;
bk0=bk0+r.*(w0+ct);
if(abs((bk0-ww)./bk0)< 1.0d-15)break; end;
ww=bk0;
end;
bk0=bk0+ct;
else;
a1(:)=[0.125d0,0.2109375d0,1.0986328125d0,1.1775970458984d01,2.1461706161499d02,5.9511522710323d03,2.3347645606175d05,1.2312234987631d07];
cb=0.5d0./x;
xr2=1.0d0./x2;
bk0=1.0d0;
for  k=1:8;
bk0=bk0+a1(k).*xr2.^k;
end;  k=8+1;
bk0=cb.*bk0./bi0;
end;
bk1=(1.0d0./x-bi1.*bk0)./bi0;
di0=bi1;
di1=bi0-bi1./x;
dk0=-bk1;
dk1=-bk0-bk1./x;
return;
end

