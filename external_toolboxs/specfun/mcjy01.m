function mcjy01
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ================================================================
%       Purpose: This program computes Bessel functions J0(z), J1(z),
%                Y0(z), Y1(z), and their derivatives for a complex
%                argument using subroutine CJY01
%       Input :  z --- Complex argument
%       Output:  CBJ0 --- J0(z)
%                CDJ0 --- J0'(z)
%                CBJ1 --- J1(z)
%                CDJ1 --- J1'(z)
%                CBY0 --- Y0(z)
%                CDY0 --- Y0'(z)
%                CBY1 --- Y1(z)
%                CDY1 --- Y1'(z)
%       Example: z =  4.0 + i  2.0
%     n     Re[Jn(z)]Im[Jn(z)]Re[Jn'(z)]Im[Jn'(z)]
%   --------------------------------------------------------------------
%     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
%     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
%     n     Re[Yn(z)]Im[Yn(z)]Re[Yn'(z)]Im[Yn'(z)]
%   --------------------------------------------------------------------
%     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
%     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
%       ================================================================
z=[];cbj0=[];cdj0=[];cbj1=[];cdj1=[];cby0=[];cdy0=[];cby1=[];cdy1=[];
fprintf(1,'%s \n','  please enter x,y(z=x+iy)');
%        READ(*,*)X,Y
x=4.0;
y=2.0;
z=complex(x,y);
fprintf(1,[repmat(' ',1,3),'z =','%6.2g',' + i','%6.2g' ' \n'],x,y);
[z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1]=cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[jn(z)]im[jn(z)]');fprintf(1,'%s \n','       re[jn''(z)]im[jn''(z)]');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','-------------------------------');
fprintf(1,[repmat(' ',1,3),'0',repmat(' ',1,2),repmat('%16.8g',1,4) ' \n'],cbj0,cdj0);
fprintf(1,[repmat(' ',1,3),'1',repmat(' ',1,2),repmat('%16.8g',1,4) ' \n'],cbj1,cdj1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  n      re[yn(z)]im[yn(z)]');fprintf(1,'%s \n','       re[yn''(z)]im[yn''(z)]');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','-------------------------------');
fprintf(1,[repmat(' ',1,3),'0',repmat(' ',1,2),repmat('%16.8g',1,4) ' \n'],cby0,cdy0);
fprintf(1,[repmat(' ',1,3),'1',repmat(' ',1,2),repmat('%16.8g',1,4) ' \n'],cby1,cdy1);
%format(3x,'0',2x,4d16.8);
%format(3x,'1',2x,4d16.8);
%format(3x,,f6.2,' + i',f6.2);
end
function [z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1]=cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1,varargin);
%       =======================================================
%       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z),
%                Y1(z), and their derivatives for a complex
%                argument
%       Input :  z --- Complex argument
%       Output:  CBJ0 --- J0(z)
%                CDJ0 --- J0'(z)
%                CBJ1 --- J1(z)
%                CDJ1 --- J1'(z)
%                CBY0 --- Y0(z)
%                CDY0 --- Y0'(z)
%                CBY1 --- Y1(z)
%                CDY1 --- Y1'(z)
%       =======================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,12);
b1=zeros(1,12);
pi=3.141592653589793d0;
el=0.5772156649015329d0;
rp2=2.0d0./pi;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z2=z.*z;
z1=z;
if(a0 == 0.0d0);
cbj0=complex(1.0d0,0.0d0);
cbj1=complex(0.0d0,0.0d0);
cdj0=complex(0.0d0,0.0d0);
cdj1=complex(0.5d0,0.0d0);
cby0=-complex(1.0d300,0.0d0);
cby1=-complex(1.0d300,0.0d0);
cdy0=complex(1.0d300,0.0d0);
cdy1=complex(1.0d300,0.0d0);
return;
end;
if(real(z)< 0.0)z1=-z; end;
if(a0 <= 12.0);
cbj0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*k);
cbj0=cbj0+cr;
if(abs(cr)< abs(cbj0).*1.0d-15)break; end;
end;
cbj1=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*(k+1.0d0));
cbj1=cbj1+cr;
if(abs(cr)< abs(cbj1).*1.0d-15)break; end;
end;
cbj1=0.5d0.*z1.*cbj1;
w0=0.0d0;
cr=complex(1.0d0,0.0d0);
cs=complex(0.0d0,0.0d0);
for  k=1:40;
w0=w0+1.0d0./k;
cr=-0.25d0.*cr./(k.*k).*z2;
cp=cr.*w0;
cs=cs+cp;
if(abs(cp)< abs(cs).*1.0d-15)break; end;
end;
cby0=rp2.*(log(z1./2.0d0)+el).*cbj0-rp2.*cs;
w1=0.0d0;
cr=complex(1.0d0,0.0d0);
cs=complex(1.0d0,0.0d0);
for  k=1:40;
w1=w1+1.0d0./k;
cr=-0.25d0.*cr./(k.*(k+1)).*z2;
cp=cr.*(2.0d0.*w1+1.0d0./(k+1.0d0));
cs=cs+cp;
if(abs(cp)< abs(cs).*1.0d-15)break; end;
end;
cby1=rp2.*((log(z1./2.0d0)+el).*cbj1-1.0d0./z1-.25d0.*z1.*cs);
else;
a(:)=[-.703125d-01,.112152099609375d+00,-.5725014209747314d+00,.6074042001273483d+01,-.1100171402692467d+03,.3038090510922384d+04,-.1188384262567832d+06,.6252951493434797d+07,-.4259392165047669d+09,.3646840080706556d+11,-.3833534661393944d+13,.4854014686852901d+15];
b(:)=[.732421875d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02,.5513358961220206d+03,-.1825775547429318d+05,.8328593040162893d+06,-.5006958953198893d+08,.3836255180230433d+10,-.3649010818849833d+12,.4218971570284096d+14,-.5827244631566907d+16];
a1(:)=[.1171875d+00,-.144195556640625d+00,.6765925884246826d+00,-.6883914268109947d+01,.1215978918765359d+03,-.3302272294480852d+04,.1276412726461746d+06,-.6656367718817688d+07,.4502786003050393d+09,-.3833857520742790d+11,.4011838599133198d+13,-.5060568503314727d+15];
b1(:)=[-.1025390625d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02,-.6038440767050702d+03,.1971837591223663d+05,-.8902978767070678d+06,.5310411010968522d+08,-.4043620325107754d+10,.3827011346598605d+12,-.4406481417852278d+14,.6065091351222699d+16];
k0=12;
if(a0 >= 35.0)k0=10; end;
if(a0 >= 50.0)k0=8; end;
ct1=z1-.25d0.*pi;
cp0=complex(1.0d0,0.0d0);
for  k=1:k0;
cp0=cp0+a(k).*z1.^(-2.*k);
end;  k=k0+1;
cq0=-0.125d0./z1;
for  k=1:k0;
cq0=cq0+b(k).*z1.^(-2.*k-1);
end;  k=k0+1;
cu=sqrt(rp2./z1);
cbj0=cu.*(cp0.*cos(ct1)-cq0.*sin(ct1));
cby0=cu.*(cp0.*sin(ct1)+cq0.*cos(ct1));
ct2=z1-.75d0.*pi;
cp1=complex(1.0d0,0.0d0);
for  k=1:k0;
cp1=cp1+a1(k).*z1.^(-2.*k);
end;  k=k0+1;
cq1=0.375d0./z1;
for  k=1:k0;
cq1=cq1+b1(k).*z1.^(-2.*k-1);
end;  k=k0+1;
cbj1=cu.*(cp1.*cos(ct2)-cq1.*sin(ct2));
cby1=cu.*(cp1.*sin(ct2)+cq1.*cos(ct2));
end;
if(real(z)< 0.0);
if(imag(z)< 0.0)cby0=cby0-2.0d0.*ci.*cbj0; end;
if(imag(z)> 0.0)cby0=cby0+2.0d0.*ci.*cbj0; end;
if(imag(z)< 0.0)cby1=-(cby1-2.0d0.*ci.*cbj1); end;
if(imag(z)> 0.0)cby1=-(cby1+2.0d0.*ci.*cbj1); end;
cbj1=-cbj1;
end;
cdj0=-cbj1;
cdj1=cbj0-1.0d0./z.*cbj1;
cdy0=-cby1;
cdy1=cby0-1.0d0./z.*cby1;
return;
end

