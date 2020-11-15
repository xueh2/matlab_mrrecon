function me1z
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes the complex exponential
%                integral E1(z)using subroutine E1Z
%       Example:
%                     z            Re[E1(z)]Im[E1(z)]
%                -----------------------------------------------
%                 3.0    2.0    -.90959209D-02   -.69001793D-02
%                 3.0   -2.0    -.90959209D-02    .69001793D-02
%                -3.0    2.0    -.28074891D+01    .59603353D+01
%                -3.0   -2.0    -.28074891D+01   -.59603353D+01
%                25.0   10.0    -.29302080D-12    .40391222D-12
%                25.0  -10.0    -.29302080D-12   -.40391222D-12
%               -25.0   10.0     .27279957D+10   -.49430610D+09
%               -25.0  -10.0     .27279957D+10    .49430610D+09
%       =========================================================
z=[];ce1=[];
fprintf(1,'%s \n','please enter x and y(z =x+iy)');
%        READ(*,*)X,Y
x=3.0;
y=2.0;
z=complex(x,y);
[z,ce1]=e1z(z,ce1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','       z           re[e1(z)]im[e1(z)]');
fprintf(1,'%s \n',' -----------------------------------------------');
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat(' ',1,2),'%5.1g',repmat(' ',1,1),repmat('%17.8g',1,2) ' \n'],x,y,ce1);
%format(1x,f5.1,2x,f5.1,1x,2d17.8);
end
function [z,ce1]=e1z(z,ce1,varargin);
%       ====================================================
%       Purpose: Compute complex exponential integral E1(z)
%       Input :  z   --- Argument of E1(z)
%       Output:  CE1 --- E1(z)
%       ====================================================
pi=3.141592653589793d0;
el=0.5772156649015328d0;
x=real(z);
a0=abs(z);
if(a0 == 0.0d0);
ce1=complex(1.0d+300,0.0d0);
elseif(a0 <= 10.0|x < 0.0&a0 < 20.0);
ce1=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:150;
cr=-cr.*k.*z./(k+1.0d0).^2;
ce1=ce1+cr;
if(abs(cr)<= abs(ce1).*1.0d-15)break; end;
end;
ce1=-el-log(z)+z.*ce1;
else;
ct0=complex(0.0d0,0.0d0);
for  k=120:-1:1;
ct0=k./(1.0d0+k./(z+ct0));
end;  k=1-1;
ct=1.0d0./(z+ct0);
ce1=exp(-z).*ct;
if(x <= 0.0&imag(z)== 0.0)ce1=ce1-pi.*complex(0.0d0,1.0d0); end;
end;
return;
end

