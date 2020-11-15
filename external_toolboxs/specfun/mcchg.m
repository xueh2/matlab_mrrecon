function mcchg
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ===========================================================
%     Purpose: This program computes confluent hypergeometric
%     function M(a,b,z)with real parameters a, b, and
%     a complex argument z using subroutine CCHG
%     Input :  a --- Parameter
%     b --- Parameter
%     z --- Complex argument
%     Output:  CHG --- M(a,b,z)
%     Examples:
%     a      b        z        Re[M(a,b,z)]Im[M(a,b,z)]
%     -------------------------------------------------------
%     3.3   4.25    10 + 0i    .61677489D+04    0
%     3.3   4.25    25 + 0i    .95781835D+10  -.15738228D-03
%     3.3   4.25     3 -  i    .75828716D+01  -.86815474D+01
%     3.3   4.25    15 +10i   -.58313765D+06  -.48195426D+05
%     ===========================================================
a=[];b=[];z=[];chg=[];
fprintf(1,'%s \n','please enter a, b, x and y(z=x+iy)');
%     READ(*,*)A,B,X,Y
a=3.3;
b=4.25;
x=15.0;
y=10.0;
fprintf(1,[repmat(' ',1,1),'a =','%5.1g',',  ','b =','%5.1g',',  ','x =','%5.1g',',  ','y =','%5.1g' ' \n'],a,b,x,y);
z=complex(x,y);
[a,b,z,chg]=cchg(a,b,z,chg);
fprintf(1,[repmat(' ',1,10),'m(a,b,z)=','%18.8g',' + i ','%18.8g' ' \n'],chg);
fprintf(1,'%0.15g \n', chg);
%format(10x,',d18.8,' + i ',d18.8);
%format(1x,',f5.1,',  ',',f5.1,',  ',',f5.1, ',  ',',f5.1);
end
function [a,b,z,chg]=cchg(a,b,z,chg,varargin);
%     ===================================================
%     Purpose: Compute confluent hypergeometric function
%     M(a,b,z)with real parameters a, b and a
%     complex argument z
%     Input :  a --- Parameter
%     b --- Parameter
%     z --- Complex argument
%     Output:  CHG --- M(a,b,z)
%     Routine called: GAMMA for computing gamma function
%     ===================================================
g1=[];g2=[];ba=[];g3=[];
chw=0.0;
pi=3.141592653589793d0;
ci=complex(0.0d0,1.0d0);
a0=a;
a1=a;
z0=z;
if(b == 0.0|b == -fix(abs(b)));
chg=complex(1.0d+300,0.0d0);
elseif(a == 0.0d0|z == 0.0d0);
chg=complex(1.0d0,0.0d0);
elseif(a == -1.0d0);
chg=1.0d0-z./b;
elseif(a == b);
chg=exp(z);
elseif(a-b == 1.0d0);
chg=(1.0d0+z./b).*exp(z);
elseif(a == 1.0d0&b == 2.0d0);
chg=(exp(z)-1.0d0)./z;
elseif(a == fix(a)&a < 0.0d0);
m=fix(-a);
cr=complex(1.0d0,0.0d0);
chg=complex(1.0d0,0.0d0);
for  k=1:m;
cr=cr.*(a+k-1.0d0)./k./(b+k-1.0d0).*z;
chg=chg+cr;
end;  k=m+1;
else;
x0=real(z);
if(x0 < 0.0d0);
a=b-a;
a0=a;
z=-z;
end;
if(a < 2.0d0)nl=0; end;
if(a >= 2.0d0);
nl=1;
la=fix(a);
a=a-la-1.0d0;
end;
for  n=0:nl;
if(a0 >= 2.0d0)a=a+1.0d0; end;
if(abs(z)< 20.0d0+abs(b)|a < 0.0d0);
chg=complex(1.0d0,0.0d0);
crg=complex(1.0d0,0.0d0);
for  j=1:500;
crg=crg.*(a+j-1.0d0)./(j.*(b+j-1.0d0)).*z;
chg=chg+crg;
if(abs((chg-chw)./chg)< 1.d-15)break; end;
chw=chg;
end;
else;
[a,g1]=gamma(a,g1);
[b,g2]=gamma(b,g2);
ba=b-a;
[ba,g3]=gamma(ba,g3);
cs1=complex(1.0d0,0.0d0);
cs2=complex(1.0d0,0.0d0);
cr1=complex(1.0d0,0.0d0);
cr2=complex(1.0d0,0.0d0);
for  i=1:8;
cr1=-cr1.*(a+i-1.0d0).*(a-b+i)./(z.*i);
cr2=cr2.*(b-a+i-1.0d0).*(i-a)./(z.*i);
cs1=cs1+cr1;
cs2=cs2+cr2;
end;  i=8+1;
x=real(z);
y=imag(z);
if(x == 0.0&y >= 0.0);
phi=0.5d0.*pi;
elseif(x == 0.0&y <= 0.0);
phi=-0.5d0.*pi;
else;
phi=atan(y./x);
end;
if(phi > -0.5.*pi&phi < 1.5.*pi)ns=1; end;
if(phi > -1.5.*pi&phi <= -0.5.*pi)ns=-1; end;
cfac=exp(ns.*ci.*pi.*a);
if(y == 0.0d0)cfac=cos(pi.*a); end;
chg1=g2./g3.*z.^(-a).*cfac.*cs1;
chg2=g2./g1.*exp(z).*z.^(a-b).*cs2;
chg=chg1+chg2;
end;
if(n == 0)cy0=chg; end;
if(n == 1)cy1=chg; end;
end;
if(a0 >= 2.0d0);
for  i=1:la-1;
chg=((2.0d0.*a-b+z).*cy1+(b-a).*cy0)./a;
cy0=cy1;
cy1=chg;
a=a+1.0d0;
end;  i=la-1+1;
end;
if(x0 < 0.0d0)chg=chg.*exp(-z); end;
end;
a=a1;
z=z0;
return;
end
function [x,ga]=gamma(x,ga,varargin);
%     ==================================================
%     Purpose: Compute gamma function ג(x)
%     Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%     Output:  GA --- ג(x)
%     ==================================================
 g=zeros(1,26);
pi=3.141592653589793d0;
if(x == fix(x));
if(x > 0.0d0);
ga=1.0d0;
m1=x-1;
for  k=2:m1;
ga=ga.*k;
end;  k=m1+1;
else;
ga=1.0d+300;
end;
else;
if(abs(x)> 1.0d0);
z=abs(x);
m=fix(z);
r=1.0d0;
for  k=1:m;
r=r.*(z-k);
end;  k=m+1;
z=z-m;
else;
z=x;
end;
g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
gr=g(26);
for  k=25:-1:1;
gr=gr.*z+g(k);
end;  k=1-1;
ga=1.0d0./(gr.*z);
if(abs(x)> 1.0d0);
ga=ga.*r;
if(x < 0.0d0)ga=-pi./(x.*ga.*sin(pi.*x)); end;
end;
end;
return;
end

