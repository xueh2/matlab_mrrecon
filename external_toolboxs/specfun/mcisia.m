function mcisia
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ========================================================
%     Purpose: This program computes the cosine and sine
%     integrals using subroutine CISIA
%     Input :  x  --- Argument of Ci(x)and Si(x)
%     Output:  CI --- Ci(x)
%     SI --- Si(x)
%     Example:
%     x         Ci(x)Si(x)
%     ------------------------------------
%     0.0     - ì             .00000000
%     5.0     -.19002975     1.54993124
%     10.0     -.04545643     1.65834759
%     20.0      .04441982     1.54824170
%     30.0     -.03303242     1.56675654
%     40.0      .01902001     1.58698512
%     ========================================================
x=[];ci=[];si=[];
  ci=0;
 si=0;
 x=0;
fprintf(1,'%s \n','please enter x ');
%     READ(*,*)X
x=5.0;
fprintf(1,'%s \n','   x         ci(x)si(x)');
fprintf(1,'%s \n','------------------------------------');
[x,ci,si]=cisia(x,ci,si);
if(x ~= 0.0d0)fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%15.8g',1,2) ' \n'],x,ci,si); end;
if(x == 0.0d0)fprintf(1,[repmat(' ',1,3),' .0',repmat(' ',1,4),' - ì',repmat(' ',1,13),'.00000000' ' \n']); end;
%format(1x,f5.1,2f15.8);
%format(3x,' .0',4x,' - ì',13x,'.00000000');
end
function [x,ci,si]=cisia(x,ci,si,varargin);
%     =============================================
%     Purpose: Compute cosine and sine integrals
%     Si(x)and Ci(x,x ò 0)
%     Input :  x  --- Argument of Ci(x)and Si(x)
%     Output:  CI --- Ci(x)
%     SI --- Si(x)
%     =============================================
 bj=zeros(1,101);
p2=1.570796326794897d0;
el=.5772156649015329d0;
eps=1.0d-15;
x2=x.*x;
if(x == 0.0d0);
ci=-1.0d+300;
si=0.0d0;
elseif(x <= 16.0d0);
xr=-.25d0.*x2;
ci=el+log(x)+xr;
for  k=2:40;
xr=-.5d0.*xr.*(k-1)./(k.*k.*(2.*k-1)).*x2;
ci=ci+xr;
if(abs(xr)< abs(ci).*eps)break; end;
end;
xr=x;
si=x;
for  k=1:40;
xr=-.5d0.*xr.*(2.*k-1)./k./(4.*k.*k+4.*k+1).*x2;
si=si+xr;
if(abs(xr)< abs(si).*eps)return; end;
end;  k=40+1;
elseif(x <= 32.0d0);
m=fix(47.2+.82.*x);
xa1=0.0d0;
xa0=1.0d-100;
for  k=m:-1:1;
xa=4.0d0.*k.*xa0./x-xa1;
bj(k)=xa;
xa1=xa0;
xa0=xa;
end;  k=1-1;
xs=bj(1);
for  k=3:2:m;
xs=xs+2.0d0.*bj(k);
end;  k=m+1;
bj(1)=bj(1)./xs;
for  k=2:m;
bj(k)=bj(k)./xs;
end;  k=m+1;
xr=1.0d0;
xg1=bj(1);
for  k=2:m;
xr=.25d0.*xr.*(2.0.*k-3.0).^2./((k-1.0).*(2.0.*k-1.0).^2).*x;
xg1=xg1+bj(k).*xr;
end;  k=m+1;
xr=1.0d0;
xg2=bj(1);
for  k=2:m;
xr=.25d0.*xr.*(2.0.*k-5.0).^2./((k-1.0).*(2.0.*k-3.0).^2).*x;
xg2=xg2+bj(k).*xr;
end;  k=m+1;
xcs=cos(x./2.0d0);
xss=sin(x./2.0d0);
ci=el+log(x)-x.*xss.*xg1+2.*xcs.*xg2-2.*xcs.*xcs;
si=x.*xcs.*xg1+2.*xss.*xg2-sin(x);
else;
xr=1.0d0;
xf=1.0d0;
for  k=1:9;
xr=-2.0d0.*xr.*k.*(2.*k-1)./x2;
xf=xf+xr;
end;  k=9+1;
xr=1.0d0./x;
xg=xr;
for  k=1:8;
xr=-2.0d0.*xr.*(2.*k+1).*k./x2;
xg=xg+xr;
end;  k=8+1;
ci=xf.*sin(x)./x-xg.*cos(x)./x;
si=p2-xf.*cos(x)./x-xg.*sin(x)./x;
end;
return;
end

