function mittika
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes the integral of[I0(t)-1]/t
%                with respect to t from 0 to x and K0(t)/t with
%                respect to t from x to ì using subroutine ITTIKA
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTI --- Integration of[I0(t)-1]/t from 0 to x
%                TTK --- Integration of K0(t)/t from x to ì
%       Example:
%                   x[1-I0(t)]/tdt     K0(t)/tdt
%                ---------------------------------------
%                  5.0   .71047763D+01   .58635626D-03
%                 10.0   .34081537D+03   .15629282D-05
%                 15.0   .25437619D+05   .59837472D-08
%                 20.0   .23673661D+07   .26790545D-10
%                 25.0   .24652751D+09   .13100706D-12
%       ============================================================
x=[];tti=[];ttk=[];
  x=0;
 tti=0;
 ttk=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x[1-i0(t)]/tdt     k0(t)/tdt');
fprintf(1,'%s \n','---------------------------------------');
[x,tti,ttk]=ittika(x,tti,ttk);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,2) ' \n'],x,tti,ttk);
%format(1x,f5.1,2d16.8);
end
function [x,tti,ttk]=ittika(x,tti,ttk,varargin);
%       =========================================================
%       Purpose: Integrate[I0(t)-1]/t with respect to t from 0
%                to x, and K0(t)/t with respect to t from x to ì
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTI --- Integration of[I0(t)-1]/t from 0 to x
%                TTK --- Integration of K0(t)/t from x to ì
%       =========================================================
 c=zeros(1,8);
pi=3.141592653589793d0;
el=.5772156649015329d0;
c(:)=[1.625d0,4.1328125d0,1.45380859375d+1,6.553353881835d+1,3.6066157150269d+2,2.3448727161884d+3,1.7588273098916d+4,1.4950639538279d+5];
if(x == 0.0d0);
tti=0.0d0;
ttk=1.0d+300;
return;
end;
if(x < 40.0d0);
tti=1.0d0;
r=1.0d0;
for  k=2:50;
r=.25d0.*r.*(k-1.0d0)./(k.*k.*k).*x.*x;
tti=tti+r;
if(abs(r./tti)< 1.0d-12)break; end;
end;
tti=tti.*.125d0.*x.*x;
else;
tti=1.0d0;
r=1.0d0;
for  k=1:8;
r=r./x;
tti=tti+c(k).*r;
end;  k=8+1;
rc=x.*sqrt(2.0d0.*pi.*x);
tti=tti.*exp(x)./rc;
end;
if(x <= 12.0d0);
e0=(.5d0.*log(x./2.0d0)+el).*log(x./2.0d0)+pi.*pi./24.0d0+.5d0.*el.*el;
b1=1.5d0-(el+log(x./2.0d0));
rs=1.0d0;
r=1.0d0;
for  k=2:50;
r=.25d0.*r.*(k-1.0d0)./(k.*k.*k).*x.*x;
rs=rs+1.0d0./k;
r2=r.*(rs+1.0d0./(2.0d0.*k)-(el+log(x./2.0d0)));
b1=b1+r2;
if(abs(r2./b1)< 1.0d-12)break; end;
end;
ttk=e0-.125d0.*x.*x.*b1;
else;
ttk=1.0d0;
r=1.0d0;
for  k=1:8;
r=-r./x;
ttk=ttk+c(k).*r;
end;  k=8+1;
rc=x.*sqrt(2.0d0./pi.*x);
ttk=ttk.*exp(-x)./rc;
end;
return;
end

