function mitika
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program evaluates the integral of modified
%                Bessel functions I0(t)and K0(t)with respect to t
%                from 0 to x using subroutine ITIKA
%       Input :  x  --- Upper limit of the integral(x ò 0)
%       Output:  TI --- Integration of I0(t)from 0 to x
%                TK --- Integration of K0(t)from 0 to x
%       Example:
%                    x         I0(t)dt         K0(t)dt
%                 --------------------------------------
%                   5.0    .31848668D+02     1.56738739
%                  10.0    .29930445D+04     1.57077931
%                  15.0    .35262048D+06     1.57079623
%                  20.0    .44758593D+08     1.57079633
%                  25.0    .58991731D+10     1.57079633
%       ============================================================
x=[];ti=[];tk=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x         i0(t)dt         k0(t)dt');
fprintf(1,'%s \n',' --------------------------------------');
[x,ti,tk]=itika(x,ti,tk);
fprintf(1,[repmat(' ',1,1),'%5.1g','%17.8g','%15.8g' ' \n'],x,ti,tk);
%format(1x,f5.1,d17.8,f15.8);
end
function [x,ti,tk]=itika(x,ti,tk,varargin);
%       =======================================================
%       Purpose: Integrate modified Bessel functions I0(t)and
%                K0(t)with respect to t from 0 to x
%       Input :  x  --- Upper limit of the integral(x ò 0)
%       Output:  TI --- Integration of I0(t)from 0 to x
%                TK --- Integration of K0(t)from 0 to x
%       =======================================================
 a=zeros(1,10);
pi=3.141592653589793d0;
el=.5772156649015329d0;
a(:)=[.625d0,1.0078125d0,2.5927734375d0,9.1868591308594d0,4.1567974090576d+1,2.2919635891914d+2,1.491504060477d+3,1.1192354495579d+4,9.515939374212d+4,9.0412425769041d+5];
if(x == 0.0d0);
ti=0.0d0;
tk=0.0d0;
return;
elseif(x < 20.0d0);
x2=x.*x;
ti=1.0d0;
r=1.0d0;
for  k=1:50;
r=.25d0.*r.*(2.*k-1.0d0)./(2.*k+1.0d0)./(k.*k).*x2;
ti=ti+r;
if(abs(r./ti)< 1.0d-12)break; end;
end;
ti=ti.*x;
else;
ti=1.0d0;
r=1.0d0;
for  k=1:10;
r=r./x;
ti=ti+a(k).*r;
end;  k=10+1;
rc1=1.0d0./sqrt(2.0d0.*pi.*x);
ti=rc1.*exp(x).*ti;
end;
if(x < 12.0d0);
e0=el+log(x./2.0d0);
b1=1.0d0-e0;
b2=0.0d0;
rs=0.0d0;
r=1.0d0;
for  k=1:50;
r=.25d0.*r.*(2.*k-1.0d0)./(2.*k+1.0d0)./(k.*k).*x2;
b1=b1+r.*(1.0d0./(2.*k+1)-e0);
rs=rs+1.0d0./k;
b2=b2+r.*rs;
tk=b1+b2;
if(abs((tk-tw)./tk)< 1.0d-12)break; end;
tw=tk;
end;
tk=tk.*x;
else;
tk=1.0d0;
r=1.0d0;
for  k=1:10;
r=-r./x;
tk=tk+a(k).*r;
end;  k=10+1;
rc2=sqrt(pi./(2.0d0.*x));
tk=pi./2.0d0-rc2.*tk.*exp(-x);
end;
return;
end

