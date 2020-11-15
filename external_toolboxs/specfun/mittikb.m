function mittikb
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
%                respect to t from x to ì using subroutine ITTIKB
%       Input :  x   --- Upper limit of the integral
%       Output:  TTI --- Integration of[I0(t)-1]/t from 0 to x
%                TTK --- Integration of K0(t)/t from x to ì
%       Example:
%                   x[1-I0(t)]/tdt      K0(t)/tdt
%                ---------------------------------------
%                  5.0     .710478D+01     .586361D-03
%                 10.0     .340811D+03     .156293D-05
%                 15.0     .254373D+05     .598363D-08
%                 20.0     .236735D+07     .267906D-10
%                 25.0     .246534D+09     .131007D-12
%       ============================================================
x=[];tti=[];ttk=[];
  x=0;
 tti=0;
 ttk=0;
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x[1-i0(t)]/tdt      k0(t)/tdt');
fprintf(1,'%s \n','---------------------------------------');
[x,tti,ttk]=ittikb(x,tti,ttk);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.6g',1,2) ' \n'],x,tti,ttk);
%format(1x,f5.1,2d16.6);
end
function [x,tti,ttk]=ittikb(x,tti,ttk,varargin);
%       =========================================================
%       Purpose: Integrate[I0(t)-1]/t with respect to t from 0
%                to x, and K0(t)/t with respect to t from x to ì
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTI --- Integration of[I0(t)-1]/t from 0 to x
%                TTK --- Integration of K0(t)/t from x to ì
%       =========================================================
pi=3.141592653589793d0;
el=.5772156649015329d0;
if(x == 0.0d0);
tti=0.0d0;
elseif(x <= 5.0d0);
x1=x./5.0d0;
t=x1.*x1;
tti=(((((((.1263d-3.*t+.96442d-3).*t+.968217d-2).*t+.06615507d0).*t+.33116853d0).*t+1.13027241d0).*t+2.44140746d0).*t+3.12499991d0).*t;
else;
t=5.0d0./x;
tti=(((((((((2.1945464d0.*t-3.5195009d0).*t-11.9094395d0).*t+40.394734d0).*t-48.0524115d0).*t+28.1221478d0).*t-8.6556013d0).*t+1.4780044d0).*t-.0493843d0).*t+.1332055d0).*t+.3989314d0;
tti=tti.*exp(x)./(sqrt(x).*x);
end;
if(x == 0.0d0);
ttk=1.0d+300;
elseif(x <= 2.0d0);
t1=x./2.0d0;
t=t1.*t1;
ttk=(((((.77d-6.*t+.1544d-4).*t+.48077d-3).*t+.925821d-2).*t+.10937537d0).*t+.74999993d0).*t;
e0=el+log(x./2.0d0);
ttk=pi.*pi./24.0d0+e0.*(.5d0.*e0+tti)-ttk;
elseif(x <= 4.0d0);
t=2.0d0./x;
ttk=(((.06084d0.*t-.280367d0).*t+.590944d0).*t -.850013d0).*t+1.234684d0;
ttk=ttk.*exp(-x)./(sqrt(x).*x);
else;
t=4.0d0./x;
ttk=(((((.02724d0.*t-.1110396d0).*t+.2060126d0).*t-.2621446d0).*t+.3219184d0).*t-.5091339d0).*t +1.2533141d0;
ttk=ttk.*exp(-x)./(sqrt(x).*x);
end;
return;
end

