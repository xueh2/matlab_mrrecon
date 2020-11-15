function mittjyb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program computes the integral of[1-J0(t)]/t
%                with respect to t from 0 to x and Y0(t)/t with
%                respect to t from x to ì using subroutine ITTJYB
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTJ --- Integration of[1-J0(t)]/t from 0 to x
%                TTY --- Integration of Y0(t)/t from x to ì
%       Example:
%                  x[1-J0(t)]/tdt       Y0(t)/tdt
%                ----------------------------------------
%                 5.0     .1540347D+01    -.4632208D-01
%                10.0     .2177866D+01    -.2298791D-01
%                15.0     .2578551D+01     .3857453D-03
%                20.0     .2877311D+01     .8503154D-02
%                25.0     .3108231D+01     .3526339D-02
%       ===========================================================
x=[];ttj=[];tty=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x[1-j0(t)]/tdt      y0(t)/tdt');
fprintf(1,'%s \n','----------------------------------------');
[x,ttj,tty]=ittjyb(x,ttj,tty);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%17.7g',1,2) ' \n'],x,ttj,tty);
%format(1x,f5.1,2d17.7);
end
function [x,ttj,tty]=ittjyb(x,ttj,tty,varargin);
%       ==========================================================
%       Purpose: Integrate[1-J0(t)]/t with respect to t from 0
%                to x, and Y0(t)/t with respect to t from x to ì
%       Input :  x   --- Variable in the limits(x ò 0)
%       Output:  TTJ --- Integration of[1-J0(t)]/t from 0 to x
%                TTY --- Integration of Y0(t)/t from x to ì
%       ==========================================================
pi=3.141592653589793d0;
el=.5772156649015329d0;
if(x == 0.0d0);
ttj=0.0d0;
tty=-1.0d+300;
elseif(x <= 4.0d0);
x1=x./4.0d0;
t=x1.*x1;
ttj=((((((.35817d-4.*t-.639765d-3).*t+.7092535d-2).*t-.055544803d0).*t+.296292677d0).*t-.999999326d0).*t+1.999999936d0).*t;
tty=(((((((-.3546d-5.*t+.76217d-4).*t-.1059499d-2).*t+.010787555d0).*t-.07810271d0).*t+.377255736d0).*t-1.114084491d0).*t+1.909859297d0).*t;
e0=el+log(x./2.0d0);
tty=pi./6.0d0+e0./pi.*(2.0d0.*ttj-e0)-tty;
elseif(x <= 8.0d0);
xt=x+.25d0.*pi;
t1=4.0d0./x;
t=t1.*t1;
f0=(((((.0145369d0.*t-.0666297d0).*t+.1341551d0).*t-.1647797d0).*t+.1608874d0).*t-.2021547d0).*t +.7977506d0;
g0=((((((.0160672d0.*t-.0759339d0).*t+.1576116d0).*t-.1960154d0).*t+.1797457d0).*t-.1702778d0).*t +.3235819d0).*t1;
ttj=(f0.*cos(xt)+g0.*sin(xt))./(sqrt(x).*x);
ttj=ttj+el+log(x./2.0d0);
tty=(f0.*sin(xt)-g0.*cos(xt))./(sqrt(x).*x);
else;
t=8.0d0./x;
xt=x+.25d0.*pi;
f0=(((((.18118d-2.*t-.91909d-2).*t+.017033d0).*t-.9394d-3).*t-.051445d0).*t-.11d-5).*t+.7978846d0;
g0=(((((-.23731d-2.*t+.59842d-2).*t+.24437d-2).*t-.0233178d0).*t+.595d-4).*t+.1620695d0).*t;
ttj=(f0.*cos(xt)+g0.*sin(xt))./(sqrt(x).*x)+el+log(x./2.0d0);
tty=(f0.*sin(xt)-g0.*cos(xt))./(sqrt(x).*x);
end;
return;
end

