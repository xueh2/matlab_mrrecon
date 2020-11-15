function mitikb
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
%                from 0 to x using subroutine ITIKB
%       Input :  x  --- Upper limit of the integral(x ò 0)
%       Output:  TI --- Integration of I0(t)from 0 to x
%                TK --- Integration of K0(t)from 0 to x
%       Example:
%                    x         I0(t)dt         K0(t)dt
%                 -------------------------------------
%                   5.0     .318487D+02       1.567387
%                  10.0     .299305D+04       1.570779
%                  15.0     .352619D+06       1.570796
%                  20.0     .447586D+08       1.570796
%                  25.0     .589919D+10       1.570796
%       ============================================================
x=[];ti=[];tk=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=25.0;
fprintf(1,'%s \n','   x         i0(t)dt         k0(t)dt');
fprintf(1,'%s \n','--------------------------------------');
[x,ti,tk]=itikb(x,ti,tk);
fprintf(1,[repmat(' ',1,1),'%5.1g','%16.6g','%15.6g' ' \n'],x,ti,tk);
%format(1x,f5.1,d16.6,f15.6);
end
function [x,ti,tk]=itikb(x,ti,tk,varargin);
%       =======================================================
%       Purpose: Integrate Bessel functions I0(t)and K0(t)
%                with respect to t from 0 to x
%       Input :  x  --- Upper limit of the integral(x ò 0)
%       Output:  TI --- Integration of I0(t)from 0 to x
%                TK --- Integration of K0(t)from 0 to x
%       =======================================================
pi=3.141592653589793d0;
if(x == 0.0d0);
ti=0.0d0;
elseif(x < 5.0d0);
t1=x./5.0d0;
t=t1.*t1;
ti=((((((((.59434d-3.*t+.4500642d-2).*t+.044686921d0).*t+.300704878d0).*t+1.471860153d0).*t+4.844024624d0).*t+9.765629849d0).*t +10.416666367d0).*t+5.0d0).*t1;
elseif(x >= 5.0&x <= 8.0d0);
t=5.0d0./x;
ti=(((-.015166d0.*t-.0202292d0).*t+.1294122d0).*t -.0302912d0).*t+.4161224d0;
ti=ti.*exp(x)./sqrt(x);
else;
t=8.0d0./x;
ti=(((((-.0073995d0.*t+.017744d0).*t-.0114858d0).*t+.55956d-2).*t+.59191d-2).*t+.0311734d0).*t +.3989423d0;
ti=ti.*exp(x)./sqrt(x);
end;
if(x == 0.0d0);
tk=0.0d0;
elseif(x <= 2.0d0);
t1=x./2.0d0;
t=t1.*t1;
tk=((((((.116d-5.*t+.2069d-4).*t+.62664d-3).*t+.01110118d0).*t+.11227902d0).*t+.50407836d0).*t +.84556868d0).*t1;
tk=tk-log(x./2.0d0).*ti;
elseif(x > 2.0&x <= 4.0d0);
t=2.0d0./x;
tk=(((.0160395d0.*t-.0781715d0).*t+.185984d0).*t -.3584641d0).*t+1.2494934d0;
tk=pi./2.0d0-tk.*exp(-x)./sqrt(x);
elseif(x > 4.0&x <= 7.0d0);
t=4.0d0./x;
tk=(((((.37128d-2.*t-.0158449d0).*t+.0320504d0).*t-.0481455d0).*t+.0787284d0).*t-.1958273d0).*t +1.2533141d0;
tk=pi./2.0d0-tk.*exp(-x)./sqrt(x);
else;
t=7.0d0./x;
tk=(((((.33934d-3.*t-.163271d-2).*t+.417454d-2).*t-.933944d-2).*t+.02576646d0).*t-.11190289d0).*t +1.25331414d0;
tk=pi./2.0d0-tk.*exp(-x)./sqrt(x);
end;
return;
end

