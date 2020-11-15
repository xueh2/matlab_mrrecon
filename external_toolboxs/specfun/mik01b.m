function mik01b
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
%                using subroutine IK01B
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
%        1.0   .126607D+01   .565159D+00   .565159D+00   .700907D+00
%       10.0   .281572D+04   .267099D+04   .267099D+04   .254862D+04
%       20.0   .435583D+08   .424550D+08   .424550D+08   .414355D+08
%       30.0   .781672D+12   .768532D+12   .768532D+12   .756055D+12
%       40.0   .148948D+17   .147074D+17   .147074D+17   .145271D+17
%       50.0   .293255D+21   .290308D+21   .290308D+21   .287449D+21
%         x      K0(x)K0'(x)K1(x)K1'(x)
%       -------------------------------------------------------------
%        1.0   .421024D+00  -.601907D+00   .601907D+00  -.102293D+01
%       10.0   .177801D-04  -.186488D-04   .186488D-04  -.196449D-04
%       20.0   .574124D-09  -.588306D-09   .588306D-09  -.603539D-09
%       30.0   .213248D-13  -.216773D-13   .216773D-13  -.220474D-13
%       40.0   .839286D-18  -.849713D-18   .849713D-18  -.860529D-18
%       50.0   .341017D-22  -.344410D-22   .344410D-22  -.347905D-22
%       =============================================================
x=[];bi0=[];di0=[];bi1=[];di1=[];bk0=[];dk0=[];bk1=[];dk1=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=1.0;
fprintf(1,[repmat(' ',1,3),'x =','%5.1g' ' \n'],x);
fprintf(1,'%s ','  x       i0(x)i0''(x)i1(x)');fprintf(1,'%s \n','          i1''(x)');
fprintf(1,'%s ','-------------------------------------------');fprintf(1,'%s \n','----------------------');
[x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1]=ik01b(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1);
fprintf(1,[repmat(' ',1,1),'%4.1g',repmat('%15.7g',1,4) ' \n'],x,bi0,di0,bi1,di1);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','  x       k0(x)k0''(x)k1(x)');fprintf(1,'%s \n','          k1''(x)');
fprintf(1,'%s ','-------------------------------------------');fprintf(1,'%s \n','----------------------');
fprintf(1,[repmat(' ',1,1),'%4.1g',repmat('%15.7g',1,4) ' \n'],x,bk0,dk0,bk1,dk1);
%format(',f5.1);
%format(1x,f4.1,4d15.7);
end
function [x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1]=ik01b(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1,varargin);
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
elseif(x <= 3.75d0);
t=x./3.75d0;
t2=t.*t;
bi0=(((((.0045813d0.*t2+.0360768d0).*t2+.2659732d0).*t2+1.2067492d0).*t2+3.0899424d0).*t2 +3.5156229d0).*t2+1.0d0;
bi1=x.*((((((.00032411d0.*t2+.00301532d0).*t2+.02658733d0).*t2+.15084934d0).*t2+.51498869d0).*t2+.87890594d0).*t2+.5d0);
else;
t=3.75d0./x;
bi0=((((((((.00392377d0.*t-.01647633d0).*t+.02635537d0).*t-.02057706d0).*t+.916281d-2).*t-.157565d-2).*t+.225319d-2).*t+.01328592d0).*t +.39894228d0).*exp(x)./sqrt(x);
bi1=((((((((-.420059d-2.*t+.01787654d0).*t-.02895312d0).*t+.02282967d0).*t-.01031555d0).*t+.163801d-2).*t-.00362018d0).*t-.03988024d0).*t +.39894228d0).*exp(x)./sqrt(x);
end;
if(x <= 2.0d0);
t=x./2.0d0;
t2=t.*t;
bk0=(((((.0000074d0.*t2+.0001075d0).*t2+.00262698d0).*t2+.0348859d0).*t2+.23069756d0).*t2+.4227842d0).*t2-.57721566d0-bi0.*log(t);
bk1=((((((-.00004686d0.*t2-.00110404d0).*t2-.01919402d0).*t2-.18156897d0).*t2-.67278579d0).*t2+.15443144d0).*t2+1.0d0)./x+bi1.*log(t);
else;
t=2.0d0./x;
t2=t.*t;
bk0=((((((.00053208d0.*t-.0025154d0).*t+.00587872d0).*t-.01062446d0).*t+.02189568d0).*t-.07832358d0).*t+1.25331414d0).*exp(-x)./sqrt(x);
bk1=((((((-.00068245d0.*t+.00325614d0).*t-.00780353d0).*t+.01504268d0).*t-.0365562d0).*t+.23498619d0).*t+1.25331414d0).*exp(-x)./sqrt(x);
end;
di0=bi1;
di1=bi0-bi1./x;
dk0=-bk1;
dk1=-bk0-bk1./x;
return;
end

