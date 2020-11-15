function mcva2
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     =============================================================
%     Purpose: This program calculates a specific characteristic
%     value of Mathieu functions using subroutine CVA2
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,...)
%     KD=2 for cem(x,q,m = 1,3,5,...)
%     KD=3 for sem(x,q,m = 1,3,5,...)
%     KD=4 for sem(x,q,m = 2,4,6,...)
%     Output:  A  --- Characteristic value
%     Example: q = 25.0, m = 0,1,2,...,12
%     Characteristic values of Mathieu functions
%     m            a                  b
%     ------------------------------------------
%     0      -40.256779547
%     1      -21.314899691      -40.256778985
%     2       -3.522164727      -21.314860622
%     3       12.964079444       -3.520941527
%     4       27.805240581       12.986489953
%     5       40.050190986       28.062765899
%     6       48.975786716       41.801071292
%     7       57.534689001       55.002957151
%     8       69.524065166       69.057988351
%     9       85.076999882       85.023356505
%     10      103.230204804      103.225680042
%     11      123.643012376      123.642713667
%     12      146.207690643      146.207674647
%     =============================================================
kd=[];m=[];q=[];a=[];
fprintf(1,'%s \n','please enter kd, m and q ');
%     READ(*,*)KD,M,Q
kd=1;
m=12;
q=25.0;
[kd,m,q,a]=cva2(kd,m,q,a);
if(kd <= 2);
fprintf(1,[repmat(' ',1,1),'characteristic value of even mathieu function, a' ' \n']);
else;
fprintf(1,[repmat(' ',1,1),'characteristic value of odd mathieu function, b' ' \n']);
end;
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  m         a(or b)');
fprintf(1,'%s \n','------------------------');
fprintf(1,[repmat(' ',1,1),'%3g','%18.8g' ' \n'],m,a);
%format(1x,i3,f18.8);
%format(1x,'characteristic value of even mathieu function, a');
%format(1x,'characteristic value of odd mathieu function, b');
end
function [kd,m,q,a]=cva2(kd,m,q,a,varargin);
%     ======================================================
%     Purpose: Calculate a specific characteristic value of
%     Mathieu functions
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,...)
%     KD=2 for cem(x,q,m = 1,3,5,...)
%     KD=3 for sem(x,q,m = 1,3,5,...)
%     KD=4 for sem(x,q,m = 2,4,6,...)
%     Output:  A  --- Characteristic value
%     Routines called:
%(1)REFINE for finding accurate characteristic
%     value using an iteration method
%(2)CV0 for finding initial characteristic
%     values using polynomial approximation
%(3)CVQM for computing initial characteristic
%     values for q ó 3*m
%(3)CVQL for computing initial characteristic
%     values for q ò m*m
%     ======================================================
q1=[];a1=[];q2=[];a2=[];qq=[];iflag=[];
if(m <= 12|q <= 3.0.*m|q > m.*m);
[kd,m,q,a]=cv0(fix(kd),fix(m),q,a);
if(q ~= 0.0d0);
[kd,m,q,a]=refine(fix(kd),fix(m),q,a,1);
end;
else;
ndiv=10;
delta=(fix(m)-3.0).*fix(m)./ndiv;
if((q-3.0.*m)<=(m.*m-q));
while (1);
nn=fix((q-3.0.*fix(m))./delta)+1;
delta=(q-3.0.*fix(m))./nn;
q1=2.0.*fix(m);
[m,q1,a1]=cvqm(fix(m),q1,a1);
q2=3.0.*fix(m);
[m,q2,a2]=cvqm(fix(m),q2,a2);
qq=3.0.*fix(m);
for  i=1:nn;
qq=qq+delta;
a=(a1.*q2-a2.*q1+(a2-a1).*qq)./(q2-q1);
iflag=1;
if(i == nn)iflag=-1; end;
[kd,m,qq,a,iflag]=refine(fix(kd),fix(m),qq,a,iflag);
q1=q2;
q2=qq;
a1=a2;
a2=a;
end;  i=nn+1;
if(iflag == -10);
ndiv=ndiv.*2;
delta=(fix(m)-3.0).*fix(m)./ndiv;
break;
end;
end;
else;
while (1);
nn=fix((fix(m).*fix(m)-q)./delta)+1;
delta=(fix(m).*fix(m)-q)./nn;
q1=fix(m).*(fix(m)-1.0);
[kd,m,q1,a1]=cvql(fix(kd),fix(m),q1,a1);
q2=fix(m).*fix(m);
[kd,m,q2,a2]=cvql(fix(kd),fix(m),q2,a2);
qq=fix(m).*fix(m);
for  i=1:nn;
qq=qq-delta;
a=(a1.*q2-a2.*q1+(a2-a1).*qq)./(q2-q1);
iflag=1;
if(i == nn)iflag=-1; end;
[kd,m,qq,a,iflag]=refine(fix(kd),fix(m),qq,a,iflag);
q1=q2;
q2=qq;
a1=a2;
a2=a;
end;  i=nn+1;
if(iflag == -10);
ndiv=ndiv.*2;
delta=(fix(m)-3.0).*fix(m)./ndiv;
break;
end;
end;
end;
end;
return;
end
function [kd,m,q,a,iflag]=refine(kd,m,q,a,iflag,varargin);
%     =====================================================
%     Purpose: calculate the accurate characteristic value
%     by the secant method
%     Input :  m --- Order of Mathieu functions
%     q --- Parameter of Mathieu functions
%     A --- Initial characteristic value
%     Output:  A --- Refineed characteristic value
%     Routine called:  CVF for computing the value of F for
%     characteristic equation
%     ========================================================
x0=[];mj=[];f0=[];x1=[];f1=[];x=[];f=[];
eps=1.0d-14;
mj=10+fix(m);
x0=a;
[kd,m,q,x0,mj,f0]=cvf(fix(kd),fix(m),q,x0,mj,f0);
x1=1.002.*a;
[kd,m,q,x1,mj,f1]=cvf(fix(kd),fix(m),q,x1,mj,f1);
for  it=1:100;
mj=mj+1;
x=x1-(x1-x0)./(1.0d0-f0./f1);
[kd,m,q,x,mj,f]=cvf(fix(kd),fix(m),q,x,mj,f);
if(abs(1.0-x1./x)< eps|f == 0.0)break; end;
x0=x1;
f0=f1;
x1=x;
f1=f;
end;
a=x;
return;
end
function [kd,m,q,a,mj,f]=cvf(kd,m,q,a,mj,f,varargin);
%     ======================================================
%     Purpose: Compute the value of F for characteristic
%     equation of Mathieu functions
%     Input :  m --- Order of Mathieu functions
%     q --- Parameter of Mathieu functions
%     A --- Characteristic value
%     Output:  F --- Value of F for characteristic equation
%     ======================================================
b=a;
ic=fix(fix(m)./2);
l=0;
l0=0;
j0=2;
jf=ic;
if(kd == 1)l0=2; end;
if(kd == 1)j0=3; end;
if(kd == 2|kd == 3)l=1; end;
if(kd == 4)jf=ic-1; end;
t1=0.0d0;
for  j=mj:-1:ic+1;
t1=-q.*q./((2.0d0.*j+l).^2-b+t1);
end;  j=ic+1-1;
if(m <= 2);
t2=0.0d0;
if(kd == 1&m == 0);
t1=t1+t1;
end;
if(kd == 1&m == 2);
t1=-2.0.*q.*q./(4.0-b+t1)-4.0;
end;
if(kd == 2&m == 1)t1=t1+q; end;
if(kd == 3&m == 1)t1=t1-q; end;
else;
if(kd == 1)t0=4.0d0-b+2.0d0.*q.*q./b; end;
if(kd == 2)t0=1.0d0-b+q; end;
if(kd == 3)t0=1.0d0-b-q; end;
if(kd == 4)t0=4.0d0-b; end;
t2=-q.*q./t0;
for  j=j0:jf;
t2=-q.*q./((2.0d0.*j-l-l0).^2-b+t2);
end;  j=jf+1;
end;
f=(2.0d0.*ic+l).^2+t1+t2-b;
return;
end
function [kd,m,q,a0]=cv0(kd,m,q,a0,varargin);
%     =====================================================
%     Purpose: Compute the initial characteristic value of
%     Mathieu functions for m ó 12  or q ó 300 or
%     q ò m*m
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     Output:  A0 --- Characteristic value
%     Routines called:
%(1)CVQM for computing initial characteristic
%     value for q ó 3*m
%(2)CVQL for computing initial characteristic
%     value for q ò m*m
%     ====================================================
q2=q.*q;
if(m == 0);
if(q <= 1.0);
a0=(((.0036392.*q2-.0125868).*q2+.0546875).*q2-.5).*q2;
elseif(q <= 10.0);
a0=((3.999267d-3.*q-9.638957d-2).*q-.88297).*q +.5542818;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 1);
if(q <= 1.0&kd == 2);
a0=(((-6.51e-4.*q-.015625).*q-.125).*q+1.0).*q+1.0;
elseif(q <= 1.0&kd == 3);
a0=(((-6.51e-4.*q+.015625).*q-.125).*q-1.0).*q+1.0;
elseif(q <= 10.0& kd == 2);
a0=(((-4.94603d-4.*q+1.92917d-2).*q-.3089229).*q+1.33372).*q+.811752;
elseif(q <= 10.0&kd == 3);
a0=((1.971096d-3.*q-5.482465d-2).*q-1.152218).*q+1.10427;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 2);
if(q <= 1.0&kd == 1);
a0=(((-.0036391.*q2+.0125888).*q2-.0551939).*q2 +.416667).*q2+4.0;
elseif(q <= 1.0&kd == 4);
a0=(.0003617.*q2-.0833333).*q2+4.0;
elseif(q <= 15&kd == 1);
a0=(((3.200972d-4.*q-8.667445d-3).*q -1.829032d-4).*q+.9919999).*q+3.3290504;
elseif(q <= 10.0&kd == 4);
a0=((2.38446d-3.*q-.08725329).*q-4.732542d-3).*q+4.00909;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 3);
if(q <= 1.0&kd == 2);
a0=((6.348e-4.*q+.015625).*q+.0625).*q2+9.0;
elseif(q <= 1.0&kd == 3);
a0=((6.348e-4.*q-.015625).*q+.0625).*q2+9.0;
elseif(q <= 20.0&kd == 2);
a0=(((3.035731d-4.*q-1.453021d-2).*q +.19069602).*q-.1039356).*q+8.9449274;
elseif(q <= 15.0&kd == 3);
a0=((9.369364d-5.*q-.03569325).*q+.2689874).*q +8.771735;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 4);
if(q <= 1.0&kd == 1);
a0=((-2.1e-6.*q2+5.012e-4).*q2+.0333333).*q2+16.0;
elseif(q <= 1.0&kd == 4);
a0=((3.7e-6.*q2-3.669e-4).*q2+.0333333).*q2+16.0;
elseif(q <= 25.0&kd == 1);
a0=(((1.076676d-4.*q-7.9684875d-3).*q +.17344854).*q-.5924058).*q+16.620847;
elseif(q <= 20.0&kd == 4);
a0=((-7.08719d-4.*q+3.8216144d-3).*q +.1907493).*q+15.744;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 5);
if(q <= 1.0&kd == 2);
a0=((6.8e-6.*q+1.42e-5).*q2+.0208333).*q2+25.0;
elseif(q <= 1.0&kd == 3);
a0=((-6.8e-6.*q+1.42e-5).*q2+.0208333).*q2+25.0;
elseif(q <= 35.0&kd == 2);
a0=(((2.238231d-5.*q-2.983416d-3).*q +.10706975).*q-.600205).*q+25.93515;
elseif(q <= 25.0&kd == 3);
a0=((-7.425364d-4.*q+2.18225d-2).*q +4.16399d-2).*q+24.897;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 6);
if(q <= 1.0);
a0=(.4d-6.*q2+.0142857).*q2+36.0;
elseif(q <= 40.0&kd == 1);
a0=(((-1.66846d-5.*q+4.80263d-4).*q +2.53998d-2).*q-.181233).*q+36.423;
elseif(q <= 35.0&kd == 4);
a0=((-4.57146d-4.*q+2.16609d-2).*q-2.349616d-2).*q +35.99251;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m == 7);
if(q <= 10.0);
[m,q,a0]=cvqm(fix(m),q,a0);
elseif(q <= 50.0&kd == 2);
a0=(((-1.411114d-5.*q+9.730514d-4).*q -3.097887d-3).*q+3.533597d-2).*q+49.0547;
elseif(q <= 40.0&kd == 3);
a0=((-3.043872d-4.*q+2.05511d-2).*q -9.16292d-2).*q+49.19035;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m >= 8);
if(q <= 3..*m);
[m,q,a0]=cvqm(fix(m),q,a0);
elseif(q > m.*m);
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
else;
if(m == 8&kd == 1);
a0=(((8.634308d-6.*q-2.100289d-3).*q+.169072).*q -4.64336).*q+109.4211;
elseif(m == 8&kd == 4);
a0=((-6.7842d-5.*q+2.2057d-3).*q+.48296).*q+56.59;
elseif(m == 9&kd == 2);
a0=(((2.906435d-6.*q-1.019893d-3).*q+.1101965).*q -3.821851).*q+127.6098;
elseif(m == 9&kd == 3);
a0=((-9.577289d-5.*q+.01043839).*q+.06588934).*q +78.0198;
elseif(m == 10&kd == 1);
a0=(((5.44927d-7.*q-3.926119d-4).*q+.0612099).*q -2.600805).*q+138.1923;
elseif(m == 10&kd == 4);
a0=((-7.660143d-5.*q+.01132506).*q-.09746023).*q +99.29494;
elseif(m == 11&kd == 2);
a0=(((-5.67615d-7.*q+7.152722d-6).*q+.01920291).*q -1.081583).*q+140.88;
elseif(m == 11&kd == 3);
a0=((-6.310551d-5.*q+.0119247).*q-.2681195).*q +123.667;
elseif(m == 12&kd == 1);
a0=(((-2.38351d-7.*q-2.90139d-5).*q+.02023088).*q -1.289).*q+171.2723;
elseif(m == 12&kd == 4);
a0=(((3.08902d-7.*q-1.577869d-4).*q+.0247911).*q -1.05454).*q+161.471;
end;
end;
end;
return;
end
function [kd,m,q,a0]=cvql(kd,m,q,a0,varargin);
%     ========================================================
%     Purpose: Compute the characteristic value of Mathieu
%     functions  for q ò 3m
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     Output:  A0 --- Initial characteristic value
%     ========================================================
if(kd == 1|kd == 2)w=2.0d0.*m+1.0d0; end;
if(kd == 3|kd == 4)w=2.0d0.*m-1.0d0; end;
w2=w.*w;
w3=w.*w2;
w4=w2.*w2;
w6=w2.*w4;
d1=5.0+34.0./w2+9.0./w4;
d2=(33.0+410.0./w2+405.0./w4)./w;
d3=(63.0+1260.0./w2+2943.0./w4+486.0./w6)./w2;
d4=(527.0+15617.0./w2+69001.0./w4+41607.0./w6)./w3;
c1=128.0;
p2=q./w4;
p1=sqrt(p2);
cv1=-2.0.*q+2.0.*w.*sqrt(q)-(w2+1.0)./8.0;
cv2=(w+3.0./w)+d1./(32.0.*p1)+d2./(8.0.*c1.*p2);
cv2=cv2+d3./(64.0.*c1.*p1.*p2)+d4./(16.0.*c1.*c1.*p2.*p2);
a0=cv1-cv2./(c1.*p1);
return;
end
function [m,q,a0]=cvqm(m,q,a0,varargin);
%     =====================================================
%     Purpose: Compute the characteristic value of Mathieu
%     functions for q ó m*m
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     Output:  A0 --- Initial characteristic value
%     =====================================================
hm1=.5.*q./(fix(m).*fix(m)-1.0);
hm3=.25.*hm1.^3./(fix(m).*fix(m)-4.0);
hm5=hm1.*hm3.*q./((fix(m).*fix(m)-1.0).*(fix(m).*fix(m)-9.0));
a0=fix(m).*fix(m)+q.*(hm1+(5.0.*fix(m).*fix(m)+7.0).*hm3+(9.0.*fix(m).^4+58.0.*fix(m).*fix(m)+29.0).*hm5);
return;
end

