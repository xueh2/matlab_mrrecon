function mitairy
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===========================================================
%       Purpose: This program computes the integrals of Airy
%                functions using subroutine ITAIRY
%       Input  : x   --- Upper limit of the integral
%       Output : APT --- Integration of Ai(t)from 0 and x
%                BPT --- Integration of Bi(t)from 0 and x
%                ANT --- Integration of Ai(-t)from 0 and x
%                BNT --- Integration of Bi(-t)from 0 and x
%       Example:
%         x      Ai(t)dt       Bi(t)dt       Ai(-t)dt     Bi(-t)dt
%        ----------------------------------------------------------
%         5    .33328759   .32147832D+03    .71788220    .15873094
%        10    .33333333   .14780980D+09    .76569840    .01504043
%        15    .33333333   .49673090D+16    .68358063    .07202621
%        20    .33333333   .47447423D+25    .71173925   -.03906173
%        25    .33333333   .78920820D+35    .70489539    .03293190
%       ===========================================================
x=[];apt=[];bpt=[];ant=[];bnt=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=5.0;
fprintf(1,[repmat(' ',1,3),'x',repmat(' ',1,8),'ai(t)dt',repmat(' ',1,7),'bi(t)dt',repmat(' ',1,9),'ai(-t)dt',repmat(' ',1,6),'bi(-t)dt' ' \n']);
fprintf(1,[repmat(' ',1,2),'----------------------------------','------------------------------' ' \n']);
[x,apt,bpt,ant,bnt]=itairy(x,apt,bpt,ant,bnt);
fprintf(1,[repmat(' ',1,1),'%5.1g','%14.8g',repmat(' ',1,2),'%15.8g',repmat('%14.8g',1,2) ' \n'],x,apt,bpt,ant,bnt);
%format(1x,f5.1,f14.8,2x,d15.8,2f14.8);
%format(3x,'x',8x,'ai(t)dt',7x,'bi(t)dt',9x, 'ai(-t)dt',6x,'bi(-t)dt');
%format(2x,'----------------------------------','------------------------------');
end
function [x,apt,bpt,ant,bnt]=itairy(x,apt,bpt,ant,bnt,varargin);
%       ======================================================
%       Purpose: Compute the integrals of Airy fnctions with
%                respect to t from 0 and x(x ò 0)
%       Input  : x   --- Upper limit of the integral
%       Output : APT --- Integration of Ai(t)from 0 and x
%                BPT --- Integration of Bi(t)from 0 and x
%                ANT --- Integration of Ai(-t)from 0 and x
%                BNT --- Integration of Bi(-t)from 0 and x
%       ======================================================
 a=zeros(1,16);
eps=1.0d-15;
pi=3.141592653589793d0;
c1=.355028053887817d0;
c2=.258819403792807d0;
sr3=1.732050807568877d0;
if(x == 0.0d0);
apt=0.0d0;
bpt=0.0d0;
ant=0.0d0;
bnt=0.0d0;
else;
if(abs(x)<= 9.25d0);
for  l=0:1;
x=(-1).^l.*x;
fx=x;
r=x;
for  k=1:40;
r=r.*(3.0.*k-2.0d0)./(3.0.*k+1.0d0).*x./(3.0.*k).*x./(3.0.*k-1.0d0).*x;
fx=fx+r;
if(abs(r)< abs(fx).*eps)break; end;
end;
gx=.5d0.*x.*x;
r=gx;
for  k=1:40;
r=r.*(3.0.*k-1.0d0)./(3.0.*k+2.0d0).*x./(3.0.*k).*x./(3.0.*k+1.0d0).*x;
gx=gx+r;
if(abs(r)< abs(gx).*eps)break; end;
end;
ant=c1.*fx-c2.*gx;
bnt=sr3.*(c1.*fx+c2.*gx);
if(l == 0);
apt=ant;
bpt=bnt;
else;
ant=-ant;
bnt=-bnt;
x=-x;
end;
end;
else;
a(:)=[.569444444444444d0,.891300154320988d0,.226624344493027d+01,.798950124766861d+01,.360688546785343d+02,.198670292131169d+03,.129223456582211d+04,.969483869669600d+04,.824184704952483d+05,.783031092490225d+06,.822210493622814d+07,.945557399360556d+08,.118195595640730d+10,.159564653040121d+11,.231369166433050d+12,.358622522796969d+13];
q2=1.414213562373095d0;
q0=.3333333333333333d0;
q1=.6666666666666667d0;
xe=x.*sqrt(x)./1.5d0;
xp6=1.0d0./sqrt(6.0d0.*pi.*xe);
su1=1.0d0;
r=1.0d0;
xr1=1.0d0./xe;
for  k=1:16;
r=-r.*xr1;
su1=su1+a(k).*r;
end;  k=16+1;
su2=1.0d0;
r=1.0d0;
for  k=1:16;
r=r.*xr1;
su2=su2+a(k).*r;
end;  k=16+1;
apt=q0-exp(-xe).*xp6.*su1;
bpt=2.0d0.*exp(xe).*xp6.*su2;
su3=1.0d0;
r=1.0d0;
xr2=1.0d0./(xe.*xe);
for  k=1:8;
r=-r.*xr2;
su3=su3+a(2.*k).*r;
end;  k=8+1;
su4=a(1).*xr1;
r=xr1;
for  k=1:7;
r=-r.*xr2;
su4=su4+a(2.*k+1).*r;
end;  k=7+1;
su5=su3+su4;
su6=su3-su4;
ant=q1-q2.*xp6.*(su5.*cos(xe)-su6.*sin(xe));
bnt=q2.*xp6.*(su5.*sin(xe)+su6.*cos(xe));
end;
end;
return;
end

