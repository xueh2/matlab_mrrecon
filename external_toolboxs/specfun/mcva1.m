function mcva1
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ============================================================
%     Purpose: This program computes a sequence of characteristic
%     values of Mathieu functions using subroutine CVA1
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,...)
%     KD=2 for cem(x,q,m = 1,3,5,...)
%     KD=3 for sem(x,q,m = 1,3,5,...)
%     KD=4 for sem(x,q,m = 2,4,6,...)
%     Output:  CV(I)--- Characteristic values; I = 1,2,3,...
%     For KD=1, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of cem for m = 0,2,4,...
%     For KD=2, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of cem for m = 1,3,5,...
%     For KD=3, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of sem for m = 1,3,5,...
%     For KD=4, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of sem for m = 0,2,4,...
%     Example: Mmax = 12,    q = 25.00
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
%     ============================================================
mmax=[];q=[];cv1=[];cv2=[];
 cv1=zeros(1,200);
cv2=zeros(1,200);
cve=zeros(1,200);
cvs=zeros(1,200);
fprintf(1,'%s \n','please enter mmax,q =?');
%     READ(*,*)MMAX,Q
mmax=12;
q=25.0;
fprintf(1,[repmat(' ',1,3),'mmax =','%3g',',    ','q =','%6.2g' ' \n'],mmax,q);
fprintf(1,'%0.15g \n');
[dumvar1,mmax,q,cv1]=cva1(1,mmax,q,cv1);
[dumvar1,mmax,q,cv2]=cva1(2,mmax,q,cv2);
for  j=1:mmax./2+1;
cve(2.*j-1)=cv1(j);
cve(2.*j)=cv2(j);
end;  j=mmax./2+1+1;
[dumvar1,mmax,q,cv1]=cva1(3,mmax,q,cv1);
[dumvar1,mmax,q,cv2]=cva1(4,mmax,q,cv2);
for  j=1:mmax./2+1;
cvs(2.*j)=cv1(j);
cvs(2.*j+1)=cv2(j);
end;  j=mmax./2+1+1;
fprintf(1,[repmat(' ',1,1),'characteristic values of mathieu functions' ' \n']);
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  m            a                  b');
fprintf(1,'%s \n','------------------------------------------');
for  j=0:mmax;
if(j == 0)fprintf(1,[repmat(' ',1,1),'%3g',repmat('%19.9g',1,2) ' \n'],j,cve(j+1)); end;
if(j ~= 0)fprintf(1,[repmat(' ',1,1),'%3g',repmat('%19.9g',1,2) ' \n'],j,cve(j+1),cvs(j+1)); end;
end;  j=mmax+1;
%format(3x,,i3,',    ',,f6.2);
%format(1x,i3,2f19.9);
%format(1x,'characteristic values of mathieu functions');
end
function [kd,m,q,cv]=cva1(kd,m,q,cv,varargin);
%     ============================================================
%     Purpose: Compute a sequence of characteristic values of
%     Mathieu functions
%     Input :  M  --- Maximum order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,תתת)
%     KD=2 for cem(x,q,m = 1,3,5,תתת)
%     KD=3 for sem(x,q,m = 1,3,5,תתת)
%     KD=4 for sem(x,q,m = 2,4,6,תתת)
%     Output:  CV(I)--- Characteristic values; I = 1,2,3,...
%     For KD=1, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of cem for m = 0,2,4,...
%     For KD=2, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of cem for m = 1,3,5,...
%     For KD=3, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of sem for m = 1,3,5,...
%     For KD=4, CV(1), CV(2), CV(3),..., correspond to
%     the characteristic values of sem for m = 0,2,4,...
%     ============================================================
 g=zeros(1,200);
h=zeros(1,200);
d=zeros(1,500);
e=zeros(1,500);
f=zeros(1,500);
eps=1.0d-14;
icm=fix(fix(m)./2)+1;
if(kd == 4)icm=m./2; end;
if(q == 0.0d0);
if(kd == 1);
for  ic=1:icm;
cv(ic)=4.0d0.*(ic-1.0d0).^2;
end;  ic=icm+1;
elseif(kd ~= 4);
for  ic=1:icm;
cv(ic)=(2.0d0.*ic-1.0d0).^2;
end;  ic=icm+1;
else;
for  ic=1:icm;
cv(ic)=4.0d0.*ic.*ic;
end;  ic=icm+1;
end;
else;
nm=fix(10+1.5.*fix(m)+0.5.*q);
e(1)=0.0d0;
f(1)=0.0d0;
if(kd == 1);
d(1)=0.0d0;
for  i=2:nm;
d(i)=4.0d0.*(i-1.0d0).^2;
e(i)=q;
f(i)=q.*q;
end;  i=nm+1;
e(2)=sqrt(2.0d0).*q;
f(2)=2.0d0.*q.*q;
elseif(kd ~= 4);
d(1)=1.0d0+(-1).^fix(kd).*q;
for  i=2:nm;
d(i)=(2.0d0.*i-1.0d0).^2;
e(i)=q;
f(i)=q.*q;
end;  i=nm+1;
else;
d(1)=4.0d0;
for  i=2:nm;
d(i)=4.0d0.*i.*i;
e(i)=q;
f(i)=q.*q;
end;  i=nm+1;
end;
xa=d(nm)+abs(e(nm));
xb=d(nm)-abs(e(nm));
nm1=nm-1;
for  i=1:nm1;
t=abs(e(i))+abs(e(i+1));
t1=d(i)+t;
if(xa < t1)xa=t1; end;
t1=d(i)-t;
if(t1 < xb)xb=t1; end;
end;  i=nm1+1;
for  i=1:icm;
g(i)=xa;
h(i)=xb;
end;  i=icm+1;
for  k=1:icm;
for  k1=k:icm;
if(g(k1)< g(k));
g(k)=g(k1);
break;
end;
end;
if(k ~= 1&h(k)< h(k-1))h(k)=h(k-1); end;
while (1);
x1=(g(k)+h(k))./2.0d0;
cv(k)=x1;
if(~(abs((g(k)-h(k))./x1)< eps));
j=0;
s=1.0d0;
for  i=1:nm;
if(s == 0.0d0)s=s+1.0d-30; end;
t=f(i)./s;
s=d(i)-t-x1;
if(s < 0.0)j=j+1; end;
end;  i=nm+1;
if(j < k);
h(k)=x1;
else;
g(k)=x1;
if(j >= icm);
g(icm)=x1;
else;
if(h(j+1)< x1)h(j+1)=x1; end;
if(x1 < g(j))g(j)=x1; end;
end;
end;
else;
cv(k)=x1;
break;
end;
end;
end;
end;
return;
end

