function mfcszo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ===========================================================
%     Purpose : This program computes the complex zeros of the
%     Fresnel integral C(z)or S(z)using subroutine
%     FCSZO
%     Input :   KF  --- Function code
%     KF=1 for C(z)or KF=2 for S(z)
%     NT  --- Total number of zeros
%     Output:   ZO(L)--- L-th zero of C(z)or S(z)
%     Example:  NT=10
%     n     Complex zeros of C(z)Complex zeros of S(z)
%     ------------------------------------------------------------
%     1    1.7436675 + i .30573506      2.0092570 + i .28854790
%     2    2.6514596 + i .25290396      2.8334772 + i .24428524
%     3    3.3203593 + i .22395346      3.4675331 + i .21849268
%     4    3.8757345 + i .20474747      4.0025782 + i .20085103
%     5    4.3610635 + i .19066973      4.4741893 + i .18768859
%     6    4.7976077 + i .17970801      4.9006784 + i .17732036
%     7    5.1976532 + i .17081930      5.2929467 + i .16884418
%     8    5.5690602 + i .16339854      5.6581068 + i .16172492
%     9    5.9172173 + i .15706585      6.0011034 + i .15562108
%     10    6.2460098 + i .15156826      6.3255396 + i .15030246
%     ===========================================================
kf=[];nt=[];zo=[];
 zo=zeros(1,100);
fprintf(1,'%s \n','please enter kf and nt ');
%     READ(*,*)KF,NT
kf=2;
nt=10;
fprintf(1,[repmat(' ',1,2),'kf=','%2g',',     ','nt=','%3g' ' \n'],kf,nt);
fprintf(1,'%s \n',' *****     please wait %     *****');
[kf,nt,zo]=fcszo(kf,nt,zo);
fprintf(1,'%0.15g \n');
if(kf == 1)fprintf(1,'%s \n','  n        complex zeros of c(z)'); end;
if(kf == 2)fprintf(1,'%s \n','  n        complex zeros of s(z)'); end;
fprintf(1,'%s \n','-----------------------------------');
for  i=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g','%13.8g',repmat(' ',1,2),'+i','%13.8g' ' \n'],i,zo(i));
end;  i=nt+1;
%format(2x,',i2,',     ',',i3);
%format(1x,i3,f13.8,2x,2h+i,f13.8);
end
function [kf,nt,zo]=fcszo(kf,nt,zo,varargin);
%     ===============================================================
%     Purpose: Compute the complex zeros of Fresnel integral C(z)
%     or S(z)using modified Newton's iteration method
%     Input :  KF  --- Function code
%     KF=1 for C(z)or KF=2 for S(z)
%     NT  --- Total number of zeros
%     Output:  ZO(L)--- L-th zero of C(z)or S(z)
%     Routines called:
%(1)CFC for computing Fresnel integral C(z)
%(2)CFS for computing Fresnel integral S(z)
%     ==============================================================
z=[];zf=[];zd=[];
w=0.0;
pi=3.141592653589793d0;
for  nr=1:nt;
if(kf == 1)psq=sqrt(4.0d0.*nr-1.0d0); end;
if(kf == 2)psq=2.0d0.*nr.^(0.5); end;
px=psq-log(pi.*psq)./(pi.*pi.*psq.^3.0);
py=log(pi.*psq)./(pi.*psq);
z=complex(px,py);
if(kf == 2);
if(nr == 2)z=complex(2.8334,0.2443); end;
if(nr == 3)z=complex(3.4674,0.2185); end;
if(nr == 4)z=complex(4.0025,0.2008); end;
end;
it=0;
while (1);
it=it+1;
if(kf == 1)[z,zf,zd]=cfc(z,zf,zd); end;
if(kf == 2)[z,zf,zd]=cfs(z,zf,zd); end;
zp=complex(1.0d0,0.0d0);
for  i=1:nr-1;
zp=zp.*(z-zo(i));
end;  i=nr-1+1;
zfd=zf./zp;
zq=complex(0.0d0,0.0d0);
for  i=1:nr-1;
zw=complex(1.0d0,0.0d0);
for  j=1:nr-1;
if(~(j == i))zw=zw.*(z-zo(j)); end;
end;  j=nr-1+1;
zq=zq+zw;
end;  i=nr-1+1;
zgd=(zd-zq.*zfd)./zp;
z=z-zfd./zgd;
w0=w;
w=abs(z);
if(~(it <= 50&abs((w-w0)./w)> 1.0d-12))break; end;
end;
zo(nr)=z;
end;
return;
end
function [z,zf,zd]=cfc(z,zf,zd,varargin);
%     =========================================================
%     Purpose: Compute complex Fresnel integral C(z)and C'(z)
%     Input :  z --- Argument of C(z)
%     Output:  ZF --- C(z)
%     ZD --- C'(z)
%     =========================================================
eps=1.0d-14;
pi=3.141592653589793d0;
w0=abs(z);
zp=0.5d0.*pi.*z.*z;
zp2=zp.*zp;
z0=complex(0.0d0,0.0d0);
if(z == z0);
c=z0;
elseif(w0 <= 2.5);
cr=z;
c=cr;
igoon=1;
for k=1:80;
cr=-.5d0.*cr.*(4.0d0.*k-3.0d0)./k./(2.0d0.*k-1.0d0)./(4.0d0.*k+1.0d0).*zp2;
c=c+cr;
wa=abs(c);
if(abs((wa-wa0)./wa)< eps&k > 10)igoon=0; end;
if(igoon == 1);
wa0=wa;
else;
break;
end;
end;
elseif(w0 > 2.5&w0 < 4.5);
m=85;
c=z0;
cf1=z0;
cf0=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=(2.0d0.*k+3.0d0).*cf0./zp-cf1;
if(k == fix(k./2).*2)c=c+cf; end;
cf1=cf0;
cf0=cf;
end;  k=0-1;
c=sqrt(2.0d0./(pi.*zp)).*sin(zp)./cf.*c;
else;
cr=complex(1.0d0,0.0d0);
cf=complex(1.0d0,0.0d0);
for  k=1:20;
cr=-.25d0.*cr.*(4.0d0.*k-1.0d0).*(4.0d0.*k-3.0d0)./zp2;
cf=cf+cr;
end;  k=20+1;
cr=1.0d0./(pi.*z.*z);
cg=cr;
for  k=1:12;
cr=-.25d0.*cr.*(4.0d0.*k+1.0d0).*(4.0d0.*k-1.0d0)./zp2;
cg=cg+cr;
end;  k=12+1;
c=.5d0+(cf.*sin(zp)-cg.*cos(zp))./(pi.*z);
end;
zf=c;
zd=cos(0.5.*pi.*z.*z);
return;
end
function [z,zf,zd]=cfs(z,zf,zd,varargin);
%     =========================================================
%     Purpose: Compute complex Fresnel Integral S(z)and S'(z)
%     Input :  z  --- Argument of S(z)
%     Output:  ZF --- S(z)
%     ZD --- S'(z)
%     =========================================================
wb0=0.0;
eps=1.0d-14;
pi=3.141592653589793d0;
w0=abs(z);
zp=0.5d0.*pi.*z.*z;
zp2=zp.*zp;
z0=complex(0.0d0,0.0d0);
if(z == z0);
s=z0;
elseif(w0 <= 2.5);
s=z.*zp./3.0d0;
cr=s;
igoon=1;
for k=1:80;
cr=-.5d0.*cr.*(4.0d0.*k-1.0d0)./k./(2.0d0.*k+1.0d0)./(4.0d0.*k+3.0d0).*zp2;
s=s+cr;
wb=abs(s);
if(abs(wb-wb0)< eps&k > 10)igoon=1; end;
if(igoon == 1);
wb0=wb;
else;
break;
end;
end;
elseif(w0 > 2.5&w0 < 4.5);
m=85;
s=z0;
cf1=z0;
cf0=complex(1.0d-100,0.0d0);
for  k=m:-1:0;
cf=(2.0d0.*k+3.0d0).*cf0./zp-cf1;
if(k ~= fix(k./2).*2)s=s+cf; end;
cf1=cf0;
cf0=cf;
end;  k=0-1;
s=sqrt(2.0d0./(pi.*zp)).*sin(zp)./cf.*s;
else;
cr=complex(1.0d0,0.0d0);
cf=complex(1.0d0,0.0d0);
for  k=1:20;
cr=-.25d0.*cr.*(4.0d0.*k-1.0d0).*(4.0d0.*k-3.0d0)./zp2;
cf=cf+cr;
end;  k=20+1;
cr=1.0d0./(pi.*z.*z);
cg=cr;
for  k=1:12;
cr=-.25d0.*cr.*(4.0d0.*k+1.0d0).*(4.0d0.*k-1.0d0)./zp2;
cg=cg+cr;
end;  k=12+1;
s=.5d0-(cf.*cos(zp)+cg.*sin(zp))./(pi.*z);
end;
zf=s;
zd=sin(0.5.*pi.*z.*z);
return;
end

