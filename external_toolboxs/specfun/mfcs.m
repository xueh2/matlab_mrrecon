function mfcs
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     =======================================================
%     Purpose: This program computes the Fresnel integrals
%     C(x)and S(x)using subroutine FCS
%     Input :  x --- Argument of C(x)and S(x)
%     Output:  C --- C(x)
%     S --- S(x)
%     Example:
%     x          C(x)S(x)
%     -----------------------------------
%     0.0      .00000000      .00000000
%     0.5      .49234423      .06473243
%     1.0      .77989340      .43825915
%     1.5      .44526118      .69750496
%     2.0      .48825341      .34341568
%     2.5      .45741301      .61918176
%     =======================================================
x=[];c=[];s=[];
fprintf(1,'%s \n','please enter x ');
%     READ(*,*)X
x=2.5;
fprintf(1,'%s \n','   x          c(x)s(x)');
fprintf(1,'%s \n',' -----------------------------------');
[x,c,s]=fcs(x,c,s);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%15.8g',1,2) ' \n'],x,c,s);
%format(1x,f5.1,2f15.8);
end
function [x,c,s]=fcs(x,c,s,varargin);
%     =================================================
%     Purpose: Compute Fresnel integrals C(x)and S(x)
%     Input :  x --- Argument of C(x)and S(x)
%     Output:  C --- C(x)
%     S --- S(x)
%     =================================================
eps=1.0d-15;
pi=3.141592653589793d0;
xa=abs(x);
px=pi.*xa;
t=.5d0.*px.*xa;
t2=t.*t;
if(xa == 0.0);
c=0.0d0;
s=0.0d0;
elseif(xa < 2.5d0);
r=xa;
c=r;
for  k=1:50;
r=-.5d0.*r.*(4.0d0.*k-3.0d0)./k./(2.0d0.*k-1.0d0)./(4.0d0.*k+1.0d0).*t2;
c=c+r;
if(abs(r)< abs(c).*eps)break; end;
end;
s=xa.*t./3.0d0;
r=s;
for  k=1:50;
r=-.5d0.*r.*(4.0d0.*k-1.0d0)./k./(2.0d0.*k+1.0d0)./(4.0d0.*k+3.0d0).*t2;
s=s+r;
if(abs(r)< abs(s).*eps)break; end;
end;
elseif(xa < 4.5d0);
m=fix(42.0+1.75.*t);
su=0.0d0;
c=0.0d0;
s=0.0d0;
f1=0.0d0;
f0=1.0d-100;
for  k=m:-1:0;
f=(2.0d0.*k+3.0d0).*f0./t-f1;
if(k == fix(k./2).*2);
c=c+f;
else;
s=s+f;
end;
su=su+(2.0d0.*k+1.0d0).*f.*f;
f1=f0;
f0=f;
end;  k=0-1;
q=sqrt(su);
c=c.*xa./q;
s=s.*xa./q;
else;
r=1.0d0;
f=1.0d0;
for  k=1:20;
r=-.25d0.*r.*(4.0d0.*k-1.0d0).*(4.0d0.*k-3.0d0)./t2;
f=f+r;
end;  k=20+1;
r=1.0d0./(px.*xa);
g=r;
for  k=1:12;
r=-.25d0.*r.*(4.0d0.*k+1.0d0).*(4.0d0.*k-1.0d0)./t2;
g=g+r;
end;  k=12+1;
t0=t-fix(t./(2.0d0.*pi)).*2.0d0.*pi;
c=.5d0+(f.*sin(t0)-g.*cos(t0))./px;
s=.5d0-(f.*cos(t0)+g.*sin(t0))./px;
end;
if(x < 0.0d0);
c=-c;
s=-s;
end;
return;
end

