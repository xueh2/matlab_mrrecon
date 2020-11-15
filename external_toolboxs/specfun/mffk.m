function mffk
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ==============================================================
%     Purpose: This program computes the modified Fresnel integrals
%     Fñ(x)and Kñ(x)using subroutine FFK
%     Input :  x   --- Argument of Fñ(x)and Kñ(x)
%     KS  --- Sign code
%     KS=0 for calculating F+(x)and K+(x)
%     KS=1 for calculating F_(x)and K_(x)
%     Output:  FR  --- Re[Fñ(x)]
%     FI  --- Im[Fñ(x)]
%     FM  --- |Fñ(x)|
%     FA  --- Arg[Fñ(x)](Degs.)
%     GR  --- Re[Kñ(x)]
%     GI  --- Im[Kñ(x)]
%     GM  --- |Kñ(x)|
%     GA  --- Arg[Kñ(x)](Degs.)
%     Example:
%     x     Re[Fñ(x)]ñIm[Fñ(x)]Mod[Fñ(x)]ñArg[Fñ(x)]
%     ----------------------------------------------------------
%     0.0    .62665707    .62665707    .88622693    45.000000
%     2.0    .16519561   -.17811942    .24293233   -47.155835
%     4.0    .03219674   -.12047678    .12470479   -75.037684
%     6.0    .08245304   -.01180212    .08329342    -8.145843
%     8.0   -.05729996    .02493542    .06249048   156.482601
%     10.0    .02553188    .04298617    .04999688    59.291561
%     x     Re[Kñ(x)]ñIm[Kñ(x)]Mod[Kñ(x)]ñArg[Kñ(x)]
%     ----------------------------------------------------------
%     0.0    .50000000    .00000000    .50000000     0.000000
%     2.0    .10702394    .08562295    .13705989    38.661047
%     4.0    .05126306    .04818949    .07035714    43.229843
%     6.0    .03368650    .03276566    .04699328    44.206095
%     8.0    .02512396    .02473472    .03525648    44.552712
%     10.0    .02004532    .01984592    .02820772    44.713609
%     ===============================================================
x=[];fr=[];fi=[];fm=[];fa=[];gr=[];gi=[];gm=[];ga=[];
fprintf(1,'%s \n','please enter x');
%     READ(*,*)X
x=10.0;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   x      re[fñ(x)]ñim[fñ(x)]');fprintf(1,'%s \n', 'mod[fñ(x)]ñarg[fñ(x)]');
fprintf(1,'%s ',' ---------------------------------------');fprintf(1,'%s \n','-----------------------');
[dumvar1,x,fr,fi,fm,fa,gr,gi,gm,ga]=ffk(0,x,fr,fi,fm,fa,gr,gi,gm,ga);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%14.8g',1,3),'%14.6g' ' \n'],x,fr,fi,fm,fa);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   x      re[kñ(x)]ñim[kñ(x)]');fprintf(1,'%s \n', 'mod[kñ(x)]ñarg[kñ(x)]');
fprintf(1,'%s ',' ---------------------------------------');fprintf(1,'%s \n','-----------------------');
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%14.8g',1,3),'%14.6g' ' \n'],x,gr,gi,gm,ga);
%format(1x,f5.1,3f14.8,f14.6);
end
function [ks,x,fr,fi,fm,fa,gr,gi,gm,ga]=ffk(ks,x,fr,fi,fm,fa,gr,gi,gm,ga,varargin);
%     =======================================================
%     Purpose: Compute modified Fresnel integrals Fñ(x)
%     and Kñ(x)
%     Input :  x   --- Argument of Fñ(x)and Kñ(x)
%     KS  --- Sign code
%     KS=0 for calculating F+(x)and K+(x)
%     KS=1 for calculating F_(x)and K_(x)
%     Output:  FR  --- Re[Fñ(x)]
%     FI  --- Im[Fñ(x)]
%     FM  --- |Fñ(x)|
%     FA  --- Arg[Fñ(x)](Degs.)
%     GR  --- Re[Kñ(x)]
%     GI  --- Im[Kñ(x)]
%     GM  --- |Kñ(x)|
%     GA  --- Arg[Kñ(x)](Degs.)
%     ======================================================
srd= 57.29577951308233d0;
eps=1.0d-15;
pi=3.141592653589793d0;
pp2=1.2533141373155d0;
p2p=.7978845608028654d0;
xa=abs(x);
x2=x.*x;
x4=x2.*x2;
if(x == 0.0d0);
fr=.5d0.*sqrt(0.5d0.*pi);
fi=(-1).^fix(ks).*fr;
fm=sqrt(0.25d0.*pi);
fa=(-1).^fix(ks).*45.0d0;
gr=.5d0;
gi=0.0d0;
gm=.5d0;
ga=0.0d0;
else;
if(xa <= 2.5d0);
xr=p2p.*xa;
c1=xr;
for  k=1:50;
xr=-.5d0.*xr.*(4.0d0.*k-3.0d0)./k./(2.0d0.*k-1.0d0)./(4.0d0.*k+1.0d0).*x4;
c1=c1+xr;
if(abs(xr./c1)< eps)break; end;
end;
s1=p2p.*xa.*xa.*xa./3.0d0;
xr=s1;
for  k=1:50;
xr=-.5d0.*xr.*(4.0d0.*k-1.0d0)./k./(2.0d0.*k+1.0d0)./(4.0d0.*k+3.0d0).*x4;
s1=s1+xr;
if(abs(xr./s1)< eps)break; end;
end;
elseif(xa < 5.5d0);
m=fix(42+1.75.*x2);
xsu=0.0d0;
xc=0.0d0;
xs=0.0d0;
xf1=0.0d0;
xf0=1d-100;
for  k=m:-1:0;
xf=(2.0d0.*k+3.0d0).*xf0./x2-xf1;
if(k == 2.*fix(k./2));
xc=xc+xf;
else;
xs=xs+xf;
end;
xsu=xsu+(2.0d0.*k+1.0d0).*xf.*xf;
xf1=xf0;
xf0=xf;
end;  k=0-1;
xq=sqrt(xsu);
xw=p2p.*xa./xq;
c1=xc.*xw;
s1=xs.*xw;
else;
xr=1.0d0;
xf=1.0d0;
for  k=1:12;
xr=-.25d0.*xr.*(4.0d0.*k-1.0d0).*(4.0d0.*k-3.0d0)./x4;
xf=xf+xr;
end;  k=12+1;
xr=1.0d0./(2.0d0.*xa.*xa);
xg=xr;
for  k=1:12;
xr=-.25d0.*xr.*(4.0d0.*k+1.0d0).*(4.0d0.*k-1.0d0)./x4;
xg=xg+xr;
end;  k=12+1;
c1=.5d0+(xf.*sin(x2)-xg.*cos(x2))./sqrt(2.0d0.*pi)./xa;
s1=.5d0-(xf.*cos(x2)+xg.*sin(x2))./sqrt(2.0d0.*pi)./xa;
end;
fr=pp2.*(.5d0-c1);
fi0=pp2.*(.5d0-s1);
fi=(-1).^fix(ks).*fi0;
fm=sqrt(fr.*fr+fi.*fi);
if(fr >= 0.0);
fa=srd.*atan(fi./fr);
elseif(fi > 0.0);
fa=srd.*(atan(fi./fr)+pi);
elseif(fi < 0.0);
fa=srd.*(atan(fi./fr)-pi);
end;
xp=x.*x+pi./4.0d0;
cs=cos(xp);
ss=sin(xp);
xq2=1.0d0./sqrt(pi);
gr=xq2.*(fr.*cs+fi0.*ss);
gi=(-1).^fix(ks).*xq2.*(fi0.*cs-fr.*ss);
gm=sqrt(gr.*gr+gi.*gi);
if(gr >= 0.0);
ga=srd.*atan(gi./gr);
elseif(gi > 0.0);
ga=srd.*(atan(gi./gr)+pi);
elseif(gi < 0.0);
ga=srd.*(atan(gi./gr)-pi);
end;
if(x < 0.0d0);
fr=pp2-fr;
fi=(-1).^fix(ks).*pp2-fi;
fm=sqrt(fr.*fr+fi.*fi);
fa=srd.*atan(fi./fr);
gr=cos(x.*x)-gr;
gi=-(-1).^fix(ks).*sin(x.*x)-gi;
gm=sqrt(gr.*gr+gi.*gi);
ga=srd.*atan(gi./gr);
end;
end;
return;
end

