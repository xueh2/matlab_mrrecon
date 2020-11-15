function mklvnb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes Kelvin functions ber x,
%                bei x, ker x and kei x, and their derivatives
%                using subroutine KLVNB
%       Input :  x   --- Argument of Kelvin functions
%       Output:  BER --- ber x
%                BEI --- bei x
%                GER --- ker x
%                GEI --- kei x
%                DER --- ber'x
%                DEI --- bei'x
%                HER --- ker'x
%                HEI --- kei'x
%       Example:
%     x       ber x         bei x         ker x         kei x
%   -------------------------------------------------------------
%     0    .100000D+01   .000000D+00    ì           -.785398D+00
%     5   -.623008D+01   .116034D+00  -.115117D-01   .111876D-01
%    10    .138840D+03   .563705D+02   .129466D-03  -.307525D-03
%    15   -.296725D+04  -.295271D+04  -.151433D-07   .796289D-05
%    20    .474894D+05   .114775D+06  -.771523D-07  -.185894D-06
%     x       ber'x         bei'x         ker'x         kei'x
%   -------------------------------------------------------------
%     0    .000000D+00   .000000D+00  - ì            .000000D+00
%     5   -.384534D+01  -.435414D+01   .171934D-01  -.819979D-03
%    10    .511952D+02   .135309D+03  -.315597D-03   .140914D-03
%    15    .910555D+02  -.408776D+04   .564468D-05  -.588222D-05
%    20   -.488032D+05   .111855D+06  -.750186D-07   .190624D-06
%       ========================================================
x=[];ber=[];bei=[];ger=[];gei=[];der=[];dei=[];her=[];hei=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=20;
fprintf(1,'%s ','   x        ber x           bei x');fprintf(1,'%s \n','           ker x           kei x');
fprintf(1,'%s ','--------------------------------');fprintf(1,'%s \n','--------------------------------------');
[x,ber,bei,ger,gei,der,dei,her,hei]=klvnb(x,ber,bei,ger,gei,der,dei,her,hei);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.6g',1,4) ' \n'],x,ber,bei,ger,gei);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   x        ber''x           bei''x');fprintf(1,'%s \n','           ker''x           kei''x');
fprintf(1,'%s ','--------------------------------');fprintf(1,'%s \n','--------------------------------------');
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.6g',1,4) ' \n'],x,der,dei,her,hei);
%format(1x,f5.1,4d16.6);
end
function [x,ber,bei,ger,gei,der,dei,her,hei]=klvnb(x,ber,bei,ger,gei,der,dei,her,hei,varargin);
%       ======================================================
%       Purpose: Compute Kelvin functions ber x, bei x, ker x
%                and kei x, and their derivatives(x > 0)
%       Input :  x   --- Argument of Kelvin functions
%       Output:  BER --- ber x
%                BEI --- bei x
%                GER --- ker x
%                GEI --- kei x
%                DER --- ber'x
%                DEI --- bei'x
%                HER --- ker'x
%                HEI --- kei'x
%       ================================================
pi=3.141592653589793d0;
if(x == 0.0d0);
ber=1.0d0;
bei=0.0d0;
ger=1.0d+300;
gei=-.25d0.*pi;
der=0.0d0;
dei=0.0d0;
her=-1.0d+300;
hei=0.0d0;
elseif(x < 8.0d0);
t=x./8.0d0;
t2=t.*t;
u=t2.*t2;
ber=((((((-.901d-5.*u+.122552d-2).*u-.08349609d0).*u+2.64191397d0).*u-32.36345652d0).*u +113.77777774d0).*u-64.0d0).*u+1.0d0;
bei=t.*t.*((((((.11346d-3.*u-.01103667d0).*u +.52185615d0).*u-10.56765779d0).*u+72.81777742d0).*u-113.77777774d0).*u+16.0d0);
ger=((((((-.2458d-4.*u+.309699d-2).*u-.19636347d0).*u+5.65539121d0).*u-60.60977451d0).*u+171.36272133d0).*u-59.05819744d0).*u-.57721566d0;
ger=ger-log(.5d0.*x).*ber+.25d0.*pi.*bei;
gei=t2.*((((((.29532d-3.*u-.02695875d0).*u +1.17509064d0).*u-21.30060904d0).*u+124.2356965d0).*u-142.91827687d0).*u +6.76454936d0);
gei=gei-log(.5d0.*x).*bei-.25d0.*pi.*ber;
der=x.*t2.*((((((-.394d-5.*u+.45957d-3).*u-.02609253d0).*u+.66047849d0).*u-6.0681481d0).*u +14.22222222d0).*u-4.0d0);
dei=x.*((((((.4609d-4.*u-.379386d-2).*u+.14677204d0).*u-2.31167514d0).*u+11.37777772d0).*u -10.66666666d0).*u+.5d0);
her=x.*t2.*((((((-.1075d-4.*u+.116137d-2).*u-.06136358d0).*u+1.4138478d0).*u-11.36433272d0).*u+21.42034017d0).*u-3.69113734d0);
her=her-log(.5d0.*x).*der-ber./x+.25d0.*pi.*dei;
hei=x.*((((((.11997d-3.*u-.926707d-2).*u+.33049424d0).*u-4.65950823d0).*u+19.41182758d0).*u-13.39858846d0).*u+.21139217d0);
hei=hei-log(.5d0.*x).*dei-bei./x-.25d0.*pi.*der;
else;
t=8.0d0./x;
for  l=1:2;
v=(-1).^l.*t;
tpr=((((.6d-6.*v-.34d-5).*v-.252d-4).*v-.906d-4).*v.*v+.0110486d0).*v;
tpi=((((.19d-5.*v+.51d-5).*v.*v-.901d-4).*v-.9765d-3).*v-.0110485d0).*v-.3926991d0;
if(l == 1);
tnr=tpr;
tni=tpi;
end;
end;  l=2+1;
yd=x./sqrt(2.0d0);
ye1=exp(yd+tpr);
ye2=exp(-yd+tnr);
yc1=1.0d0./sqrt(2.0d0.*pi.*x);
yc2=sqrt(pi./(2.0d0.*x));
csp=cos(yd+tpi);
ssp=sin(yd+tpi);
csn=cos(-yd+tni);
ssn=sin(-yd+tni);
ger=yc2.*ye2.*csn;
gei=yc2.*ye2.*ssn;
fxr=yc1.*ye1.*csp;
fxi=yc1.*ye1.*ssp;
ber=fxr-gei./pi;
bei=fxi+ger./pi;
for  l=1:2;
v=(-1).^l.*t;
ppr=(((((.16d-5.*v+.117d-4).*v+.346d-4).*v+.5d-6).*v-.13813d-2).*v-.0625001d0).*v+.7071068d0;
ppi=(((((-.32d-5.*v-.24d-5).*v+.338d-4).*v+.2452d-3).*v+.13811d-2).*v-.1d-6).*v+.7071068d0;
if(l == 1);
pnr=ppr;
pni=ppi;
end;
end;  l=2+1;
her=gei.*pni-ger.*pnr;
hei=-(gei.*pnr+ger.*pni);
der=fxr.*ppr-fxi.*ppi-hei./pi;
dei=fxi.*ppr+fxr.*ppi+her./pi;
end;
return;
end

