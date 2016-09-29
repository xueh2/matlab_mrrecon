function [w,f,r2] = fitMagnLS_1R2_StartInPhase( images, t, fF,winit,finit,r2init,THRESHOLD )

if nargin<7
  THRESHOLD = 0;
end

% Initialize variables
sx = size(images,1);
sy = size(images,2);
N = size(images,3);

if nargin < 4
    r2 = zeros(sx,sy);
    w = zeros(sx,sy);
    f = zeros(sx,sy);
else
   r2 = r2init;
   w = winit;
   f = finit;
end
    
smagn = sqrt(sum(sum(abs(images).^2,3),4));
smagn = smagn/max(smagn(:));

% Loop over voxels
for kx=1:sx
%  disp([num2str(kx)]);
    for ky=1:sy
% $$$       disp([num2str([kx ky])]);
      scur = reshape(images(kx,ky,:),[],1);
      
      if norm(smagn(kx,ky))<THRESHOLD
        w(kx,ky) = 0;
        f(kx,ky) = 0;
        r2(kx,ky) = 0;
      else
        
        % Initial guesstimate
        if nargin <4
          w0 = (min(abs(scur)) + max(abs(scur)))/2;
          f0 = max(abs(scur)) - w0; 
          r20 = 10;
        else 
          w0 = abs(w(kx,ky));
          f0 = abs(f(kx,ky));
          r20 = r2(kx,ky);
        end
        x0 = [w0; f0;r20];       
        x = x0;


        % ML estimate
        if kx==1 & ky==1001
          DEBUG = 1;
        else
          DEBUG = 0;
        end
        
        MAXR2 = 300;
        MAXS = max(abs(scur(:)));
        lb = [-2*MAXS*exp(MAXR2*t(1))*ones(2,1);0];
        ub = [2*MAXS*exp(MAXR2*t(1))*ones(2,1);MAXR2];
        %[x,fval,exitflag,output] = fmincon(@(a)fitMagnLS_1R2_ricianLogLikelihood( a, t,fS, sigma, abs(scur),DEBUG ),x0,[],[],[],[],lb,ub);
        
        options = optimset('Display','off','Jacobian','on', ...
                           'DerivativeCheck','off','MaxIter', 500, ...
                           'MaxFunEvals', 1e10 , 'LargeScale', ...
                           'on','TolFun',1e-14,'TolX',1e-14);
%        x = lsqnonlin( @(y)fitError_Fieldmap_NoR2_MP(y,s,t,N,M,C,fs,relAmps,DEBUG),x0,lb,ub,options);
        
        [x,fval,exitflag] = lsqnonlin(@(a)fitMagnLS_1R2( a, t,fF, abs(scur),DEBUG ),x0,lb,ub,options) ;
        
% $$$           [r2w(kx,ky),indw] = min(x(3:4));
% $$$           [r2f(kx,ky),indf] = max(x(3:4));                
% $$$           w(kx,ky) = abs(x(indw));
% $$$           f(kx,ky) = abs(x(indf));

        r2(kx,ky) = x(3);  
        w(kx,ky) = x(1);
        f(kx,ky) = x(2);

      end
      
    end
end

        
function [fiter,J] = fitMagnLS_1R2( x, t,fS, s, DEBUG_ON )

if nargin<5
  DEBUG_ON = 0;
end

N = length(s);

rhoW = abs(x(1));
rhoF = abs(x(2));
r2 = x(3);


s = reshape(abs(s),[],1);
t = reshape(t,[],1);

% Signal magnitude
shat = abs(rhoW.*exp(-r2*t) + rhoF*exp(j*(2*pi*fS*t)).*exp(-r2*t));

fiter = shat(:) - s(:);

% $$$ disp(['fiter = ' num2str(norm(fiter))])
% $$$ disp(['r2 = ' num2str([r2])]);
% $$$ disp(['amps = ' num2str([rhoW rhoF])])


if DEBUG_ON>0
%  plot(s);hold on;plot(shat,'r');hold off;drawnow;pause(0.3);
  disp(['r2 = ' num2str([r2])]);
  disp(['amps = ' num2str([rhoW rhoF])])
  disp(['fiter = ' num2str(norm(fiter))])
end

if nargout>1
  J = zeros(3,N);

  ds1 = exp(-r2*t).*(rhoW.*exp(-r2*t) + rhoF*cos(2*pi*fS*t).*exp(-r2*t))./(shat);
  ds2 = (rhoW*cos(2*pi*fS*t).*exp(-r2*t).*exp(-r2*t) + rhoF.*exp(-2*r2*t))./(shat);
  ds3 = -rhoW*t.*exp(-r2*t).*(rhoW.*exp(-r2*t) + rhoF*cos(2*pi*fS*t).*exp(-r2*t))./(shat);
  ds4 = (-rhoF^2*t.*exp(-2*r2*t) - rhoW*rhoF*t.*exp(-(r2+r2)*t).*cos(2*pi*fS*t) )./(shat);
  
  J = [ds1(:),ds2(:), ds3(:)+ds4(:)];
    
end

