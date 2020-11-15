function [YRWI, YRWSTDH] = function_RWI3D(Z,V,EY,sigmanoise_arr,RWreg2,SObserved,SReconstructed);

% The function performs Regularized Wiener Inverse (RWI) reconstruction of 
% slices of a 3D object at specified positions.
%
% SYNTAX 
%   [YRWI, YRWSTDH] = function_RWI(Z,V,EY,RWreg2,SObserved,SReconstructed)
%
% DESCRIPTION
%   The function performs Regularized Wiener Inverse (RWI) reconstruction 
%   of a 3D object Y from its blurred and noisy observations Z obtained at 
%   positions specified at SObserved. The reconstructed is done for slices
%   which are specified to be at locations SReconstructed. Wiener filter
%   uses the estimate EY of a true signal.
% 
%   Z - is an set of blurred and noisy observations of a true object.
%
%   V - is a matrix of PSFs (Vij) which correspond to the observation of the 
%   object slice j from focusing at the position i.
%
%   EY - is an estimate of a true signal.
%
%   RWreg2 - is a regularization parameter r^2.
%
%   SObserved - is a set of values i, which specify the positions from
%   which observations Z are done.
%
%   SReconstructed - is a set of values which specify positions of slices 
%   which should be reconstructed. Observed and reconstructed slices are
%   not necessary should be the same.
%
% RETURNS
%   YRWI - is a reconstructed object
%
%   YRWSTDH - is a standard deviation for each slice estimate.
% 
% For more details read 9.6 '3D Inverse' part of the book.
%
% Dmitriy Paliy, Tampere University of Technology,
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

if nargin<1,
disp('The function performs Regularized Wiener Inverse (RWI) reconstruction of ');
disp('slices of a 3D object at specified positions.');
disp(' ');
disp('SYNTAX ');
disp('   [YRWI, YRWSTDH] = function_RWI(Z,V,EY,RWreg2,SObserved,SReconstructed)');
disp(' ');
disp('DESCRIPTION');
disp('   The function performs Regularized Wiener Inverse (RWI) reconstruction ');
disp('   of a 3D object Y from its blurred and noisy observations Z obtained at ');
disp('   positions specified at SObserved. The reconstructed is done for slices');
disp('   which are specified to be at locations SReconstructed. Wiener filter');
disp('   uses the estimate EY of a true signal.'); 
disp(' ');
disp('   Z - is an set of blurred and noisy observations of a true object.');
disp(' ');
disp('   V - is a matrix of PSFs (Vij) which correspond to the observation of the ');
disp('   object slice j from focusing at the position i.');
disp(' ');
disp('   EY - is an estimate of a true signal.');
disp(' ');
disp('   RIreg2 - is a regularization parameter r^2.');
disp(' ');
disp('   SObserved - is a set of values i, which specify the positions from');
disp('   which observations Z are done.');
disp(' ');
disp('   SReconstructed - is a set of values which specify positions of slices ');
disp('   which should be reconstructed. Observed and reconstructed slices are');
disp('   not necessary should be the same.');
disp(' ');
disp('RETURNS');
disp('   YRWI - is a reconstructed object');
disp(' ');
disp('   YRWSTDH - is a standard deviation for each slice estimate.');
disp(' ');
disp('For more details read 9.6 "3D Inverse" part of the book.');
disp(' ');
disp('Dmitriy Paliy, Tampere University of Technology');
disp('Updated 31-01-2008');
disp('dmitriy.paliy@tut.fi');
return;
end;

[yN, xN] = size(Z{1});

Vtransp = V';
for sobs=1:length(SObserved), ZZ(sobs,:,:) = fft2(Z{sobs}); end;

%%
%%%%%%%%%% RW %%%%%%%%%%%%%
% Y*Y'
% Ytransp = Y';
% tic
for sorig1=1:length(SReconstructed),
    V1 = zeros(xN,xN);
    V1(:,:) = EY{sorig1};
    
    for sorig2=1:length(SReconstructed),
        V2 = zeros(xN,xN);
        V2(:,:) = EY{sorig2};
            
        a = zeros(1,1,xN,xN);
        a = fft2(V1).*conj(fft2(V2));
                
        YY{sorig1,sorig2} = a;
    end;
end;


%%
% disp('============== YYV RW ==============')
% (Y*Y')*V'
for sorig1=1:length(SReconstructed),
    for sobs=1:length(SObserved),
        YYV(sorig1,sobs,1:xN,1:xN) = zeros(xN,xN);
    
        for sorig2=1:length(SReconstructed),
            V1 = zeros(xN,xN);
            V1(:,:) = YY{sorig1,sorig2};

            V2 = zeros(xN,xN);
            V2(:,:) = Vtransp{sorig2,sobs};
            
            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = V1.*conj(V2);
            % a(1,1,:,:) = V1.*conj(fft2(V2));
        
            YYV(sorig1,sobs,1:xN,1:xN) = YYV(sorig1,sobs,1:xN,1:xN) + a;
        end;
    end;
end;

%%
% disp('============== VYYV RW ==============')
% V*(Y*Y')*V'
clear VYYVtmp;
for sobs1=1:length(SObserved),
    for sobs2=1:length(SObserved),
        VYYV(sobs1,sobs2,1:xN,1:xN) = zeros(xN,xN);
    
        for sorig=1:length(SReconstructed),
            V1 = zeros(xN,xN);
            V1(:,:) = V{sobs1,sorig};

            V2 = zeros(xN,xN);
            V2(:,:) = squeeze(YYV(sorig,sobs2,1:xN,1:xN));
            
            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = V1.*V2;
        
            VYYV(sobs1,sobs2,1:xN,1:xN) = VYYV(sobs1,sobs2,1:xN,1:xN) + a;
        end;
        
        VYYVtmp(sobs1,sobs2,1:xN*xN) = reshape(VYYV(sobs1,sobs2,:,:),1,xN*xN);
    end;
    VYYVtmp(sobs1,sobs1,:) = VYYVtmp(sobs1,sobs1,:) + (xN*yN)*RWreg2*sigmanoise_arr(sobs)^2;
end;

% xN*yN*


%%
% (V*(Y*Y')*V' +eps*I)^-1
for i=1:xN*xN,
        a = VYYVtmp(:,:,i);
        
        % detRW(i) = cond(real(a));
        
        a = inv(a);
               
        VYYVtmp(:,:,i) = a;
end;

%%
% reshape (V*(Y*Y')*V' +eps*I)^-1
for sobs1=1:length(SObserved),
    for sobs2=1:length(SObserved),
        VYYV(sobs1,sobs2,:,:) = reshape(VYYVtmp(sobs1,sobs2,1:xN*xN),xN,xN);
    end;
end;


%%
for sorig1=1:length(SReconstructed),
    for sobs1=1:length(SObserved),
        YYV_VYYV(sorig1,sobs1,1:xN,1:xN) = zeros(xN,xN);
        % VVVtransp(sobs1,sorig1,1:xN,1:xN) = zeros(xN,xN);
    
        for sobs2=1:length(SObserved),
            % compute (V'*V +eps*I)^-1*V 
            V1 = zeros(xN,xN);
            V1(:,:) = squeeze(YYV(sorig1,sobs2,:,:));

            V2 = zeros(xN,xN);
            V2(:,:) = squeeze(VYYV(sobs2,sobs1,:,:));

            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = V1.*V2;
            
            YYV_VYYV(sorig1,sobs1,:,:) = YYV_VYYV(sorig1,sobs1,:,:) + a;
        end;
    end;
end;


%%
% disp('============== Reconstruction RW ==============')
% warning off all
% compute V'*Z/(V'*V +eps*I)
YRWestFFT = zeros(length(SReconstructed),xN,xN);

for i=1:xN,
    for j=1:xN,
        % a = YYV(:,:,i,j);
        
        % b = VYYV(:,:,i,j)  + sigmanoise*eye(length(SObserved)) + reg2*eye(length(SObserved));
        % b = VYYV(:,:,i,j);
        
        % HW = a*b;        
        
        c = squeeze(ZZ(:,i,j));
        

        HW = YYV_VYYV(:,:,i,j);
        
        YRWestFFT(:,i,j) = HW*c(:);
    end;
end;


%%
% STDH PsiVVPsi
for sorig1=1:length(SReconstructed),
    for sorig2=1:length(SReconstructed),
        PsiVVPsiRW(sorig1,sorig2,1:xN,1:xN) = zeros(xN,xN);

        for sobs1=1:length(SObserved),
            V1 = squeeze(YYV_VYYV(sorig1,sobs1,:,:));

            V2 = squeeze(YYV_VYYV(sorig2,sobs1,:,:));

            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = V1.*conj(V2);

            PsiVVPsiRW(sorig1,sorig2,:,:) = PsiVVPsiRW(sorig1,sorig2,:,:) + a;
        end;
    end;
end;

YRWI = YRWestFFT;
YRWSTDH = PsiVVPsiRW;
