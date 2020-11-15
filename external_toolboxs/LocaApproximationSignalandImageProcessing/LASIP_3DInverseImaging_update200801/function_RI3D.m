function [YRI, YSTDH] = function_RI3D(Z,V,RIreg2,SObserved,SReconstructed)

% The function performs Regularized Inverse (RI) reconstruction of slices 
% of a 3D object at specified positions.
%
% SYNTAX 
%   [YRI, YSTDH] = function_RI(Z,V,RIreg2,SObserved,SReconstructed)
%
% DESCRIPTION
%   The function performs Regularized Inverse (RI) reconstruction of a 3D 
%   object Y from its blurred and noisy observations Z obtained at 
%   positions specified at SObserved. The reconstructed is done for slices
%   which are specified to be at locations SReconstructed.
% 
%   Z - is an set of blurred and noisy observations of a true object.
%
%   V - is a matrix of PSFs (Vij) which correspond to the observation of the 
%   object slice j from focusing at the position i.
%
%   RIreg2 - is a regularization parameter r^2.
%
%   SObserved - is a set of values i, which specify the positions from
%   which observations Z are done.
%
%   SReconstructed - is a set of values which specify positions of slices 
%   which should be reconstructed. Observed and reconstructed slices are
%   not necessary should be the same.
%
% RETURNS
%   YRI - is a reconstructed object
%
%   YSTDH - is a standard deviation for each slice estimate.
% 
% For more details read 9.6 '3D Inverse' part of the book.
%
% Dmitriy Paliy, Tampere University of Technology, 
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

%%
if nargin<1,
disp('The function performs Regularized Inverse (RI) reconstruction of slices');
disp('of a 3D object at specified positions.');
disp(' ');
disp('SYNTAX ');
disp('   [YRI, YSTDH] = function_RI(Z,V,RIreg2,SObserved,SReconstructed)');
disp(' ');
disp('DESCRIPTION');
disp('   The function performs Regularized Inverse (RI) reconstruction of a 3D ');
disp('   object Y from its blurred and noisy observations Z obtained at ');
disp('   positions specified at SObserved. The reconstructed is done for slices');
disp('   which are specified to be at locations SReconstructed.');
disp(' ');
disp('   Z - is an set of blurred and noisy observations of a true object.');
disp(' ');
disp('   V - is a matrix of PSFs (Vij) which correspond to the observation of the ');
disp('   object slice j from focusing at the position i.');
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
disp('   YRI - is a reconstructed object');
disp(' ');
disp('   YSTDH - is a standard deviation for each slice estimate.');
disp(' ');
disp('For more details read 9.6 "3D Inverse" part of the book.');
disp(' ');
disp('Dmitriy Paliy, Tampere University of Technology.');
disp('Updated 31-01-2008');
disp('dmitriy.paliy@tut.fi');
return;
end;

[yN, xN] = size(Z{1});

Vtransp = V';

%%
% compute (V'*V)
reg1 = 0.0;
for sorig1=1:length(SReconstructed),
    for sorig2=1:length(SReconstructed),
        VVtmp(sorig1,sorig2,1:xN,1:xN) = reg1*ones(xN,xN);

        for sobs=1:length(SObserved),
            V1 = Vtransp{sorig1,sobs};

            V2 = V{sobs,sorig2};

            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = conj(V1).*V2;

            VVtmp(sorig1,sorig2,:,:) = VVtmp(sorig1,sorig2,:,:) + a;
        end;

        VV(sorig1,sorig2,1:xN*xN) = reshape(VVtmp(sorig1,sorig2,:,:),1,xN*xN);
    end;
    VV(sorig1,sorig1,1:xN*xN) = VV(sorig1,sorig1,1:xN*xN) + repmat(RIreg2,[1,1,xN*xN]);
end;



%%
% (V'*V +eps*I)^-1
for i=1:xN*xN,
    a = VV(:,:,i);

    a = inv(a);

    VV(:,:,i) = a;
end;

% VV = 1./VV;



%%
% reshape (V'*V +eps*I)^-1
for sorig1=1:length(SReconstructed),
    for sorig2=1:length(SReconstructed),
        VVtmp(sorig1,sorig2,:,:) = reshape(VV(sorig1,sorig2,1:xN*xN),xN,xN);
    end;
end;

% disp('============== VVV RI ==============')
for sorig1=1:length(SReconstructed),
    for sobs1=1:length(SObserved),
        VVV(sorig1,sobs1,1:xN,1:xN) = zeros(xN,xN);
        VVVtransp(sobs1,sorig1,1:xN,1:xN) = zeros(xN,xN);

        for sorig2=1:length(SReconstructed),
            V1(:,:) = squeeze(VVtmp(sorig1,sorig2,:,:));

            V2(:,:) = Vtransp{sorig2,sobs1};

            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = V1.*conj(V2);

            VVV(sorig1,sobs1,:,:) = VVV(sorig1,sobs1,:,:) + a;
        end;
    end;
end;


%%
% disp('============== Reconstruction RI ==============')
% warning off all
% compute V'*(V'*V +eps*I)^-1*Z
% disp('============== VVVZ RI ==============')
for sobs=1:length(SObserved), ZZ(sobs,:,:) = fft2(Z{sobs}); end;
YRIestFFT = zeros(length(SReconstructed),xN,xN);
for i=1:xN,
    for j=1:xN,
        a = VVV(:,:,i,j);

        b = ZZ(:,i,j);

        c = a*b;

        YRIestFFT(:,i,j) = c(:);
    end;
end;

%%
% STDH PsiVVPsi
for sorig1=1:length(SReconstructed),
    for sorig2=1:length(SReconstructed),
        PsiVVPsi(sorig1,sorig2,1:xN,1:xN) = zeros(xN,xN);

        for sobs1=1:length(SObserved),
            V1 = squeeze(VVV(sorig1,sobs1,:,:));

            V2 = squeeze(VVV(sorig2,sobs1,:,:));

            a = zeros(1,1,xN,xN);
            a(1,1,:,:) = V1.*conj(V2);

            PsiVVPsi(sorig1,sorig2,:,:) = PsiVVPsi(sorig1,sorig2,:,:) + a;
        end;
    end;
end;

YRI = YRIestFFT;
YSTDH = PsiVVPsi;
