% Demo of 3D Inverse Imaging
%
% A nonparametric approach to reconstruction of three-dimensional (3D)
% objects from 2D blurred and noisy observations. The restoration scheme
% incorporates regularized inverse (RI) and regularized Wiener inverse (RWI)
% filters in combination with an adaptive denoising technique LPA-ICI. 
%
% Dmitriy Paliy. Tampere University of Technology. 
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

warning off all
clear all

init= 0;
randn('seed',init);

init1= 2055615588;
rand('seed',init1);

% ---------------------------------------------------
% 3D model parameters
% ---------------------------------------------------
MODEL = 0; % MODEL = 0; corresponds to 128x128x3 3D object which consits of 
           %'Cameraman', 'Lena', and 'Testpat' strata consiquently. 
           % MODEL = 1; corresponds to 64x64x64 3D object which consists of 
           % 5 nonoverlapping spheres of different colors forming
           % 'V'-shape.
           % MODEL = 2; 'forum.png' with Gaussian blur

% ---------------------------------------------------
% creates 3D model
% ---------------------------------------------------
[Fxyz] = function_create3DModel(MODEL);

% ---------------------------------------------------
% paramaters for different models
% ---------------------------------------------------
switch MODEL
    case 1 % 'V'-spheres object
        % ---------------------------------------------------
        % level of noise in dB
        % ---------------------------------------------------
        BSNRdB = 30;
        
        % ---------------------------------------------------
        % do filtering LPA-ICI
        % ---------------------------------------------------
        DoFiltering = 0; % '0' states for no LPA-ICI filtering is done in the scheme
                         % '1' means the LPA-ICI filtering in the scheme after
                         % inversion.
                         
        RIreg2 = 0.014^2; % regularization parameter for RI with LPA-ICI
        RWreg2 = 0.11;   % regularization parameter for RWI with LPA-ICI
        
        
        % coordinates of original stratas, i.e. strata we observe
        SOriginal      = [14 62 110];

        % coordinates of observed noisy stratas, i.e. focal planes where the strata
        % observed
        SObserved      = [14 62 110];

        % coordinates of reconstructed stratas, i.e. at which planes the object is
        % reconstructed.
        SReconstructed = [14 62 110];
    case 2 % forum
        % ---------------------------------------------------
        % level of noise in dB
        % ---------------------------------------------------
        BSNRdB = 40;
        
        % ---------------------------------------------------
        % do filtering LPA-ICI
        % ---------------------------------------------------
        DoFiltering = 1; % '0' states for no LPA-ICI filtering is done in the scheme
                         % '1' means the LPA-ICI filtering in the scheme after
                         % inversion.
        if DoFiltering==0,
            % cameraman
            % RIreg2 = 0.030^2; % regularization parameter for RI
            % RWreg2 = 0.6;   % regularization parameter for RWI
            
            % forum
            RIreg2 = 0.023^2; % regularization parameter for RI
            RWreg2 = 0.65;   % regularization parameter for RWI
            
        else            
            RIreg2 = 0.015^2; % regularization parameter for RI with LPA-ICI
            RWreg2 = 0.25;   % regularization parameter for RWI with LPA-ICI
            
            GammaParameterRI = 1.45;  % Gamma parameter for ICI after RI inversion
            GammaParameterRW = 1.55; % Gamma parameter for ICI after RWI inversion
        end;

        % coordinates of original stratas, i.e. strata we observe
        SOriginal      = [1];

        % coordinates of observed noisy stratas, i.e. focal planes where the strata
        % observed
        SObserved      = [2.5]; % 2.5 corresponds to Gauss sigma=1.5
                                % 1.5 corresponds to Gauss sigma=0.5

        % coordinates of reconstructed stratas, i.e. at which planes the object is
        % reconstructed.
        SReconstructed = [1];
        
    otherwise % 'Cameraman', 'Lena' images
        % ---------------------------------------------------
        % level of noise in dB
        % ---------------------------------------------------
        BSNRdB = 40;
        
        % ---------------------------------------------------
        % do filtering LPA-ICI
        % ---------------------------------------------------
        DoFiltering = 1; % '0' states for no LPA-ICI filtering is done in the scheme
                         % '1' means the LPA-ICI filtering in the scheme after
                         % inversion.
        if DoFiltering==0,
            RIreg2 = 0.009^2; % regularization parameter for RI with LPA-ICI
            RWreg2 = 90;   % regularization parameter for RWI with LPA-ICI            
        else            
            RIreg2 = 0.006^2; % regularization parameter for RI with LPA-ICI
            RWreg2 = 0.05;   % regularization parameter for RWI with LPA-ICI

            
            GammaParameterRI = 1.2;  % Gamma parameter for ICI after RI inversion
            GammaParameterRW = 4.5; % Gamma parameter for ICI after RWI inversion
        end;

        % coordinates of original stratas, i.e. strata we observe
        SOriginal      = [1 2];

        % coordinates of observed noisy stratas, i.e. focal planes where the strata
        % observed
        SObserved      = [0 1.5 3];

        % coordinates of reconstructed stratas, i.e. at which planes the object is
        % reconstructed.
        SReconstructed = [1 2];        
end;


N = size(Fxyz,1);
xN = N; yN = N;


% generates the observations
utility_GenerateObservations3D;

V = VR;

figure('name','Observation:');
for sobs=1:length(SObserved),
    subplot(1,length(SObserved),sobs), imshow(Z{sobs},[]),
    title(['Observ. at depth j = ',num2str(SObserved(sobs))]);
end;

drawnow;

% Regularized Inverse (RI) reconstruction
utility_DeconvolutionRI3D;

% Regularized Wiener Inverse (RWI) reconstruction
utility_DeconvolutionRWI3D;
