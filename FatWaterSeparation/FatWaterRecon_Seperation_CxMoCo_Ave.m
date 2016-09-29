
function [Water, Fat, InPhase, OppPhase, B0Map, T2StarMap] = FatWaterRecon_Seperation_CxMoCo_Ave(measDatName, unwarppedImCombined, dataSize, ...
                                                                                                voxelsize, performMoCo, PDIRCoRegFlag, keyFrame, iters, sigma, ... 
                                                                                                interpolator, gyroRatio, B0, configFile)
% function [Water, Fat, InPhase, OppPhase, B0Map, T2StarMap] = FatWaterRecon_Seperation_CxMoCo_Ave(unwarppedImCombined, dataSize, ...
%                                                                                                voxelsize, performMoCo, PDIRCoRegFlag, keyFrame, iters, sigma, ... 
%                                                                                                interpolator, gyroRatio, B0, configFile)
%
% this function performs the fat water seperation, complex moco and averaging
%
% Inputs:
%   measDatName : meas .dat file name
%   unwarppedImCombined : reconstructed complex image after coil combination [COL LIN ECO REP SET]
%   dataSize : a vector storing the dimensions [Nfe Npe numOfCoil numOfEcho numOfRep numOfSet]
%   voxelsize : the reconstructed voxel size
%   performMoCo : if 1, the CxMoCo will be performed; 0, will not be performed
%   PDIRCoRegFlag : if 1, perform PD co-reg to IR; if 0, not
%   keyFrame : a [numOfSet 1] vector, the keyFrame for every set; first frame is frame 0
%   iters, sigma : number of iterations and sigma for non-rigid registration
%   interpolator : method of interpolation, 'BSpline' or 'Linear' or 'NN'
%   gyroRatio : for 1.5T, it is 42575575 (LarmorConstant in the protocol)
%   B0 : main field strength, for Siemens 1.5T, it is 1.494000 (NominalB0 in the protocol)
%   configFile : the configuration for fat water seperation
%
% Output:
%    unwrappedIm : reconstructed uncombined complex image [COL LIN CHA ECO REP SET]
%    fullkspace : reconstructed full kspace [COL LIN CHA ECO REP SET]
%    sensitivityMap : estimated coil sensitivity [COL LIN CHA ECO REP SET]
%    unwarppedImCombined : reconstructed complex image after coil combination [COL LIN ECO REP SET]
%    dataSize : a vector storing the dimensions [Nfe Npe numOfCoil numOfEcho numOfRep numOfSet]
%    voxelsize : the reconstructed voxel size, the zero-filling has been taken care of
%
%     ***************************************
%     *  Hui Xue  (hui-xue@siemens.com)     *
%     ***************************************

[headers,protocol]=read_dat_headers(measDatName);

Nfe = dataSize(1);
Npe = dataSize(2);
numOfCoil = dataSize(3);
numOfEcho = dataSize(4);
numOfRep = dataSize(5);
numOfSet = dataSize(6);

alTE = zeros(numel(protocol.alTE), 1);
for tt = 1:numel(alTE)
    alTE(tt) = protocol.alTE{tt};
end       
TEs = alTE(1:numOfEcho)/1e6 % 1/1000 ms to second
       
if ( performMoCo )
    Fat = zeros(Nfe, Npe, 1, numOfSet);
    Water = zeros(Nfe, Npe, 1, numOfSet);
    InPhase = zeros(Nfe, Npe, 1, numOfSet);
    OppPhase = zeros(Nfe, Npe, 1, numOfSet);
    B0Map = zeros(Nfe, Npe, 1, numOfSet);
    T2StarMap = zeros(Nfe, Npe, 1, numOfSet);
    mocoAveRes = zeros(Nfe, Npe, numOfEcho, 1, numOfSet);
else       
    Fat = zeros(Nfe, Npe, numOfRep, numOfSet);
    Water = zeros(Nfe, Npe, numOfRep, numOfSet);
    InPhase = zeros(Nfe, Npe, numOfRep, numOfSet);
    OppPhase = zeros(Nfe, Npe, numOfRep, numOfSet);
    B0Map = zeros(Nfe, Npe, numOfRep, numOfSet);
    T2StarMap = zeros(Nfe, Npe, numOfRep, numOfSet);
    mocoAveRes = zeros(Nfe, Npe, numOfEcho, numOfRep, numOfSet);
end

% moco
for set=1:numOfSet
    unwarppedImCombinedSet = unwarppedImCombined(:,:,:,:,set);    
    %% complex moco and average
    if ( numOfRep > 1 )

        strategy = 'Direct';
        inverse = 1; 
        initial = 0; 
        numOfPre = 0; 
%         iters = [32 32 32]; 
%         sigma = 12.0; 
        neighbor = 2.0; 
        stepDiv =  3.0; 
        moreIterInv = 1; 
        algo = 'GLCC'; 
        volumePreserving = 0;
        % interpolator = 'NN';

        % moco results
        mocoRes = unwarppedImCombined(:,:,:,:,set);

        % take the first echo
        firstEcho = squeeze(unwarppedImCombined(:,:,1,:,set));
        header = CreateFtkHeaderInfo(firstEcho, voxelsize);

        [moco, dx, dy, invDx, invDy] = PerformTemporalMotionCorrectionComplex(firstEcho, header, keyFrame(set), strategy, inverse, ...
                                initial, numOfPre, iters, sigma, neighbor, stepDiv, moreIterInv, algo, volumePreserving, interpolator);

        mocoRes(:,:,1,:) = moco;
        header = CreateFtkHeaderInfo(moco, voxelsize);

        % apply the deformation fields to other echo
        for echo=1:numOfEcho
            warpped = PerformComplexImageWarpping(squeeze(unwarppedImCombined(:,:,echo,:,set)), header, dx, dy, keyFrame(set), interpolator);
            mocoRes(:,:,echo,:) = warpped;
        end

        % compute averaging
        mocoAveRes(:,:,:,1,set) = mean(mocoRes, 4);
        
        % co-register the PD to IR
        if ( set == 2 & PDIRCoRegFlag )
            set1Echo1Im = mocoAveRes(:,:,1,1,1);
            set2Echo1Im = mocoAveRes(:,:,1,1,2);
            
            data(:,:,1) = set1Echo1Im;
            data(:,:,2) = set2Echo1Im;
            
            header = CreateFtkHeaderInfo(data, voxelsize);
            [moco, dx, dy, invDx, invDy] = PerformTemporalMotionCorrectionComplex(data, header, 0, strategy, inverse, ...
                                initial, numOfPre, iters, 24.0, neighbor, stepDiv, moreIterInv, algo, volumePreserving, interpolator);

            % warp other echo
            header.sizeZ = 1;
            for echo=1:numOfEcho
                warpped = PerformComplexImageWarpping(mocoAveRes(:,:,echo,1,2), header, dx(:,:,2), dy(:,:,2), 0, interpolator);
                mocoAveRes(:,:,echo,1,2) = warpped;
            end
        end
    else
        mocoAveRes(:,:,:,:,set) = unwarppedImCombinedSet;
    end
end

% fat water seperation
% for set=1:numOfSet
% 
%     unwarppedImCombinedSet = mocoAveRes(:,:,:,:,set);
%     
%     for rep=1:size(unwarppedImCombinedSet, 4)
% 
%         MultiEchoImage = unwarppedImCombinedSet(:,:,:,rep);
% 
%         % calling the fat water seperation
%         MultiEchoImage = MultiEchoImage ./ max(abs(MultiEchoImage(:)));
% 
%         header = CreateFtkHeaderInfo(MultiEchoImage, voxelsize)
%         [Fat_Rep, Water_Rep, InPhase_Rep, OppPhase_Rep, B0Map_Rep, T2StarMap_Rep] = Matlab_PerformFatWaterSeparation(single(MultiEchoImage), header, gyroRatio, B0, numOfEcho, TEs, configFile);
% 
%         Fat(:,:,rep,set) = Fat_Rep';
%         Water(:,:,rep,set) = Water_Rep';
%         InPhase(:,:,rep,set) = InPhase_Rep';
%         OppPhase(:,:,rep,set) = OppPhase_Rep';
%         B0Map(:,:,rep,set) = B0Map_Rep';
%         T2StarMap(:,:,rep,set) = T2StarMap_Rep';
%     end
% end

for rep=1:size(mocoAveRes,4)

    unwarppedImCombinedRep = mocoAveRes(:,:,:,rep,:);
    
    MultiEchoImage = squeeze(unwarppedImCombinedRep);
    MultiEchoImage = permute(MultiEchoImage, [1 2 4 3]);
    
    % calling the fat water seperation
    MultiEchoImage = MultiEchoImage ./ max(abs(MultiEchoImage(:)));

    header = CreateFtkHeaderInfo(MultiEchoImage, voxelsize)
    [Fat_Rep, Water_Rep, InPhase_Rep, OppPhase_Rep, B0Map_Rep, T2StarMap_Rep] = Matlab_PerformFatWaterSeparation(single(MultiEchoImage), header, gyroRatio, B0, numOfEcho, TEs, configFile);

    f = permute(Fat_Rep, [2 3 1]);
    w = permute(Water_Rep, [2 3 1]);
    in = permute(InPhase_Rep, [2 3 1]);
    opp = permute(OppPhase_Rep, [2 3 1]);
    
    Fat(:,:,rep,:) = f;
    Water(:,:,rep,:) = w;
    InPhase(:,:,rep,:) = in;
    OppPhase(:,:,rep,:) = opp;
    B0Map(:,:,rep,:) = repmat(B0Map_Rep', [1 1 numOfSet]);
    T2StarMap(:,:,rep,:) = repmat(T2StarMap_Rep', [1 1 numOfSet]);
end