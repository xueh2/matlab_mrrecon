
currSNR = cell(1);
currSNR = [{subdir{d}}];
currSNR = [currSNR {sDir}];

[R, isStress] = findDataInfoR3R4R5(subdir{d});
% 
% if ( (isempty(strfind(currDir, 'STRESS'))==1) & (isempty(strfind(currDir, 'stress'))==1) )
%     isStress = 0;
% else
%     isStress = 1;
% end
%     
% R = -1;
% if ( isempty(strfind(currDir, 'R3_'))~=1 )        
%     R = 3;
% end
% 
% if ( isempty(strfind(currDir, '_R3_'))~=1 )        
%     R = 3;
% end
% 
% if ( isempty(strfind(currDir, '_3_'))~=1 )        
%     R = 3;
% end
% 
% if ( isempty(strfind(currDir, '_R4_'))~=1 )        
%     R = 4;
% end
% 
% if ( isempty(strfind(currDir, '_4_'))~=1 )        
%     R = 4;
% end
% 
% if ( isempty(strfind(currDir, '_R5_'))~=1 )
%     R = 5;
% end
% 
% if ( isempty(strfind(currDir, '_r5_'))~=1 )
%     R = 5;
% end
% 
% if ( isempty(strfind(currDir, '_5_'))~=1 )        
%     R = 5;
% end
% R 

currSNR = [currSNR R isStress];

BW = [];

if ( isFileExist('HeartMask.mat') )
    load('HeartMask.mat')
end

%% ICEPAT
[IcePAT, header] = Matlab_LoadAnalyze('IcePAT_AverageAll_All_kernalIndFE2^-2_tikhonov_0p0001_preWhitening1_30_PEUsed1_Gaussian_Weak_Hanning_Strong.hdr');
score_icepat = Compute_Cine_Recon_Aliasing_Score(IcePAT, R);

if ( isempty(BW) )
    [Cutoff, Variance, SNR_icepat, gfactor, BW] = Measure_SNR_3D(IcePAT, defineROIFlag, noiseMPorEigenValue, [], noiseFromBW);
else
    [Cutoff, Variance, SNR_icepat, gfactor, BW] = Measure_SNR_3D(IcePAT, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
end
SNR_icepat

%% TSPIRIT        
[spirit, header] = Matlab_LoadAnalyze('SPIRiT_kernel9_ThresReg0p01_wavW0p0025_cgLambda0p0001_TVW0p0001_KLTSen1_Itnlim5_fft_AllGPU.hdr');        
[Cutoff, Variance, SNR_spirit, gfactor] = Measure_SNR_3D(spirit, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
SNR_spirit

score_spirit = Compute_Cine_Recon_Aliasing_Score(spirit, R);

[h, p_spirit ] = ttest2(score_icepat, score_spirit);

%% SPIRIT - average all
if ( isStress )
    filename = 'SPIRiT_AverageAll_kernel9_ThresReg0p01_wavW0p0025_cgLambda0p0001_TVW0p0001_KLTSen1_Itnlim5_fft_AllGPU.hdr';
    if ( isFileExist(filename) )
        [spirit2, header] = Matlab_LoadAnalyze('SPIRiT_AverageAll_kernel9_ThresReg0p01_wavW0p0025_cgLambda0p0001_TVW0p0001_KLTSen1_Itnlim5_fft_AllGPU.hdr');        
    else
        filename = 'SPIRiT_Recon_noAdaptive_AverageAll_ChaThres-1_Wav0.0025_TV0.0001_DW-1.hdr';
        if ( isFileExist(filename) )
            [spirit2, header] = Matlab_LoadAnalyze('SPIRiT_Recon_noAdaptive_AverageAll_ChaThres-1_Wav0.0025_TV0.0001_DW-1.hdr');
        else
            [spirit2, header] = Matlab_LoadAnalyze('SPIRiT_Recon_noAdaptive_AverageAll_ChaThres-1_Wav0.0025_TV1e-005_DW-1.hdr');
        end
    end
    [Cutoff2, Variance2, SNR_spirit2, gfactor2] = Measure_SNR_3D(spirit2, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
    SNR_spirit2

    score_spirit2 = Compute_Cine_Recon_Aliasing_Score(spirit2, R);

    [h, p_spirit2 ] = ttest2(score_icepat, score_spirit2);
    
    currSNR = [currSNR {[SNR_icepat SNR_spirit SNR_spirit2 mean(score_icepat) max(score_icepat) mean(score_spirit) max(score_spirit) p_spirit mean(score_spirit2) max(score_spirit2) p_spirit2]} {[score_icepat score_spirit score_spirit2]}];
else
    currSNR = [currSNR {[SNR_icepat SNR_spirit mean(score_icepat) max(score_icepat) mean(score_spirit) max(score_spirit) p_spirit]} {[score_icepat score_spirit]}];
end
save(SNRMat, 'currSNR');       
