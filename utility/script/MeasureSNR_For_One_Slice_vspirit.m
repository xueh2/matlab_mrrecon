
currSNR = cell(1);
currSNR = [{subdir{d}}];
currSNR = [currSNR {sDir}];

if ( isempty(strfind(currDir, 'STRESS'))==1 )
    isStress = 0;
else
    isStress = 1;
end
    
if ( isempty(strfind(currDir, '_R3_'))~=1 )        
    R = 3;
end

if ( isempty(strfind(currDir, '_3_'))~=1 )        
    R = 3;
end

if ( isempty(strfind(currDir, '_R4_'))~=1 )        
    R = 4;
end

if ( isempty(strfind(currDir, '_r4_'))~=1 )        
    R = 4;
end

if ( isempty(strfind(currDir, '_4_'))~=1 )        
    R = 4;
end

if ( isempty(strfind(currDir, '_R5_'))~=1 )
    R = 5;
end

if ( isempty(strfind(currDir, '_r5_'))~=1 )
    R = 5;
end

if ( isempty(strfind(currDir, '_5_'))~=1 )        
    R = 5;
end
R

currSNR = [currSNR R isStress];

BW = [];

if ( isFileExist('HeartMask_vspirit.mat') )
    load('HeartMask_vspirit.mat')
end

%% linear
[data, header] = Matlab_LoadAnalyze('Linear_Recon.hdr');
if ( isempty(BW) )
    [Cutoff, Variance, SNR_Linear, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, [], noiseFromBW);
else
    [Cutoff, Variance, SNR_Linear, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
end
SNR_Linear

%% vlinear        
[data, header] = Matlab_LoadAnalyze('Linear_Recon_DstThres0.01.hdr');
[Cutoff, Variance, SNR_vLinear, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
SNR_vLinear

%% TSPIRIT        
[data, header] = Matlab_LoadAnalyze('SPIRiT_Recon_noAdaptive_ChaThres0.01_Wav0.0025_TV0.0001_DW-1.hdr');
[Cutoff, Variance, SNR_tSpirit, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
SNR_tSpirit

%% SPIRIT - average all
[data, header] = Matlab_LoadAnalyze('.hdr');
[Cutoff, Variance, SNR_SpiritAveragedAll, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
SNR_SpiritAveragedAll

%% TSPIRIT virtual
[data, header] = Matlab_LoadAnalyze('VSPIRiT_Recon_noAdaptive_ChaThres0.01_Wav0.0025_TV0.0001_DW-1.hdr');
[Cutoff, Variance, SNR_vSpirit, gfactor, BW] = Measure_SNR_3D(data, defineROIFlag, noiseMPorEigenValue, BW, noiseFromBW);
SNR_vSpirit
    
currSNR = [currSNR {[SNR_Linear SNR_vLinear SNR_tSpirit SNR_SpiritAveragedAll SNR_vSpirit]} ];

save(SNRMat, 'currSNR');
