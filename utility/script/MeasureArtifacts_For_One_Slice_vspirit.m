
currAS = cell(1);
currAS = [{subdir{d}}];
currAS = [currAS {sDir}];

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

if ( isempty(strfind(currDir, '_4_'))~=1 )        
    R = 4;
end

if ( isempty(strfind(currDir, '_R5_'))~=1 )
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
score_Linear = Compute_Cine_Recon_Aliasing_Score(data, R);

%% vlinear        
[data, header] = Matlab_LoadAnalyze('Linear_Recon_DstThres0.01.hdr');
score_vLinear = Compute_Cine_Recon_Aliasing_Score(data, R);

%% TSPIRIT        
[data, header] = Matlab_LoadAnalyze('SPIRiT_Recon_noAdaptive_ChaThres0.01_Wav0.0025_TV0.0001_DW-1.hdr');
score_tSpirit = Compute_Cine_Recon_Aliasing_Score(data, R);

%% SPIRIT - average all
[data, header] = Matlab_LoadAnalyze('.hdr');
score_SpiritAveragedAll = Compute_Cine_Recon_Aliasing_Score(data, R);

%% TSPIRIT virtual
[data, header] = Matlab_LoadAnalyze('VSPIRiT_Recon_noAdaptive_ChaThres0.01_Wav0.0025_TV0.0001_DW-1.hdr');
score_vSpirit = Compute_Cine_Recon_Aliasing_Score(data, R);
    
currAS = [currAS {[score_Linear score_vLinear score_tSpirit score_SpiritAveragedAll score_vSpirit]} ];

save(ASMat, 'currAS');       
