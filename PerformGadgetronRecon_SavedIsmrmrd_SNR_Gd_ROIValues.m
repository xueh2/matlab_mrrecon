
function res_table = PerformGadgetronRecon_SavedIsmrmrd_SNR_Gd_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, report_column, processing_always)
% res_table = PerformGadgetronRecon_SavedIsmrmrd_Matlab_ROIValues(PerfTable, resDir, contourDir, stress_column, rest_column, ischemia_column, hct_column, fixed_HCT, report_column, processing_always)

if(nargin<11)
    processing_always = 1;
end

scanInd = [];
patientID = [];
scanDate = [];
scanTime = [];
age = [];
gender = [];
stressHB = [];
restHB = [];
hematocrit = [];

sSNR = [];
rSNR = [];

sGd = [];
rGd = [];

num_column = size(PerfTable, 2);
num = size(PerfTable, 1)-1;
for n=1:num
       
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    if(ischemia_column<=num_column)
        Category = PerfTable{n+1, ischemia_column};
    else
        Category = 1;
    end
    
    disp([num2str(n) ' out of ' num2str(num) ' - ' PerfTable{n+1, report_column} ' - Processing : '  ' - ' PerfTable{n+1, rest_column}]); 
    disp(['==================================================================']);  
       
    [configName, scannerID, pID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);

    figDir = fullfile(resDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_Figure']) 
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' pID '_' studyID '_' study_dates '_ROI'])
    
    stressDir = fullfile(resDir, study_dates, stressCase)
    restDir = fullfile(resDir, study_dates, restCase)
    
    if(~isempty(fixed_HCT))
        if(fixed_HCT>=0 & fixed_HCT<=1)
            HCT = fixed_HCT;
        end
    else
        HCT = PerfTable{n+1, hct_column};
        if(isnan(HCT))
            HCT = 0;
        else
            if(HCT>1)
                HCT = HCT/100;
            end
        end
    end    
    
    disp([num2str(n) ' out of ' num2str(num) ' - HCT : ' num2str(HCT)]); 
   
    hct_str = num2str(HCT);
    ind = find(hct_str=='.');
    if(~isempty(ind))
        hct_str(ind(:)) = 'p';
    end
    suffix = [study_dates '_hct' hct_str];
    
    if(isFileExist(fullfile(roiDir, 'myo_stress1.mat')))
        s1_roi = 'myo_stress1.mat';
        s2_roi = 'myo_stress2.mat';
        s3_roi = 'myo_stress3.mat';
        
        r1_roi = 'myo_rest1.mat';
        r2_roi = 'myo_rest2.mat';
        r3_roi = 'myo_rest3.mat';
    else
        s1_roi = 's1.mat';
        s2_roi = 's2.mat';
        s3_roi = 's3.mat';
        
        r1_roi = 'r1.mat';
        r2_roi = 'r2.mat';
        r3_roi = 'r3.mat';
    end
        
    if(isFileExist(fullfile(roiDir, s1_roi)) & isFileExist(fullfile(roiDir, r1_roi)))
        
        scanInd = [scanInd; n];
        patientID = [patientID; {pID}];
        scanDate = [scanDate; study_dates];
        scanTime = [scanTime; study_time];
        age = [age; PerfTable{n+1, 14}];
        gender = [gender; PerfTable{n+1, 15}];
        stressHB = [stressHB; PerfTable{n+1, 9}];
        restHB = [restHB; PerfTable{n+1, 13}];
        hematocrit = [hematocrit; HCT];
        
        disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]);       
        
        two_ROI = 0;
       
        SNRResult_file = fullfile(figDir, ['SNR_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
        Gd_Result_file = fullfile(figDir, ['Gd_PerfResult_' scannerID '_' pID '_' studyID '_' study_dates '.mat']);
        
        has_result_file = 0;
        if(~processing_always & isFileExist(SNRResult_file) & isFileExist(Gd_Result_file))

            snr_v = load(SNRResult_file);
                
            s_SNR = snr_v.s_SNR;
            r_SNR = snr_v.r_SNR;

            sSNR = [sSNR; max(s_SNR(4:end, :))];
            rSNR = [rSNR; max(r_SNR(4:end, :))];
            
            gd_v = load(Gd_Result_file);
                
            s_Gd = gd_v.s_Gd;
            r_Gd = gd_v.r_Gd;

            sGd = [sGd; max(s_Gd(1:end, :))];
            rGd = [rGd; max(r_Gd(1:end, :))];
        else        
            cd(roiDir)
            v = PerfTable(n+1, :);

            nV = numel(v);

            if(isFileExist(fullfile(roiDir, s1_roi)))
                s1 = load(fullfile(roiDir, s1_roi));
            else
                s1 = [];
            end
            
            if(isFileExist(fullfile(roiDir, s2_roi)))
                s2 = load(fullfile(roiDir, s2_roi));
            else
                s2 = [];
            end
            
            
            if(isFileExist(fullfile(roiDir, s3_roi)))
                s3 = load(fullfile(roiDir, s3_roi));
            else
                s3 = [];
            end

            if(isFileExist(fullfile(roiDir, r1_roi)))
                r1 = load(fullfile(roiDir, r1_roi));
            else
                r1 = [];
            end 
            
            if(isFileExist(fullfile(roiDir, r2_roi)))
                r2 = load(fullfile(roiDir, r2_roi));
            else
                r2 = [];
            end
            
            if(isFileExist(fullfile(roiDir, r3_roi)))
                r3 = load(fullfile(roiDir, r3_roi));
            else
                r3 = [];
            end

            % ---------------------------------------------
            % SNR
            
            if(~isFileExist(SNRResult_file))
                sdata1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'moco_0_MAG.hdr'));
                sdata2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'moco_1_MAG.hdr'));
                sdata3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'moco_2_MAG.hdr'));

                cd(stressDir)
                s_gfactor = readGTPlusExportImageSeries_Squeeze(300);
                s_gfactor = flipdim(s_gfactor, 2);

                rdata1 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_0_MAG.hdr'));
                rdata2 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_1_MAG.hdr'));
                rdata3 = analyze75read(fullfile(restDir, 'DebugOutput', 'moco_2_MAG.hdr'));

                cd(restDir)
                r_gfactor = readGTPlusExportImageSeries_Squeeze(300);
                r_gfactor = flipdim(r_gfactor, 2);

                sdata1 = flipdim(sdata1, 2);
                sdata2 = flipdim(sdata2, 2);
                sdata3 = flipdim(sdata3, 2);
                rdata1 = flipdim(rdata1, 2);
                rdata2 = flipdim(rdata2, 2);
                rdata3 = flipdim(rdata3, 2);

                snr_s1 = 25 * sdata1 ./ squeeze(s_gfactor(:,:,1,:));
                snr_s2 = 25 * sdata2 ./ squeeze(s_gfactor(:,:,2,:));
                snr_s3 = 25 * sdata3 ./ squeeze(s_gfactor(:,:,3,:));

                snr_r1 = 25 * rdata1 ./ squeeze(r_gfactor(:,:,1,:));
                snr_r2 = 25 * rdata2 ./ squeeze(r_gfactor(:,:,2,:));
                snr_r3 = 25 * rdata3 ./ squeeze(r_gfactor(:,:,3,:));

                cd(roiDir)
                if(~isempty(s1))
                    BW1=zeros(size(sdata1(:,:,1)));
                    BW1=roipoly(sdata1(:,:,1), s1(1).ROI_info_table(1).ROI_x_coordinates, s1(1).ROI_info_table(1).ROI_y_coordinates);
                    index1=find(BW1 >0);
                end

                if(~isempty(s2))
                    BW2=zeros(size(sdata2(:,:,1)));
                    BW2=roipoly(sdata2(:,:,1), s2(1).ROI_info_table(1).ROI_x_coordinates, s2(1).ROI_info_table(1).ROI_y_coordinates);
                    index2=find(BW2 >0);
                end

                if(~isempty(s3))
                    BW3=zeros(size(sdata3(:,:,1)));
                    BW3=roipoly(sdata3(:,:,1), s3(1).ROI_info_table(1).ROI_x_coordinates, s3(1).ROI_info_table(1).ROI_y_coordinates);
                    index3=find(BW3 >0);
                end

                if(~isempty(r1))
                    rBW1=zeros(size(rdata1(:,:,1)));
                    rBW1=roipoly(rdata1(:,:,1), r1(1).ROI_info_table(1).ROI_x_coordinates, r1(1).ROI_info_table(1).ROI_y_coordinates);
                    rindex1=find(rBW1 >0);
                end

                if(~isempty(r2))
                    rBW2=zeros(size(rdata2(:,:,1)));
                    rBW2=roipoly(rdata2(:,:,1), r2(1).ROI_info_table(1).ROI_x_coordinates, r2(1).ROI_info_table(1).ROI_y_coordinates);
                    rindex2=find(rBW2 >0);
                end

                if(~isempty(r3))
                    rBW3=zeros(size(rdata3(:,:,1)));
                    rBW3=roipoly(rdata3(:,:,1), r3(1).ROI_info_table(1).ROI_x_coordinates, r3(1).ROI_info_table(1).ROI_y_coordinates);
                    rindex3=find(rBW3 >0);
                end
            
                % gfactor scaled by 100
                % data is scaled by 4
                nRep = size(sdata1, 3);

                s_SNR = zeros(nRep, 3);
                r_SNR = zeros(nRep, 3);

                for rr=1:nRep
                    sd1 = snr_s1(:,:,rr);
                    sd2 = snr_s2(:,:,rr);
                    sd3 = snr_s3(:,:,rr);

                    rd1 = snr_r1(:,:,rr);
                    rd2 = snr_r2(:,:,rr);
                    rd3 = snr_r3(:,:,rr);

                    v1 = 0;
                    if(~isempty(s1))
                        v1 = mean(sd1(index1));
                    end

                    v2 = 0;
                    if(~isempty(s2))
                        v2 = mean(sd2(index2));
                    end

                    v3 = 0;
                    if(~isempty(s3))
                        v3 = mean(sd3(index3));
                    end

                    s_SNR(rr, :) = [v1 v2 v3];

                    v1 = 0;
                    if(~isempty(r1))
                        v1 = mean(rd1(rindex1));
                    end

                    v2 = 0;
                    if(~isempty(r2))
                        v2 = mean(rd2(rindex2));
                    end

                    v3 = 0;
                    if(~isempty(r3))
                        v3 = mean(rd3(rindex3));
                    end
                    r_SNR(rr, :) = [v1 v2 v3];                
                end

                sSNR = [sSNR; max(s_SNR(4:end, :))];
                rSNR = [rSNR; max(r_SNR(4:end, :))];

                s_gfactor_map = s_gfactor(:,:,:,1) /100;
                r_gfactor_map = r_gfactor(:,:,:,1) /100;
                save(SNRResult_file, 's_SNR', 'r_SNR', 'snr_s1', 'snr_s2', 'snr_s3', 'snr_r1', 'snr_r2', 'snr_r3', 's_gfactor_map', 'r_gfactor_map');
                
                copyfile(SNRResult_file, roiDir);
            else
                snr_v = load(SNRResult_file);
                
                s_SNR = snr_v.s_SNR;
                r_SNR = snr_v.r_SNR;
                
                sSNR = [sSNR; max(s_SNR(4:end, :))];
                rSNR = [rSNR; max(r_SNR(4:end, :))];
            end

            % ---------------------------------------------
            % Gd concentration
            if(~isFileExist(Gd_Result_file))
                sdata1 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_PSIR_0.hdr'));
                sdata2 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_PSIR_1.hdr'));
                sdata3 = analyze75read(fullfile(stressDir, 'DebugOutput', 'CASignal_Perf_PSIR_2.hdr'));

                rdata1 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_PSIR_0.hdr'));
                rdata2 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_PSIR_1.hdr'));
                rdata3 = analyze75read(fullfile(restDir, 'DebugOutput', 'CASignal_Perf_PSIR_2.hdr'));

                sdata1 = flipdim(sdata1, 2);
                sdata2 = flipdim(sdata2, 2);
                sdata3 = flipdim(sdata3, 2);
                rdata1 = flipdim(rdata1, 2);
                rdata2 = flipdim(rdata2, 2);
                rdata3 = flipdim(rdata3, 2);

                cd(roiDir)
                if(~isempty(s1))
                    BW1=zeros(size(sdata1(:,:,1)));
                    BW1=roipoly(sdata1(:,:,1), s1(1).ROI_info_table(1).ROI_x_coordinates, s1(1).ROI_info_table(1).ROI_y_coordinates);
                    index1=find(BW1 >0);
                end

                if(~isempty(s2))
                    BW2=zeros(size(sdata2(:,:,1)));
                    BW2=roipoly(sdata2(:,:,1), s2(1).ROI_info_table(1).ROI_x_coordinates, s2(1).ROI_info_table(1).ROI_y_coordinates);
                    index2=find(BW2 >0);
                end

                if(~isempty(s3))
                    BW3=zeros(size(sdata3(:,:,1)));
                    BW3=roipoly(sdata3(:,:,1), s3(1).ROI_info_table(1).ROI_x_coordinates, s3(1).ROI_info_table(1).ROI_y_coordinates);
                    index3=find(BW3 >0);
                end

                if(~isempty(r1))
                    rBW1=zeros(size(rdata1(:,:,1)));
                    rBW1=roipoly(rdata1(:,:,1), r1(1).ROI_info_table(1).ROI_x_coordinates, r1(1).ROI_info_table(1).ROI_y_coordinates);
                    rindex1=find(rBW1 >0);
                end

                if(~isempty(r2))
                    rBW2=zeros(size(rdata2(:,:,1)));
                    rBW2=roipoly(rdata2(:,:,1), r2(1).ROI_info_table(1).ROI_x_coordinates, r2(1).ROI_info_table(1).ROI_y_coordinates);
                    rindex2=find(rBW2 >0);
                end

                if(~isempty(r3))
                    rBW3=zeros(size(rdata3(:,:,1)));
                    rBW3=roipoly(rdata3(:,:,1), r3(1).ROI_info_table(1).ROI_x_coordinates, r3(1).ROI_info_table(1).ROI_y_coordinates);
                    rindex3=find(rBW3 >0);
                end
                
                nRep = size(sdata1, 3);

                s_Gd = zeros(nRep, 3);
                r_Gd = zeros(nRep, 3);

                for rr=1:nRep
                    sd1 = sdata1(:,:,rr);
                    sd2 = sdata2(:,:,rr);
                    sd3 = sdata3(:,:,rr);

                    rd1 = rdata1(:,:,rr);
                    rd2 = rdata2(:,:,rr);
                    rd3 = rdata3(:,:,rr);

                    v1 = 0;
                    if(~isempty(s1))
                        v1 = mean(sd1(index1));
                    end

                    v2 = 0;
                    if(~isempty(s2))
                        v2 = mean(sd2(index2));
                    end

                    v3 = 0;
                    if(~isempty(s3))
                        v3 = mean(sd3(index3));
                    end

                    s_Gd(rr, :) = [v1 v2 v3];

                    v1 = 0;
                    if(~isempty(r1))
                        v1 = mean(rd1(rindex1));
                    end

                    v2 = 0;
                    if(~isempty(r2))
                        v2 = mean(rd2(rindex2));
                    end

                    v3 = 0;
                    if(~isempty(r3))
                        v3 = mean(rd3(rindex3));
                    end
                    r_Gd(rr, :) = [v1 v2 v3];                
                end

                sGd = [sGd; max(s_Gd(1:end, :))];
                rGd = [rGd; max(r_Gd(1:end, :))];

                save(Gd_Result_file, 's_Gd', 'r_Gd');
                
                copyfile(Gd_Result_file, roiDir);
            else
                gd_v = load(Gd_Result_file);
                
                s_Gd = gd_v.s_Gd;
                r_Gd = gd_v.r_Gd;
                
                sGd = [sGd; max(s_Gd(1:end, :))];
                rGd = [rGd; max(r_Gd(1:end, :))];
            end                                      
        end
    else
        disp(['Cannot find the ROI ... ']);
        winopen(figDir);
        mkdir(roiDir);
        winopen(roiDir);
        pause;
    end
end

res_table = table(scanInd, patientID, scanDate, scanTime, age, gender, stressHB, restHB, hematocrit, ... 
                sSNR, rSNR, sGd, rGd);

disp('=======================================================================');
disp(['Stress SNR - ' num2str(mean(sSNR(:))) '+/-' num2str(std(sSNR(:)))]);
disp(['Rest SNR - ' num2str(mean(rSNR(:))) '+/-' num2str(std(rSNR(:)))]);

end

function res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, a1, a2, a3, HCT)
% res = PerformGadgetronRecon_Matlab_ROIValues_OneCase(s1, s2, s3, a1, a2, a3)
% given the ROIs s1/s2/s3, get the flow and other values
% res has [flow, Ki, E, Visf, Vp, PS, SD, Ki_MF, Ki_Fermi, Ki_TwoCompExp, Ki_BTEX]

if(isfield(a1, 'flowmaps_grappa_PSIR'))
    
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR, a2.flowmaps_grappa_PSIR, a3.flowmaps_grappa_PSIR, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    EMap1 = a1.Ki_whole_grappa_PSIR(:,:,4) ./ (a1.flowmaps_grappa_PSIR(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR(:,:,4) ./ (a2.flowmaps_grappa_PSIR(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR(:,:,4) ./ (a3.flowmaps_grappa_PSIR(:,:,end)+eps);
    
    [f1, f2, f3] = get_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR, a2.grappa_interVolumeMap_grappa_PSIR, a3.grappa_interVolumeMap_grappa_PSIR, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR, a2.blood_volume_maps_grappa_PSIR, a3.blood_volume_maps_grappa_PSIR, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);

    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR, a2.PS_maps_grappa_PSIR, a3.PS_maps_grappa_PSIR, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR, a2.SD_maps_grappa_PSIR, a3.SD_maps_grappa_PSIR, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR, a2.flowmaps_grappa_PSIR, a3.flowmaps_grappa_PSIR, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR, a2.grappa_interVolumeMap_grappa_PSIR, a3.grappa_interVolumeMap_grappa_PSIR, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR, a2.blood_volume_maps_grappa_PSIR, a3.blood_volume_maps_grappa_PSIR, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);

        [f1, f2, f3] = get_2nd_roi_values(a1.PS_maps_grappa_PSIR, a2.PS_maps_grappa_PSIR, a3.PS_maps_grappa_PSIR, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.SD_maps_grappa_PSIR, a2.SD_maps_grappa_PSIR, a3.SD_maps_grappa_PSIR, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,1)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,2)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,3)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR(:,:,4)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
elseif(isfield(a1, 'flowmaps_grappa_PSIR_OnlyGlobalSearch'))
    
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR_OnlyGlobalSearch, a2.flowmaps_grappa_PSIR_OnlyGlobalSearch, a3.flowmaps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    EMap1 = a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a1.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a2.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4) ./ (a3.flowmaps_grappa_PSIR_OnlyGlobalSearch(:,:,end)+eps);
    
    [f1, f2, f3] = get_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a2.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a3.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a2.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a3.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);

    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR_OnlyGlobalSearch, a2.PS_maps_grappa_PSIR_OnlyGlobalSearch, a3.PS_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR_OnlyGlobalSearch, a2.SD_maps_grappa_PSIR_OnlyGlobalSearch, a3.SD_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR_OnlyGlobalSearch, a2.flowmaps_grappa_PSIR_OnlyGlobalSearch, a3.flowmaps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a2.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, a3.grappa_interVolumeMap_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a2.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, a3.blood_volume_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);

        [f1, f2, f3] = get_2nd_roi_values(a1.PS_maps_grappa_PSIR_OnlyGlobalSearch, a2.PS_maps_grappa_PSIR_OnlyGlobalSearch, a3.PS_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.SD_maps_grappa_PSIR_OnlyGlobalSearch, a2.SD_maps_grappa_PSIR_OnlyGlobalSearch, a3.SD_maps_grappa_PSIR_OnlyGlobalSearch, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,1)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,2)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,3)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_OnlyGlobalSearch(:,:,4)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
elseif(isfield(a1, 'flowmaps_grappa_PSIR_without_R2Star'))
    
    [f1, f2, f3] = get_roi_values(a1.flowmaps_grappa_PSIR_without_R2Star, a2.flowmaps_grappa_PSIR_without_R2Star, a3.flowmaps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.flow = [f1.m f2.m f3.m];

    EMap1 = a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a1.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    EMap2 = a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a2.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    EMap3 = a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4) ./ (a3.flowmaps_grappa_PSIR_without_R2Star(:,:,end)+eps);
    
    [f1, f2, f3] = get_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
    res.E = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.Visf = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.blood_volume_maps_grappa_PSIR_without_R2Star, a2.blood_volume_maps_grappa_PSIR_without_R2Star, a3.blood_volume_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.Vp = [f1.m f2.m f3.m] * (1-HCT);

    [f1, f2, f3] = get_roi_values(a1.PS_maps_grappa_PSIR_without_R2Star, a2.PS_maps_grappa_PSIR_without_R2Star, a3.PS_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.PS = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values(a1.SD_maps_grappa_PSIR_without_R2Star, a2.SD_maps_grappa_PSIR_without_R2Star, a3.SD_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
    res.SD = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), s1, s2, s3);
    res.Ki_MF = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), s1, s2, s3);
    res.Ki_Fermi = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), s1, s2, s3);
    res.Ki_TwoCompExp = [f1.m f2.m f3.m];

    [f1, f2, f3] = get_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), s1, s2, s3);
    res.Ki_BTEX = [f1.m f2.m f3.m];

    %% if having second roi
    if(numel(s1.ROI_info_table)==2)

        [f1, f2, f3] = get_2nd_roi_values(a1.flowmaps_grappa_PSIR_without_R2Star, a2.flowmaps_grappa_PSIR_without_R2Star, a3.flowmaps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.flow_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(EMap1, EMap2, EMap3, s1, s2, s3);
        res.E_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a2.grappa_interVolumeMap_grappa_PSIR_without_R2Star, a3.grappa_interVolumeMap_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.Visf_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.blood_volume_maps_grappa_PSIR_without_R2Star, a2.blood_volume_maps_grappa_PSIR_without_R2Star, a3.blood_volume_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.Vp_i = [f1.m f2.m f3.m] * (1-HCT);

        [f1, f2, f3] = get_2nd_roi_values(a1.PS_maps_grappa_PSIR_without_R2Star, a2.PS_maps_grappa_PSIR_without_R2Star, a3.PS_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.PS_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values(a1.SD_maps_grappa_PSIR_without_R2Star, a2.SD_maps_grappa_PSIR_without_R2Star, a3.SD_maps_grappa_PSIR_without_R2Star, s1, s2, s3);
        res.SD_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,1)), s1, s2, s3);
        res.Ki_MF_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,2)), s1, s2, s3);
        res.Ki_Fermi_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,3)), s1, s2, s3);
        res.Ki_TwoCompExp_i = [f1.m f2.m f3.m];

        [f1, f2, f3] = get_2nd_roi_values( squeeze(a1.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a2.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), squeeze(a3.Ki_whole_grappa_PSIR_without_R2Star(:,:,4)), s1, s2, s3);
        res.Ki_BTEX_i = [f1.m f2.m f3.m];
    end
end

end

function [f1, f2, f3] = get_roi_values(a1, a2, a3, s1, s2, s3)
    if(isempty(s1))
        f1.m = -1;
    else
        f1 = roi_statistics(flipdim(a1(:,:,end), 2), s1.ROI_info_table(1,1));
    end
    
    if(isempty(s2))
        f2.m = -1;
    else
        f2 = roi_statistics(flipdim(a2(:,:,end), 2), s2.ROI_info_table(1,1));
    end
    
    if(isempty(s3))
        f3.m = -1;
    else
        f3 = roi_statistics(flipdim(a3(:,:,end), 2), s3.ROI_info_table(1,1));
    end
end

function [f1, f2, f3] = get_2nd_roi_values(a1, a2, a3, s1, s2, s3)

    if(~isempty(s1) & numel(s1.ROI_info_table)==2)
        f1 = roi_statistics(flipdim(a1(:,:,end), 2), s1.ROI_info_table(2,1));
    else
        f1.m = -1;
    end
    
    if(~isempty(s2) & numel(s2.ROI_info_table)==2)
        f2 = roi_statistics(flipdim(a2(:,:,end), 2), s2.ROI_info_table(2,1));
    else
        f2.m = -1;
    end
    
    if(~isempty(s3) & numel(s3.ROI_info_table)==2)
        f3 = roi_statistics(flipdim(a3(:,:,end), 2), s3.ROI_info_table(2,1));
    else
        f3.m = -1;
    end
end

