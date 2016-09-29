
clear all
close all

set_UT_Dir('E')
UTDir = getenv('GTPLUS_UT_DIR')

% ------------------------------------------------------------

home = fullfile(UTDir, 'DualBolus\SSFP');

normalCases = {'20150901_A04005_4ml', ... 
            '20150911_A04027_4ml', ... 
            '20150916_A04036_4ml', ... 
            '20151001_A04062_4ml', ... 
            '20151007_A04067_4ml', ... 
            '20151013_A04076_4ml', ...
            '20151028_A04102_4ml', ...
            '20151112_A04125_4ml', ...
            '20151119_A04136_4ml', ...
            '20151124_A04141_4ml', ...
            '20151203_A04152_4ml', ...
            '20151215_A04163_4ml', ... 
            '20151218_A04173_4ml' };

infarctCases = {'20150915_A04032_4ml'};

lowEF_or_LowResponse_Cases = {'20151020_A04085_4ml', };
microVascular_or_multivessel_DiseaseCases = {};

cardiomyopathyCases = {'20151028_A04102_4ml', '20151112_A04125_4ml'};

numel(normalCases)

% -------------------------

home  = fullfile(UTDir, 'DualBolus\FLASH');

normalCases = {'20150825_A03986_4ml_s', ...
            '20150925_A04051_4ml_s', ...
            '20151023_A04093_4ml_s', ...
            '20151103_A04111_4ml_s', ...
            '20151118_A04133_4ml_s', ... 
            '20151119_A04137_4ml_s', ...
            '20151125_A04145_4ml_s', ...
            '20151217_A04169_4ml_s', ... 
            '20151222_A04177_4ml_s', ...
            '20160105_A04185_4ml_s', ...
            '20160107_A04194_4ml_s', ...
            '20160115_A04211_4ml_s' };

infarctCases = {'20150908_A04017_4ml_s', '20151006_A04065_4ml_s', '20151009_A04073_4ml_s'};

lowEF_or_LowResponse_Cases = {'20150915_A04029_4ml_s', '20151001_A04061_4ml_s', '20151216_A04165_4ml_s'};
microVascular_or_multivessel_DiseaseCases = {};

myocarditisCases = {'20150925_A04051_4ml_s'};
cardiomyopathyCases = {'20151216_A04165_4ml_s'};

numel(normalCases)

% -------------------------

home  = fullfile(UTDir, 'DualBolus\SUB_Perf\NewSeq');

normalCases = {'20150923_E07574_4ml', ...
                '20150930_E07591_4ml', ...
                '20151001_E07594_2ml', ...
                '20151013_E07616_4ml', ...
                '20151014_E07618_2ml', ...
                '20151016_E07626_4ml', ...
                '20151020_E07628_2ml', ...
                '20151028_E07643_2ml', ...
                '20151103_E07655_4ml', ...
                '20151105_E07663_2ml', ...
                '20151112_E07675_2ml', ...
                '20151117_E07684_4ml', ... 
                '20151124_E07704_2ml', ...
                '20151202_E07720_4ml'};
            
infarctCases = {'20151021_E07631_4ml', '20151027_E07640_4ml', '20151110_E07673_4ml', '20151124_E07705_4ml', '20151211_E07727_4ml', '20151215_E07733_4ml'};

lowEF_or_LowResponse_Cases = {'20150917_E07562_2ml', '20151110_E07674_2ml', '20151124_E07705_4ml'};
microVascular_or_multivessel_DiseaseCases = {'20150929_E07586_4ml', '20151110_E07674_2ml', '20151117_E07684_4ml', };

numel(normalCases)

%% start process

home_ut = home(length(UTDir)+2:end)
UTCase = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2.xml', 'grappa_res',              'grappa_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_NonLinear.xml',           'slep_res',          'slep_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Cloud = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_NonLinear_Gateway.xml',           'slep_cloud_res',          'slep_cloud_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}

UTCase_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping.xml',                                   'grappa_flow_res_new2',              'grappa_flow_ref_new2',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear.xml',                      'slep_flow_res_new',              'slep_flow_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}
UTCase_NL_Flow_Cloud = {home_ut, '20140408_13h17m13s_10066', 'VD',   'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear_Gateway.xml',        'slep_flow_cloud_res_new',              'slep_flow_cloud_ref_new',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ;}

% ------------------------------------------------------------

exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2'};
[subdir, num] = FindSubDirs(home)

debugFolder = 'E:\gtuser\mrprogs\install\DebugOutput';

isVD = 1;
isVD11 = 1;

configName = 'GT_2DT_Perfusion_NonLinear.xml';
res = 'grappa_res';

xslName = 'IsmrmrdParameterMap_Siemens_Perfusion.xsl';
deleteh5 = 0;

%% for linear recon
resDir = 'grappa_flow_res_BTEX20_TwoCompWithShifts';
resDir = 'Dotarem_r1_r2_reduced';
resDir = 'Dotarem_r1_r2_increased';

plotFlag = 1;

ResultRecordRest = [];
ResultRecordStress = [];
ResultRecordStressIschemia = [];

Ki_map_rest_normal = [];
Ki_map_stress_normal = [];

Flow_map_rest_normal = [];
Flow_map_stress_normal = [];

sickCases = [];

Ki_map_rest_sick = [];
Ki_map_stress_sick = [];

Flow_map_rest_sick = [];
Flow_map_stress_sick = [];

rest_pt = [];
stress_pt = [];

rest_pde_pt = [];
stress_pde_pt = [];

rest_t1_myo = [];
rest_t1_blood = [];

stress_t1_myo = [];
stress_t1_blood = [];

normal_cases_maps = [];

ind = 1;
indS = 1;
indS_ischemia = 1;
numNorm = numel(normalCases);

sampleinterval = 0.5;
sigmas = [1.6 4.0 5.3];
sigmaMeasure = 0.1;
thresGradient = 0.5;
   
for ii=1:num 
      
    ii
    cd(fullfile(home, subdir{ii}))
    
    if ( ~isempty(strfind( lower(subdir{1}), 't1map')) )
        continue;
    end
    
    normal = 0;
    for jj=1:numNorm
        if(strcmp(normalCases{jj}, subdir{ii})==1)
            normal = 1;
            break;
        end
    end
    
    [subdirs, nums] = FindSubDirs(fullfile(home, subdir{ii}, 'rest'))    
    if ( ~isempty(strfind( lower(subdirs{1}), 'mini')) | ~isempty(strfind(subdirs{1}, 'MINI')) )
        rest = subdirs{2};
    else
        rest = subdirs{1};
    end
    
    [subdirs, nums] = FindSubDirs(fullfile(home, subdir{ii}, 'stress'))
    if ( ~isempty(strfind( lower(subdirs{1}), 'mini')) | ~isempty(strfind(subdirs{1}, 'MINI')) )
        stress = subdirs{2};
    else
        stress = subdirs{1};
    end
    
    try
        pp = load(fullfile(home, subdir{ii}, 'EF.mat'));
        EF = pp.EF;
    catch
        EF=0;
    end
    
    normal = 0;
    for nn=1:numNorm
        if(strcmp(normalCases{nn}, subdir{ii})==1)
            normal = 1;
            break;
        end
    end
    
    infarct = 0;
    for nn=1:numel(infarctCases)
        if(strcmp(infarctCases{nn}, subdir{ii})==1)
            infarct = 1;
            break;
        end
    end
    
    restDat = fullfile(home, subdir{ii}, 'rest', rest, [rest '.h5']);
    restDir = fullfile(home, subdir{ii}, 'rest', rest, resDir)
    a = readGTPlusExportImageSeries_Squeeze(restDir, 120);
    fa = readGTPlusExportImageSeries_Squeeze(restDir, 115);
    fa = fa/100;
    fa_pde = readGTPlusExportImageSeries_Squeeze(restDir, 124);
    fa_pde = fa_pde/100;
    aif_rest = analyze75read(fullfile(restDir, 'DebugOutput', 'aif_cin_all.hdr'));
    aif_rest_without_R2Star = analyze75read(fullfile(restDir, 'DebugOutput', 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
    
    [slope_rest, timeToPeak_rest, peakTime_rest, areaUnderCurve_rest, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(aif_rest(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 0, 'Feature detection of AIF LV signal');

    aif_rest2 = aif_rest(round(peakTime_rest):end);
    [slope, timeToPeak, peakTime, areaUnderCurve, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(-aif_rest2(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 0, 'Feature detection of AIF LV signal');   
    valley_rest = round(peakTime_rest) + peakTime - 0.5;
    foot_rest = round(peakTime_rest - timeToPeak);
    
    v = mean(aif_rest_without_R2Star(1:foot_rest-1));
    aif_rest_without_R2Star = aif_rest_without_R2Star - v;
    
    dset = ismrmrd.Dataset(restDat, 'dataset');
    rest_header = ismrmrd.xml.deserialize(dset.readxml());
    acq = dset.readAcquisition(1);
    rest_t = double(acq.head.acquisition_time_stamp);
    
    sfa = size(fa);
    
    % ---------------------------
    
    stressDat = fullfile(home, subdir{ii}, 'stress', stress, [stress '.h5']);
    stressDir = fullfile(home, subdir{ii}, 'stress', stress, resDir)
    b = readGTPlusExportImageSeries_Squeeze(stressDir, 120);
    fb = readGTPlusExportImageSeries_Squeeze(stressDir, 115);
    fb = fb/100;
    fb_pde = readGTPlusExportImageSeries_Squeeze(stressDir, 124);
    fb_pde = fb_pde/100;
    aif_stress = analyze75read(fullfile(stressDir, 'DebugOutput', 'aif_cin_all.hdr'));
    aif_stress_without_R2Star = analyze75read(fullfile(stressDir, 'DebugOutput', 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
    
    [slope_stress, timeToPeak_stress, peakTime_stress, areaUnderCurve_stress, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(aif_stress(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 0, 'Feature detection of AIF LV signal');
    
    slope_stress2 = aif_stress(round(peakTime_rest):end);
    [slope, timeToPeak, peakTime, areaUnderCurve, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(-slope_stress2(:), sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 0, 'Feature detection of AIF LV signal');   
    valley_stress = round(peakTime_stress) + peakTime - 0.5;
    
    foot_stress = round(peakTime_stress - timeToPeak_stress);   
    v = mean(aif_stress_without_R2Star(1:foot_stress-1));
    aif_stress_without_R2Star = aif_stress_without_R2Star - v;
    
    dset = ismrmrd.Dataset(stressDat, 'dataset');
    stress_header = ismrmrd.xml.deserialize(dset.readxml());
    acq = dset.readAcquisition(1);
    stress_t = double(acq.head.acquisition_time_stamp);
    
    if(normal)
        if(isempty(Ki_map_rest_normal) || sfa(2)==size(Ki_map_rest_normal, 2) )
            Ki_map_rest_normal = cat(3, Ki_map_rest_normal, fa(:,:,2,3));
            Ki_map_stress_normal = cat(3, Ki_map_stress_normal, fb(:,:,2,3));

            Flow_map_rest_normal = cat(3, Flow_map_rest_normal, fa_pde(:,:,2,end));
            Flow_map_stress_normal = cat(3, Flow_map_stress_normal, fb_pde(:,:,2,end));
            
            normal_cases_maps = [normal_cases_maps; {subdir{ii}}];
        end
    else
        if(isempty(Ki_map_rest_sick) || sfa(2)==size(Ki_map_rest_sick, 2) )
            Ki_map_rest_sick = cat(3, Ki_map_rest_sick, fa(:,:,2,3));
            Ki_map_stress_sick = cat(3, Ki_map_stress_sick, fb(:,:,2,3));

            Flow_map_rest_sick = cat(3, Flow_map_rest_sick, fa_pde(:,:,2,end));
            Flow_map_stress_sick = cat(3, Flow_map_stress_sick, fb_pde(:,:,2,end));
            
            sickCases = [sickCases; {subdir{ii}}];
        end
    end
    
    disp('--------------------------------------------------------------');

    a = permute(a, [2 1 3]);
    b = permute(b, [2 1 3]);

    [namesDcm, numDcm] = findFILE(fullfile(home, subdir{ii}), '*.dcm');
    if(numDcm==0)
        [namesDcm, numDcm] = findFILE(fullfile(home, subdir{ii}), '*.IMA');
    end
    
    if(numDcm==2)
        info = dicominfo(namesDcm{1});
        t1 = info.AcquisitionTime;

        info2 = dicominfo(namesDcm{2});
        t2 = info2.AcquisitionTime;

        t1s = getSecondsFromTimeString(t1)
        t2s = getSecondsFromTimeString(t2)

        tDiff = (rest_t - stress_t)*2.5/1000/60;

        if(str2double(t1)>str2double(t2))
         % t2 is pre-stress
         pre_stress = namesDcm{2};
         pre_rest = namesDcm{1};
        else
         % t1 is  pre-sterss
         pre_stress = namesDcm{1};
         pre_rest = namesDcm{2};
        end         

        pre_stress_t1map = dicomread(pre_stress); pre_stress_t1map = double(pre_stress_t1map);
        pre_rest_t1map = dicomread(pre_rest); pre_rest_t1map = double(pre_rest_t1map);

        if(~isFileExist(fullfile(home, subdir{ii}, 'pre_stress_roi.mat')))
             figure; imagescn(pre_stress_t1map);
             figure; imagescn(pre_rest_t1map, [0 1400]);

             pause    
             closeall
        else
            r = roi_timeseries(pre_stress_t1map, fullfile(home, subdir{ii}, 'pre_stress_roi.mat'), 1, 1);
            stress_t1_blood = [stress_t1_blood; r.m tDiff];

            r = roi_timeseries(pre_stress_t1map, fullfile(home, subdir{ii}, 'pre_stress_roi.mat'), 2, 1);
            stress_t1_myo = [stress_t1_myo; r.m tDiff];

            r = roi_timeseries(pre_rest_t1map, fullfile(home, subdir{ii}, 'pre_rest_roi.mat'), 1, 1);
            rest_t1_blood = [rest_t1_blood; r.m tDiff ];

            r = roi_timeseries(pre_rest_t1map, fullfile(home, subdir{ii}, 'pre_rest_roi.mat'), 2, 1);
            rest_t1_myo = [rest_t1_myo; r.m tDiff];
        end
    end
    
    if(numDcm==1)
        info = dicominfo(namesDcm{1});
    end
    
    if(isFileExist(fullfile(home, subdir{ii}, 'rest_roi1.mat')))
        r1 = roi_timeseries(fa(:,:,1,3), fullfile(home, subdir{ii}, 'rest_roi1.mat'), 1, 1); 
        r2 = roi_timeseries(fa(:,:,2,3), fullfile(home, subdir{ii}, 'rest_roi2.mat'), 1, 1); 
        r3 = roi_timeseries(fa(:,:,3,3), fullfile(home, subdir{ii}, 'rest_roi3.mat'), 1, 1); 
        
        rest_pt = [rest_pt; r1.data(:); r2.data(:); r3.data(:)];
        
        r1_pde = roi_timeseries(fa_pde(:,:,1,end), fullfile(home, subdir{ii}, 'rest_roi1.mat'), 1, 1); 
        r2_pde = roi_timeseries(fa_pde(:,:,2,end), fullfile(home, subdir{ii}, 'rest_roi2.mat'), 1, 1); 
        r3_pde = roi_timeseries(fa_pde(:,:,3,end), fullfile(home, subdir{ii}, 'rest_roi3.mat'), 1, 1); 

        rest_pde_pt = [rest_pde_pt; r1_pde.data(:); r2_pde.data(:); r3_pde.data(:)];
        
        age = str2double(info.PatientAge(2:end-1));
        
        ResultRecordRest = [ResultRecordRest; {ind subdir{ii} 'rest' 1 r1.m r1_pde.m max(aif_rest) aif_rest timeToPeak_rest peakTime_rest valley_rest info.PatientSex age EF max(aif_rest_without_R2Star) aif_rest_without_R2Star}]; ind = ind + 1;
        ResultRecordRest = [ResultRecordRest; {ind subdir{ii} 'rest' 2 r2.m r2_pde.m max(aif_rest) aif_rest timeToPeak_rest peakTime_rest valley_rest info.PatientSex age EF max(aif_rest_without_R2Star) aif_rest_without_R2Star }]; ind = ind + 1;
        ResultRecordRest = [ResultRecordRest; {ind subdir{ii} 'rest' 3 r3.m r3_pde.m max(aif_rest) aif_rest timeToPeak_rest peakTime_rest valley_rest info.PatientSex age EF max(aif_rest_without_R2Star) aif_rest_without_R2Star}]; ind = ind + 1;

        s1 = roi_timeseries(fb(:,:,1,3), fullfile(home, subdir{ii}, 'stress_roi1.mat'), 1, 1); 
        s2 = roi_timeseries(fb(:,:,2,3), fullfile(home, subdir{ii}, 'stress_roi2.mat'), 1, 1); 
        s3 = roi_timeseries(fb(:,:,3,3), fullfile(home, subdir{ii}, 'stress_roi3.mat'), 1, 1); 
        stress_pt = [stress_pt; s1.data(:); s2.data(:); s3.data(:)];
        
        if(~normal & ~infarct)
            % ischemia
            s1_ischemia = roi_timeseries(fb(:,:,1,3), fullfile(home, subdir{ii}, 'stress_roi1.mat'), 2, 1); 
            s2_ischemia = roi_timeseries(fb(:,:,2,3), fullfile(home, subdir{ii}, 'stress_roi2.mat'), 2, 1); 
            s3_ischemia = roi_timeseries(fb(:,:,3,3), fullfile(home, subdir{ii}, 'stress_roi3.mat'), 2, 1); 
            stress_pt = [stress_pt; s1_ischemia.data(:); s2_ischemia.data(:); s3_ischemia.data(:)];
        end
         
        s1_pde = roi_timeseries(fb_pde(:,:,1,end), fullfile(home, subdir{ii}, 'stress_roi1.mat'), 1, 1); 
        s2_pde = roi_timeseries(fb_pde(:,:,2,end), fullfile(home, subdir{ii}, 'stress_roi2.mat'), 1, 1); 
        s3_pde = roi_timeseries(fb_pde(:,:,3,end), fullfile(home, subdir{ii}, 'stress_roi3.mat'), 1, 1); 
        stress_pde_pt = [stress_pde_pt; s1_pde.data(:); s2_pde.data(:); s3_pde.data(:)];

        if(~normal & ~infarct)
            % ischemia
            s1_pde_ischemia = roi_timeseries(fb_pde(:,:,1,end), fullfile(home, subdir{ii}, 'stress_roi1.mat'), 2, 1); 
            s2_pde_ischemia = roi_timeseries(fb_pde(:,:,2,end), fullfile(home, subdir{ii}, 'stress_roi2.mat'), 2, 1); 
            s3_pde_ischemia = roi_timeseries(fb_pde(:,:,3,end), fullfile(home, subdir{ii}, 'stress_roi3.mat'), 2, 1); 
            stress_pde_pt = [stress_pde_pt; s1_pde_ischemia.data(:); s2_pde_ischemia.data(:); s3_pde_ischemia.data(:)];
        end
        
        ResultRecordStress = [ResultRecordStress; {indS subdir{ii} 'stress' 1 s1.m s1_pde.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS = indS + 1;
        ResultRecordStress = [ResultRecordStress; {indS subdir{ii} 'stress' 2 s2.m s2_pde.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS = indS + 1;
        ResultRecordStress = [ResultRecordStress; {indS subdir{ii} 'stress' 3 s3.m s3_pde.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS = indS + 1;
        
        if(~normal & ~infarct)
            ResultRecordStressIschemia = [ResultRecordStressIschemia; {indS_ischemia subdir{ii} 'stress' 1 s1.m s1_pde.m s1_ischemia.m s1_pde_ischemia.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS_ischemia = indS_ischemia + 1;
            ResultRecordStressIschemia = [ResultRecordStressIschemia; {indS_ischemia subdir{ii} 'stress' 1 s2.m s2_pde.m s2_ischemia.m s2_pde_ischemia.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS_ischemia = indS_ischemia + 1;
            ResultRecordStressIschemia = [ResultRecordStressIschemia; {indS_ischemia subdir{ii} 'stress' 1 s3.m s3_pde.m s3_ischemia.m s3_pde_ischemia.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS_ischemia = indS_ischemia + 1;
        end
        
        noroi = 0;
    else
        disp(['Do not have roi files :' subdir{ii}]);
        noroi = 1;
    end
    
    if(plotFlag & noroi)
        
        disp(['normal = ' num2str(normal) ' - infarct = ' num2str(infarct)]);
        
        h = figure; 
        imagescn(cat(4, a, b), [], [], 25)
        saveas(h, fullfile(home, subdir{ii}, 'AIF_FIG'), 'fig')

        % figure; imagescn(cat(3, fa(:,:,:,1), fb(:,:,:,1)), [0 6]); PerfColorMap;
        figure; imagescn(fa(:,:,1,3), [0 6]); PerfColorMap;
        figure; imagescn(fa(:,:,2,3), [0 6]); PerfColorMap;
        figure; imagescn(fa(:,:,3,3), [0 6]); PerfColorMap;

        figure; imagescn(fb(:,:,1,3), [0 6]); PerfColorMap;
        figure; imagescn(fb(:,:,2,3), [0 6]); PerfColorMap;
        figure; imagescn(fb(:,:,3,3), [0 6]); PerfColorMap;

        figure; imagescn(fa_pde(:,:,1,end), [0 6]); PerfColorMap;
        figure; imagescn(fa_pde(:,:,2,end), [0 6]); PerfColorMap;
        figure; imagescn(fa_pde(:,:,3,end), [0 6]); PerfColorMap;

        figure; imagescn(fb_pde(:,:,1,end), [0 6]); PerfColorMap;
        figure; imagescn(fb_pde(:,:,2,end), [0 6]); PerfColorMap;
        figure; imagescn(fb_pde(:,:,3,end), [0 6]); PerfColorMap;
        
        pause    
        closeall
    end       
end

cd(home)
save BigBolusResult_water_r1_r2 ResultRecordRest ResultRecordStress normalCases Ki_map_rest_normal Ki_map_stress_normal Flow_map_rest_normal Flow_map_stress_normal sickCases Ki_map_rest_sick Ki_map_stress_sick Flow_map_rest_sick Flow_map_stress_sick rest_pt stress_pt rest_pde_pt stress_pde_pt rest_t1_myo rest_t1_blood stress_t1_myo stress_t1_blood normal_cases_maps

cd(home)
save BigBolusResult_water_reduced_r1_r2 ResultRecordRest ResultRecordStress normalCases Ki_map_rest_normal Ki_map_stress_normal Flow_map_rest_normal Flow_map_stress_normal sickCases Ki_map_rest_sick Ki_map_stress_sick Flow_map_rest_sick Flow_map_stress_sick rest_pt stress_pt rest_pde_pt stress_pde_pt rest_t1_myo rest_t1_blood stress_t1_myo stress_t1_blood normal_cases_maps

cd(home)
save BigBolusResult_water_increased_r1_r2 ResultRecordRest ResultRecordStress normalCases Ki_map_rest_normal Ki_map_stress_normal Flow_map_rest_normal Flow_map_stress_normal sickCases Ki_map_rest_sick Ki_map_stress_sick Flow_map_rest_sick Flow_map_stress_sick rest_pt stress_pt rest_pde_pt stress_pde_pt rest_t1_myo rest_t1_blood stress_t1_myo stress_t1_blood normal_cases_maps


normalCases
sickCases
size(sickCases)
size(normalCases)

%% rest 

cd(home)
a = load('BigBolusResult_water_r1_r2');
b = load('BigBolusResult_water_reduced_r1_r2');
b = load('BigBolusResult_water_increased_r1_r2');

% ResultRecordStress = [ResultRecordStress; {indS subdir{ii} 'stress' 1 s1.m s1_pde.m max(aif_stress) aif_stress timeToPeak_stress peakTime_stress valley_stress info.PatientSex age EF max(aif_stress_without_R2Star) aif_stress_without_R2Star}]; indS = indS + 1;

N = size(a.ResultRecordRest, 1)

rest_flow_a = zeros(N, 1);
rest_flow_b = zeros(N, 1);

aif_peak_a = zeros(N, 1);
aif_peak_b = zeros(N, 1);

aif_duration_a = zeros(N, 1);
aif_duration_b = zeros(N, 1);

for ii=1:N
    rest_flow_a(ii) = a.ResultRecordRest{ii, 5};
    aif_peak_a(ii) = a.ResultRecordRest{ii, 7};
    aif_duration_a(ii) = a.ResultRecordRest{ii, 11} - a.ResultRecordRest{ii, 10} + a.ResultRecordRest{ii, 9};
    
    rest_flow_b(ii) = b.ResultRecordRest{ii, 5};
    aif_peak_b(ii) = b.ResultRecordRest{ii, 7};
    aif_duration_b(ii) = b.ResultRecordRest{ii, 11} - b.ResultRecordRest{ii, 10} + b.ResultRecordRest{ii, 9};
end

legendStr1 = 'water r1 and r2'
legendStr2 = 'reduced r1 and r1'
legendStr2 = 'increased r1 and r1'

plot_x_y_with_linefit(rest_flow_a, rest_flow_b, ['x : flow, ' legendStr1], ['y : flow, ' legendStr2], 'rest flow')

plot_x_y_with_linefit(aif_peak_a, aif_peak_b, ['x : aif peak, ' legendStr1], ['y : aif peak, ' legendStr2], 'rest, aif peak Gd')

plot_x_y_with_linefit(aif_duration_a, aif_duration_b, ['x : aif duration, ' legendStr1], ['y : aif duration, ' legendStr2], 'rest, aif duration in seconds')

%% stress
closeall

N = size(a.ResultRecordStress, 1)

flow_a = zeros(N, 1);
flow_b = zeros(N, 1);

aif_peak_a = zeros(N, 1);
aif_peak_b = zeros(N, 1);

aif_duration_a = zeros(N, 1);
aif_duration_b = zeros(N, 1);

for ii=1:N
    flow_a(ii) = a.ResultRecordStress{ii, 5};
    aif_peak_a(ii) = a.ResultRecordStress{ii, 7};
    aif_duration_a(ii) = a.ResultRecordStress{ii, 11} - a.ResultRecordStress{ii, 10} + a.ResultRecordStress{ii, 9};
    
    flow_b(ii) = b.ResultRecordStress{ii, 5};
    aif_peak_b(ii) = b.ResultRecordStress{ii, 7};
    aif_duration_b(ii) = b.ResultRecordStress{ii, 11} - b.ResultRecordStress{ii, 10} + b.ResultRecordStress{ii, 9};
end

ind = [1:3 7:45];
ind = [1:30];
ind = [1:N];

legendStr1 = 'water r1 and r2'
legendStr2 = 'reduced r1 and r1'
legendStr2 = 'increased r1 and r1'

plot_x_y_with_linefit(flow_a(ind), flow_b(ind), ['x : flow, ' legendStr1], ['y : flow, ' legendStr2], 'stress flow')

plot_x_y_with_linefit(aif_peak_a(ind), aif_peak_b(ind), ['x : aif peak, ' legendStr1], ['y : aif peak, ' legendStr2], 'stress, aif peak Gd')

plot_x_y_with_linefit(aif_duration_a(ind), aif_duration_b(ind), ['x : aif duration, ' legendStr1], ['y : aif duration, ' legendStr2], 'stress, aif duration in seconds')


figure
hold on
plot(flow_a(ind), flow_b(ind), '.');
hold off
box on
grid on
xlabel('flow, water r1 and r2')
ylabel('flow, reduced r1 and r1');
axis equal
title('stress flow')

figure
hold on
plot(aif_peak_a(ind), aif_peak_b(ind), '.');
hold off
box on
grid on
xlabel('aif peak, water r1 and r2')
ylabel('aif peak, reduced r1 and r1');
axis equal
title('stress, aif peak Gd')

figure
hold on
plot(aif_duration_a(ind), aif_duration_b(ind), '.');
hold off
box on
grid on
xlabel('aif duration, water r1 and r2')
ylabel('aif duration, reduced r1 and r1');
axis equal
title('stress, aif duration in seconds')
