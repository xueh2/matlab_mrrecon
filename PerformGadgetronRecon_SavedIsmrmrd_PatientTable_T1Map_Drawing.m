
function PerformGadgetronRecon_SavedIsmrmrd_PatientTable_T1Map_Drawing(PerfTable, resDir, contourDir, T1Dir, SelectedType, splenic_column, selected_column, stress_column, rest_column, report_column)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_T1Map_Drawing(PerfTable, resDir, contourDir, T1Dir, SelectedType, splenic_column, selected_column, stress_column, rest_column, report_column)

if(nargin<4)
    SelectedType = [];
end

if(exist(contourDir)~=7)
    mkdir(contourDir);
end

[subdirs, numT1] = FindSubDirs(T1Dir);

num = size(PerfTable, 1)-1;
for n=1:num
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column} ' - ' PerfTable{n+1, splenic_column} ' - ' PerfTable{n+1, selected_column}]); 
    disp(['==================================================================']);  
    
    splenic_cut_off = PerfTable{n+1, splenic_column};
    selected = PerfTable{n+1, selected_column};
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    report = PerfTable{n+1, report_column};
    
    process = 0;
    if(isempty(SelectedType)~=1)
        for ii=1:numel(SelectedType)
            if( ~isempty(strfind(selected, SelectedType{ii})) )
                process = 1;
                break;
            end
        end
    else
        process = 1;
    end
    
    if( (strcmp(splenic_cut_off, 'Yes') ~= 1) & (isempty( strfind(splenic_cut_off, 'yes') ) == 1) )
        continue;
    end
    
    if(~process)
        continue;
    end
    
    disp(['Selected Type : ' selected]); 
    
%     onlyReview = 0;
%     [h_flow_stress, h_flow_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  stressCase, restCase, [0 6], onlyReview);        
    
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(stressCase);
    
    roiDir = fullfile(contourDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' scannerID '_' patientID '_' studyID '_' study_dates '_ROI'])
    
    if(isFileExist(roiDir))
        dir(roiDir);        
    else    
        mkdir(roiDir)
    end
    cd(roiDir)
    
    findT1 = 0;
    for kk=1:numT1
        if(~isempty(strfind(report, subdirs{kk})))
            findT1 = 1;
            break;
        end
    end
    
    disp([report ' - ' subdirs{kk}]);
    
    if(findT1 & ~isFileExist(fullfile(roiDir, 'T1Maps.mat')))
        copyfile(fullfile(T1Dir, subdirs{kk}, 'pre'), fullfile(roiDir, 'pre'));
        copyfile(fullfile(T1Dir, subdirs{kk}, 'post'), fullfile(roiDir, 'post'));
        
        [pre_names, num_pre] = findFILE(fullfile(roiDir, 'pre'), '*.dcm');
        [post_names, num_post] = findFILE(fullfile(roiDir, 'post'), '*.dcm');
        
        preT1 = [];
        for kk=1:num_pre
            preT1(:,:,kk) = dicomread(pre_names{kk});
        end
        
        postT1 = [];
        for kk=1:num_post
            postT1(:,:,kk) = dicomread(post_names{kk});
        end
        
        figure; imagescn(preT1, [0 2000], [], 16);
        figure; imagescn(postT1, [0 600], [], 16);
        
        save(fullfile(roiDir, 'T1Maps.mat'), 'preT1', 'postT1');
    else
        v = load(fullfile(roiDir, 'T1Maps.mat'));
        preT1 = v.preT1;
        postT1 = v.postT1;
        
        figure; imagescn(preT1, [0 2000], [], 16);
        figure; imagescn(postT1, [0 600], [], 16);

    end
    
    pause;
    closeall;
end
