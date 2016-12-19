
function LGE = PerformGadgetronRecon_FindLGE_PerfusionCase(dataDir, resDir, perf_case, showLGE)
% LGE = PerformGadgetronRecon_FindLGE_PerfusionCase(dataDir, resDir, perf_case)
% LGE = PerformGadgetronRecon_FindLGE_PerfusionCase('I:\BARTS', 'I:\ReconResults\BARTS', perf_cases)
% given a perfusion case, find corresponding LGE cases

if(nargin<4)
    showLGE = 0;
end

[configName, scannerID_rest, patientID_rest, studyID_rest, measurementID_rest, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(perf_case);

studyDir = fullfile(resDir, study_dates);

[subdirs, num] = FindSubDirs(studyDir);

LGE = [];

for ii=1:num
    
    v1 = strfind(subdirs{ii}, 'LGE_MOCO');
    v2 = strfind(subdirs{ii}, patientID_rest);
    
    if(~isempty(v1) & ~isempty(v2))
        
        v3 = strfind(subdirs{ii}, 'dicom');
       if(isempty(v3))
           a = readGTPlusExportImageSeries_Squeeze(fullfile(resDir, study_dates, subdirs{ii}), 111);
           LGE = [LGE; {a}];
           if(showLGE)
               figure; imagescn(a, [2048 1.4*4086], [], 10);
           end
       end
    end    
end
