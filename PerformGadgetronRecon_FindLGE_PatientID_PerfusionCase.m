
function LGE = PerformGadgetronRecon_FindLGE_PatientID_PerfusionCase(dataDir, resDir, patientID, study_date, showLGE)
% LGE = PerformGadgetronRecon_FindLGE_PatientID_PerfusionCase(dataDir, resDir, patientID, study_date, showLGE)
% LGE = PerformGadgetronRecon_FindLGE_PerfusionCase('I:\BARTS', 'I:\ReconResults\BARTS', patientID, study_date)
% given a perfusion case, find corresponding LGE cases

if(nargin<5)
    showLGE = 0;
end

studyDir = fullfile(resDir, study_date);

[subdirs, num] = FindSubDirs(studyDir);

LGE = [];

for ii=1:num
    
    v1 = strfind(subdirs{ii}, 'LGE_MOCO');
    v2 = strfind(subdirs{ii}, patientID);
    
    if(~isempty(v1) & ~isempty(v2))
        
        v3 = strfind(subdirs{ii}, 'dicom');
       if(isempty(v3))
           a = readGTPlusExportImageSeries_Squeeze(fullfile(resDir, study_date, subdirs{ii}), 111);
           LGE = [LGE; {a}];
           if(showLGE)
               figure; imagescn(a, [2048 1.4*4086], [], 10);
           end
       end
    end    
end
