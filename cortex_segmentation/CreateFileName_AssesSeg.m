function [Global_dsc_resultfile, Global_PVs_dsc_resultfile, Global_Local_PVs_dsc_resultfile] = ...
    CreateFileName_AssesSeg(home, Labeled_slices, classnumber, GlobalResultDir, LocalResultDir)
% create file name

% ========================================================================================= %
% Global_dsc
GMM_PVs_flag = 0;
GMM_PVs_flag_Local = 0;

if (classnumber==5)
    Four_classes_flag = 0;
    Five_classes_flag = 1;
end
if (classnumber==4)
    Four_classes_flag = 1;
    Five_classes_flag = 0;
end
            
Global_flag = 1; 
Local_flag = 0;

Atlas_flag_Global = 0;
Atlas_flag_Local = 0;

Kmeans_flag_Global = 1;
Kmeans_flag_Local = 1;

MRF_flag_Global = 1;
MRF_flag_Local = 1;

Save_flag = 1;

partsNumber = Labeled_slices{3};

prefix = CreatePrefix(home, Global_flag, Local_flag, Atlas_flag_Global, Atlas_flag_Local,...
    Kmeans_flag_Global, Kmeans_flag_Local, MRF_flag_Global, MRF_flag_Local, Four_classes_flag, ...
    Five_classes_flag, partsNumber, GMM_PVs_flag, GMM_PVs_flag_Local);
prefix;

if (classnumber==5)
    Global_dsc_resultfile = [Labeled_slices{1} prefix '_segResult_5classes_MRF_uint8.hdr'];
end
if (classnumber==4)
    Global_dsc_resultfile = [Labeled_slices{1} prefix '_segResult_MRF_4classes.hdr'];
end
Global_dsc_resultfile = fullfile(GlobalResultDir, Global_dsc_resultfile);

% ========================================================================================= %
% Global_PVs_dsc
GMM_PVs_flag = 0;
GMM_PVs_flag_Local = 0;

if (classnumber==5)
    Four_classes_flag = 0;
    Five_classes_flag = 1;
end
if (classnumber==4)
    Four_classes_flag = 1;
    Five_classes_flag = 0;
end
            
Global_flag = 1; 
Local_flag = 0;

Atlas_flag_Global = 0;
Atlas_flag_Local = 0;

Kmeans_flag_Global = 1;
Kmeans_flag_Local = 1;

MRF_flag_Global = 1;
MRF_flag_Local = 1;

Save_flag = 1;

partsNumber = Labeled_slices{3};

prefix = CreatePrefix(home, Global_flag, Local_flag, Atlas_flag_Global, Atlas_flag_Local,...
    Kmeans_flag_Global, Kmeans_flag_Local, MRF_flag_Global, MRF_flag_Local, Four_classes_flag, ...
    Five_classes_flag, partsNumber, GMM_PVs_flag, GMM_PVs_flag_Local);
prefix;

if (classnumber==5)
%     Global_PVs_dsc_resultfile = [Labeled_slices{1} prefix '_segResult_5classes_MRF_uint8.hdr'];
    Global_PVs_dsc_resultfile = [Labeled_slices{1} prefix '_segResult_5classes_MRF_uint8_PVs.hdr'];
end
if (classnumber==4)
%     Global_PVs_dsc_resultfile = [Labeled_slices{1} prefix '_segResult_MRF_4classes.hdr'];
    Global_PVs_dsc_resultfile = [Labeled_slices{1} prefix '_segResult_MRF_4classes_PVs.hdr'];
end
Global_PVs_dsc_resultfile = fullfile(GlobalResultDir, Global_PVs_dsc_resultfile);

% ========================================================================================= %
% Global_Local_PVs_dsc_resultfile
GMM_PVs_flag = 0;
GMM_PVs_flag_Local = 0;

if (classnumber==5)
    Four_classes_flag = 0;
    Five_classes_flag = 1;
end
if (classnumber==4)
    Four_classes_flag = 1;
    Five_classes_flag = 0;
end
            
Global_flag = 1; 
Local_flag = 1;

Atlas_flag_Global = 0;
Atlas_flag_Local = 0;

Kmeans_flag_Global = 1;
Kmeans_flag_Local = 1;

MRF_flag_Global = 1;
MRF_flag_Local = 1;

Save_flag = 1;

partsNumber = Labeled_slices{3};

prefix = CreatePrefix(home, Global_flag, Local_flag, Atlas_flag_Global, Atlas_flag_Local,...
    Kmeans_flag_Global, Kmeans_flag_Local, MRF_flag_Global, MRF_flag_Local, Four_classes_flag, ...
    Five_classes_flag, partsNumber, GMM_PVs_flag, GMM_PVs_flag_Local);
prefix;

if (classnumber==5)
    Global_Local_PVs_dsc_resultfile = [Labeled_slices{1} prefix '_newGlobal_SegResult_PVs.hdr'];
end
if (classnumber==4)
    Global_Local_PVs_dsc_resultfile = [Labeled_slices{1} prefix '_newGlobal_SegResult_PVs.hdr'];
end
Global_Local_PVs_dsc_resultfile = fullfile(LocalResultDir, Global_Local_PVs_dsc_resultfile);

return;