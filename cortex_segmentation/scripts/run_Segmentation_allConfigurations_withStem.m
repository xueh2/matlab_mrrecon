
cd(Upperhome);
[subdirs, num] = FindAllDirectory(Upperhome)
% ============================================================== %
% common setting
templateDir = './warped_templates/';
brainMaskfile = [templateDir 'warped_brainmask.hdr'];

imageDir = './N3Brain/';
imagefile = [imageDir 'withStemBrain_N3.hdr'];

Global_SegResult = './GlobalSeg_withStem_N3/';
Local_SegResult = './LocalSeg_withStem_N3/';
% algorithm parameters

TryNumber = 10;
trynumber = TryNumber;

Kmeans_InitialCentres = [];
KmeansDir = './Kmeans_InitialCentres/';

Low_Thres = 0.8;
volumeThres = 8;
partsNumber = 7;
initialCentres = [];

% ============================================================== %
for i = 1:num
    home = [Upperhome '\' subdirs{i}];
        
    disp(home);
    disp('=====================================================');

    run_Neonatal_Segmentation_allConfigurations

end

% ====================================================================== %