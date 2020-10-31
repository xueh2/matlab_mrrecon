
cd(Upperhome);
[subdirs, num] = FindAllDirectory(Upperhome)
% ============================================================== %
% common setting
templateDir = './warped_templates/';
brainMaskfile = [templateDir 'warped_brainmask_nostem.hdr'];

imageDir = './N3Brain/';
imagefile = [imageDir 'noStemBrain_N3.hdr'];

Global_SegResult = './GlobalSeg_N3/';
Local_SegResult = './LocalSeg_N3/';
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
for tt = 10:num
    tt = 3;
    home = [Upperhome '\' subdirs{tt}];
        
    disp(home);
    disp('=====================================================');

    run_Neonatal_Segmentation_allConfigurations

end

% ====================================================================== %