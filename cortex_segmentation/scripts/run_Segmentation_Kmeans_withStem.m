
cd(Upperhome);
[subdirs, num] = FindAllDirectory(Upperhome)
% ============================================================== %
% common setting
templateDir = './warped_templates/';
brainMaskfile = [templateDir 'warped_brainmask_nostem.hdr'];

imageDir = './N3Brain/';
imagefile = [imageDir 'withStemBrain_N3.hdr'];

% algorithm parameters
TryNumber = 10;
trynumber = TryNumber;
Kmeans_InitialCentres = [];
KmeansDir = './Kmeans_InitialCentres/';
initialCentres = [];
% ============================================================== %
for tt = 1:num
    home = [Upperhome '\' subdirs{tt}];
        
    disp(home);
    disp('=====================================================');
    cd(home)

    classnumber = 4;
    run_KmeansClustering;

    classnumber = 5;
    run_KmeansClustering;
end

% ====================================================================== %