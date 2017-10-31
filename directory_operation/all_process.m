

% home directory stores all subjects' subdirectory...
home = 'J:\backup\from_pc\data\mra\term_infants';
[subdirs, num] = FindAllDirectory(home)

tofDir = 'sense';
t2Dir = 't2';
% ------------------------------------ PART 1 : vessel extraction -------%
% sense: store a tof image
% t2: store a t2 image
% tofROI: store cutted tof image and extractiong results

tofROI = 'tofROI';
CreateDirs(home, tofROI); %%

tofCut = 'tofCut';
CreateDirs(home, tofCut);

tofResampled = 'tofResampled';
CreateDirs(home, tofResampled);

extractionResults = 'extractionResults';
CreateDirs(home, extractionResults);

resolutions = [0.2 0.2 0.2];

background = 35;

for i = 1:num
    
    disp(currentDir);
    disp('=====================================================');

    currentDir = [home '\' subdirs{i}];
    
    % cut the image and save the results
    CutImageRun(currentDir, tofDir, tofCut, background); %% ---> output includes a txt file to show leftup and rightdown limits
   
    ResampleRun(currentDir, tofCut, tofResampled, resolutions);
end

for i = 1:num
    currentDir = [home '\' subdirs{i}];
        
    disp(currentDir);
    disp('=====================================================');

    % performing extraction
    VesselExtractionRun(currentDir, tofResampled, extractionResults);
end

% ------------------------------------ PART 2 : tof-t2 registration -------%
% registration dof files between tof-t2 are stored in registrationresults
registrationresults_t2_tof = 'registrationresults_t2_tof';
CreateDirs(home, registrationresults_t2_tof);

parameterfile_rreg = '';
% tof ---> t2
for i = 1:num
    currentDir = [home '\' subdirs{i}];
    disp(currentDir);
    disp('=====================================================');
    RigidRegistrationRun(home, t2Dir, tofDir, registrationresults_t2_tof, parameterfile_rreg);
end

% ------------------------------------ PART 3 : extracted vessels transformed into t2 coordinates -------%
vesselTransformed_in_t2 = 'vesselTransformed_in_t2';
CreateDirs(home, vesselTransformed_in_t2);

for i = 1:num
    currentDir = [home '\' subdirs{i}];
    disp(currentDir);
    disp('=====================================================');
    VesselTransformRun(home, extractionResults, tofCut, registrationresults_t2_tof, vesselTransformed_in_t2);
end