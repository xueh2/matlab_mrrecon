
function perf_create_CMR_View_NN_data(caseName, resDir, in_CMR_View)
% perf_create_CMR_View_NN_data(caseName, resDir)

if (nargin < 3)
    in_CMR_View = 0
end

[configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(caseName);

caseDir = fullfile(resDir, study_date, caseName);
debugDir = fullfile(resDir, study_date, caseName, 'DebugOutput');
trainingDataDir = fullfile(resDir, study_date, caseName, 'Training_Data');
trainingCaseDir = trainingDataDir;

[aif_scan_geometry_info, scan_geometry_info] = read_in_GT_Perf_DebugOutput_scan_geometry_info(caseDir);

slc = size(scan_geometry_info.slice_dir,1);
disp(['Total ' num2str(slc) ' is found ...']);

good_slices = [];        
for s=1:size(scan_geometry_info.slice_dir,1)

    slice_dir = scan_geometry_info.slice_dir(s,:)

    diff = norm(aif_scan_geometry_info.aif_slice_dir - slice_dir);
    if(diff<0.1)
        good_slices = [good_slices s];
    end
end
good_slices

fmap = readNPY(fullfile(trainingCaseDir, 'fmap.npy'));
fmap_resized = readNPY(fullfile(trainingCaseDir, 'fmap_resized.npy'));
fmap_resized_training = readNPY(fullfile(trainingCaseDir, 'fmap_resized_training.npy'));
SLC = size(fmap, 3);

if(in_CMR_View)
    fmap_norm = permute(fmap, [2 1 3]);
    fmap_resized_norm = permute(fmap_resized, [2 1 3]);
    fmap_resized_training_norm = permute(fmap_resized_training, [2 1 3]);
else
    fmap_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fmap, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    fmap_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fmap_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    fmap_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fmap_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

fSD = readNPY(fullfile(trainingCaseDir, 'fSD.npy'));
fSD_resized = readNPY(fullfile(trainingCaseDir, 'fSD_resized.npy'));
fSD_resized_training = readNPY(fullfile(trainingCaseDir, 'fSD_resized_training.npy'));

if(in_CMR_View)
    fSD_norm = permute(fSD, [2 1 3]);
    fSD_resized_norm = permute(fSD_resized, [2 1 3]);
    fSD_resized_training_norm = permute(fSD_resized_training, [2 1 3]);
else
    fSD_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fSD, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    fSD_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fSD_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    fSD_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fSD_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

PS_map = readNPY(fullfile(trainingCaseDir, 'PS_map.npy'));
PS_map_resized = readNPY(fullfile(trainingCaseDir, 'PS_map_resized.npy'));
PS_map_resized_training = readNPY(fullfile(trainingCaseDir, 'PS_map_resized_training.npy'));

if(in_CMR_View)
    PS_map_norm = permute(PS_map, [2 1 3]);
    PS_map_resized_norm = permute(PS_map_resized, [2 1 3]);
    PS_map_resized_training_norm = permute(PS_map_resized_training, [2 1 3]);
else
    PS_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(PS_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    PS_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(PS_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    PS_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(PS_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

vb_map = readNPY(fullfile(trainingCaseDir, 'vb_map.npy'));
vb_map_resized = readNPY(fullfile(trainingCaseDir, 'vb_map_resized.npy'));
vb_map_resized_training = readNPY(fullfile(trainingCaseDir, 'vb_map_resized_training.npy'));

if(in_CMR_View)
    vb_map_norm = permute(vb_map, [2 1 3]);
    vb_map_resized_norm = permute(vb_map_resized, [2 1 3]);
    vb_map_resized_training_norm = permute(vb_map_resized_training, [2 1 3]);
else
    vb_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(vb_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    vb_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(vb_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    vb_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(vb_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

visf_map = readNPY(fullfile(trainingCaseDir, 'visf_map.npy'));
visf_map_resized = readNPY(fullfile(trainingCaseDir, 'visf_map_resized.npy'));
visf_map_resized_training = readNPY(fullfile(trainingCaseDir, 'visf_map_resized_training.npy'));

if(in_CMR_View)
    visf_map_norm = permute(visf_map, [2 1 3]);
    visf_map_resized_norm = permute(visf_map_resized, [2 1 3]);
    visf_map_resized_training_norm = permute(visf_map_resized_training, [2 1 3]);
else
    visf_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(visf_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    visf_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(visf_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    visf_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(visf_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

Tc_map = readNPY(fullfile(trainingCaseDir, 'Tc_map.npy'));
Tc_map_resized = readNPY(fullfile(trainingCaseDir, 'Tc_map_resized.npy'));
Tc_map_resized_training = readNPY(fullfile(trainingCaseDir, 'Tc_map_resized_training.npy'));

if(in_CMR_View)
    Tc_map_norm = permute(Tc_map, [2 1 3]);
    Tc_map_resized_norm = permute(Tc_map_resized, [2 1 3]);
    Tc_map_resized_training_norm = permute(Tc_map_resized_training, [2 1 3]);
else
    Tc_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Tc_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    Tc_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Tc_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    Tc_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Tc_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

delay_map = readNPY(fullfile(trainingCaseDir, 'delay_map.npy'));
delay_map_resized = readNPY(fullfile(trainingCaseDir, 'delay_map_resized.npy'));
delay_map_resized_training = readNPY(fullfile(trainingCaseDir, 'delay_map_resized_training.npy'));

if(in_CMR_View)
    delay_map_norm = permute(delay_map, [2 1 3]);
    delay_map_resized_norm = permute(delay_map_resized, [2 1 3]);
    delay_map_resized_training_norm = permute(delay_map_resized_training, [2 1 3]);
else
    delay_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(delay_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    delay_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(delay_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    delay_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(delay_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

Gd = readNPY(fullfile(trainingCaseDir, 'Gd.npy'));
Gd_resized = readNPY(fullfile(trainingCaseDir, 'Gd_resized.npy'));
Gd_resized_training = readNPY(fullfile(trainingCaseDir, 'Gd_resized_training.npy'));

if(in_CMR_View)
    Gd_norm = permute(Gd, [2 1 3 4]);
    Gd_resized_norm = permute(Gd_resized, [2 1 3 4]);
    Gd_resized_training_norm = permute(Gd_resized_training, [2 1 3 4]);
else
    Gd_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Gd, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    Gd_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Gd_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
    Gd_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Gd_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
end

h_fmap = figure; imagescn(fmap_resized_norm, [0 8], [1 SLC], 12); PerfColorMap;
figure; imagescn(Gd_resized_norm, [0 1.5], [1 SLC], [], 3);

saveas(h_fmap, fullfile(trainingCaseDir, 'fmap_resized_norm'), 'fig');

writeNPY(fmap_norm, fullfile(trainingCaseDir, 'fmap_norm.npy'));
writeNPY(fmap_resized_norm, fullfile(trainingCaseDir, 'fmap_resized_norm.npy'));
writeNPY(fmap_resized_training_norm, fullfile(trainingCaseDir, 'fmap_resized_training_norm.npy'));

writeNPY(PS_map_norm, fullfile(trainingCaseDir, 'PS_map_norm.npy'));
writeNPY(PS_map_resized_norm, fullfile(trainingCaseDir, 'PS_map_resized_norm.npy'));
writeNPY(PS_map_resized_training_norm, fullfile(trainingCaseDir, 'PS_map_resized_training_norm.npy'));

writeNPY(vb_map_norm, fullfile(trainingCaseDir, 'vb_map_norm.npy'));
writeNPY(vb_map_resized_norm, fullfile(trainingCaseDir, 'vb_map_resized_norm.npy'));
writeNPY(vb_map_resized_training_norm, fullfile(trainingCaseDir, 'vb_map_resized_training_norm.npy'));

writeNPY(visf_map_norm, fullfile(trainingCaseDir, 'visf_map_norm.npy'));
writeNPY(visf_map_resized_norm, fullfile(trainingCaseDir, 'visf_map_resized_norm.npy'));
writeNPY(visf_map_resized_training_norm, fullfile(trainingCaseDir, 'visf_map_resized_training_norm.npy'));

writeNPY(fSD_norm, fullfile(trainingCaseDir, 'fSD_norm.npy'));
writeNPY(fSD_resized_norm, fullfile(trainingCaseDir, 'fSD_resized_norm.npy'));
writeNPY(fSD_resized_training_norm, fullfile(trainingCaseDir, 'fSD_resized_training_norm.npy'));

writeNPY(Tc_map_norm, fullfile(trainingCaseDir, 'Tc_map_norm.npy'));
writeNPY(Tc_map_resized_norm, fullfile(trainingCaseDir, 'Tc_map_resized_norm.npy'));
writeNPY(Tc_map_resized_training_norm, fullfile(trainingCaseDir, 'Tc_map_resized_training_norm.npy'));

writeNPY(delay_map_norm, fullfile(trainingCaseDir, 'delay_map_norm.npy'));
writeNPY(delay_map_resized_norm, fullfile(trainingCaseDir, 'delay_map_resized_norm.npy'));
writeNPY(delay_map_resized_training_norm, fullfile(trainingCaseDir, 'delay_map_resized_training_norm.npy'));

writeNPY(Gd_norm, fullfile(trainingCaseDir, 'Gd_norm.npy'));
writeNPY(Gd_resized_norm, fullfile(trainingCaseDir, 'Gd_resized_norm.npy'));
writeNPY(Gd_resized_training_norm, fullfile(trainingCaseDir, 'Gd_resized_training_norm.npy'));
