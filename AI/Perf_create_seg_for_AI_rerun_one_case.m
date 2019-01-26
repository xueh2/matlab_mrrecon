[configName, scannerID, patientID, studyID, measurementID, study_date, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(caseName);

caseDir = fullfile(resDir, study_date, caseName);
debugDir = fullfile(resDir, study_date, caseName, 'DebugOutput');
labelDir = fullfile(resDir, study_date, caseName, 'Label');
trainingDataDir = fullfile(resDir, study_date, caseName, 'Training_Data');
contourDir = fullfile(resDir, study_date, caseName, 'Contour_result');
trainingCaseDir = fullfile(trainingDir, caseName)

fmap = readNPY(fullfile(trainingDataDir, 'fmap.npy'));
RO = size(fmap, 1);
E1 = size(fmap, 2);

fmap_resized = readNPY(fullfile(trainingDataDir, 'fmap_resized.npy'));
RO_resized = size(fmap_resized, 1);
E1_resized = size(fmap_resized, 2);

fmap_resized_training = readNPY(fullfile(trainingDataDir, 'fmap_resized_training.npy'));
Gd_resized_training = readNPY(fullfile(trainingDataDir, 'Gd_resized_training.npy'));
fSD_resized_training = readNPY(fullfile(trainingDataDir, 'fSD_resized_training.npy'));

contour_roi = readNPY(fullfile(trainingDataDir, 'roi.npy'));

NN_res = load(fullfile(contourDir, 'NN_res'));

RO_resized_training = size(fmap_resized_training, 1);
E1_resized_training = size(fmap_resized_training, 2);

cd(contourDir)

SLC = size(Gd_resized_training, 4);
h = figure; imagescn(Gd_resized_training, [0 1.5], [1 SLC], scale_factor, 3);
h2 = figure; imagescn(fmap_resized_training, [0 8], [1 SLC], scale_factor); PerfColorMap;
h3 = figure; imagescn(fSD_resized_training, [0 1], [1 SLC], scale_factor); PerfColorMap;

pos = get(h, 'Position');
set(h,'Position',[5 8 1.5*pos(3) 1.5*pos(4)])
set(h2,'Position',[3 -1 2*pos(3) 2*pos(4)])
set(h3,'Position',[3+pos(3) 0 pos(3) pos(4)])

disp(['please edit contours : ']);    

pause
closeall

roi_name = fullfile(contourDir, 'roi_new.mat');
tic
copyfile(fullfile(trainingDataDir, '*.npy'), trainingCaseDir);
copyfile(fullfile(contourDir, 'roi*.mat'), trainingCaseDir);
toc

% generate mask
roi = load(roi_name);

% roi = [ps_x pe_x; ps_y pe_y; aif_s aif_e ];

Seg = [];
for slc=1:SLC        
    startROI = 1 + 4*(slc-1);        
    S = perf_generate_seg_from_manual_roi(roi, contour_roi, fmap(:,:,slc), fmap_resized(:,:,slc), fmap_resized_training(:,:,slc), startROI, slc, 1);       
    Seg = [Seg S];
end

save(fullfile(contourDir, 'Seg.mat'), 'Seg');
copyfile(fullfile(contourDir, 'Seg.mat'), trainingCaseDir);   

[h_Gd, h_fmap, h_fSD] = perf_plot_NN_after_editing(trainingCaseDir, Seg);

saveas(h_Gd, fullfile(trainingCaseDir, 'seg'), 'fig');
saveas(h_fmap, fullfile(trainingCaseDir, 'seg_fmap'), 'fig');
saveas(h_fSD, fullfile(trainingCaseDir, 'seg_fSD'), 'fig');    

pos = get(h_fSD, 'Position');
set(h_Gd,'Position',[5 8 1.5*pos(3) 1.5*pos(4)])
set(h_fmap,'Position',[3 -1 2*pos(3) 2*pos(4)])
set(h_fSD,'Position',[3+pos(3) 0 pos(3) pos(4)])

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

fmap_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fmap, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
fmap_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fmap_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
fmap_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fmap_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

fSD = readNPY(fullfile(trainingCaseDir, 'fSD.npy'));
fSD_resized = readNPY(fullfile(trainingCaseDir, 'fSD_resized.npy'));
fSD_resized_training = readNPY(fullfile(trainingCaseDir, 'fSD_resized_training.npy'));
fSD_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fSD, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
fSD_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fSD_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
fSD_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(fSD_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

PS_map = readNPY(fullfile(trainingCaseDir, 'PS_map.npy'));
PS_map_resized = readNPY(fullfile(trainingCaseDir, 'PS_map_resized.npy'));
PS_map_resized_training = readNPY(fullfile(trainingCaseDir, 'PS_map_resized_training.npy'));
PS_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(PS_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
PS_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(PS_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
PS_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(PS_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

vb_map = readNPY(fullfile(trainingCaseDir, 'vb_map.npy'));
vb_map_resized = readNPY(fullfile(trainingCaseDir, 'vb_map_resized.npy'));
vb_map_resized_training = readNPY(fullfile(trainingCaseDir, 'vb_map_resized_training.npy'));
vb_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(vb_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
vb_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(vb_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
vb_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(vb_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

visf_map = readNPY(fullfile(trainingCaseDir, 'visf_map.npy'));
visf_map_resized = readNPY(fullfile(trainingCaseDir, 'visf_map_resized.npy'));
visf_map_resized_training = readNPY(fullfile(trainingCaseDir, 'visf_map_resized_training.npy'));
visf_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(visf_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
visf_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(visf_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
visf_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(visf_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

Tc_map = readNPY(fullfile(trainingCaseDir, 'Tc_map.npy'));
Tc_map_resized = readNPY(fullfile(trainingCaseDir, 'Tc_map_resized.npy'));
Tc_map_resized_training = readNPY(fullfile(trainingCaseDir, 'Tc_map_resized_training.npy'));
Tc_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Tc_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
Tc_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Tc_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
Tc_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Tc_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

delay_map = readNPY(fullfile(trainingCaseDir, 'delay_map.npy'));
delay_map_resized = readNPY(fullfile(trainingCaseDir, 'delay_map_resized.npy'));
delay_map_resized_training = readNPY(fullfile(trainingCaseDir, 'delay_map_resized_training.npy'));
delay_map_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(delay_map, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
delay_map_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(delay_map_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
delay_map_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(delay_map_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

Gd = readNPY(fullfile(trainingCaseDir, 'Gd.npy'));
Gd_resized = readNPY(fullfile(trainingCaseDir, 'Gd_resized.npy'));
Gd_resized_training = readNPY(fullfile(trainingCaseDir, 'Gd_resized_training.npy'));
Gd_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Gd, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
Gd_resized_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Gd_resized, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    
Gd_resized_training_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Gd_resized_training, 1), 2), aif_scan_geometry_info.aif_slice_dir, aif_scan_geometry_info.aif_read_dir, aif_scan_geometry_info.aif_phase_dir);    

h_fmap = figure; imagescn(fmap_resized_norm, [0 8], [1 SLC], 12); PerfColorMap;
%     figure; imagescn(fSD_resized_norm, [0 2], [1 SLC]); PerfColorMap;
%     figure; imagescn(PS_map_resized_norm, [0 8], [1 SLC]); PerfColorMap;
%     figure; imagescn(Tc_map_resized_norm, [0 10], [1 SLC]);
%     figure; imagescn(delay_map_resized_norm, [0 2], [1 SLC]); PerfColorMap;
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

Seg = load(fullfile(trainingCaseDir, 'Seg.mat'));

Seg = Seg.Seg;

slc_dir = aif_scan_geometry_info.aif_slice_dir;
read_dir = aif_scan_geometry_info.aif_read_dir;
phase_dir = aif_scan_geometry_info.aif_phase_dir;

aif_cin_foot_peak_valley = analyze75read(fullfile(debugDir, 'aif_cin_foot_peak_valley.hdr'))
cin = analyze75read(fullfile(debugDir, 'Input_cin_all_computeFlowMap.hdr'));
aif_mask = analyze75read(fullfile(debugDir, 'aif_LV_mask_final.hdr'));
aif_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(aif_mask, 1), 2), slc_dir, read_dir, phase_dir);

N = size(Gd, 3);

[aif_X, aif_Y] = find(aif_mask_norm>0);
cx = mean(aif_X);
cy = mean(aif_Y);

RO = size(Gd_resized_norm, 1);
E1 = size(Gd_resized_norm, 2);

rx = size(Gd_resized_norm,1)/size(aif_mask_norm,1);
ry = size(Gd_resized_norm,2)/size(aif_mask_norm,2);

cx = cx*rx;
cy = cy*ry;

cx = round(cx);
cy = round(cy);

ps_x = cx-dstX*rS;
ps_y = cy-dstY*rS;

if(ps_x<1)
    ps_x = 1;
end
if(ps_y<1)
    ps_y = 1;
end

pe_x = ps_x+dstX*2*rS;
pe_y = ps_y+dstY*2*rS;

if(pe_x>rS*RO)
    d = pe_x - rS*RO;
    pe_x = rS*RO;
    ps_x = pe_x - 2*dstX*rS; 
end

if(pe_y>rS*E1)
    d = pe_y - rS*E1;
    pe_y = rS*E1;
    ps_y = pe_y - 2*dstY*rS; 
end 

aif_s = aif_cin_foot_peak_valley(1)-8;
if(aif_s<1)
    aif_s = 1;
end
aif_e = aif_s+48;

if(aif_e>N)
    aif_e = N;
    aif_s = aif_e - 48;
end

endo_epi_rv_rvi_resized_mask_norm = [];
for slc=1:SLC
    Seg(slc).endo_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).epi_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).epi_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).rv_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).rv_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).myo_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).myo_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).rvi_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).rvi_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rvi_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rvi_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rv_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rv_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rv_rvi_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rv_rvi_mask, 1), 2), slc_dir, read_dir, phase_dir);
    [I, J] = find(Seg(slc).rvi_mask_norm>0);
    Seg(slc).rvi_norm = [mean(J) mean(I)];

    Seg(slc).endo_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).epi_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).epi_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).rv_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).rv_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).myo_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).myo_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).rvi_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).rvi_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rvi_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rvi_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rv_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rv_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rv_rvi_resized_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rv_rvi_resized_mask, 1), 2), slc_dir, read_dir, phase_dir);
    [I, J] = find(Seg(slc).rvi_resized_mask_norm>0);
    Seg(slc).rvi_resized_norm = [mean(J) mean(I)];

    if(slc==1)
        endo_epi_rv_rvi_resized_mask_norm = Seg(slc).endo_epi_rv_rvi_resized_mask_norm;
    else
        endo_epi_rv_rvi_resized_mask_norm = cat(3, endo_epi_rv_rvi_resized_mask_norm, Seg(slc).endo_epi_rv_rvi_resized_mask_norm);
    end

    Seg(slc).endo_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).epi_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).epi_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).rv_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).rv_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).myo_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).myo_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).rvi_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).rvi_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rvi_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rvi_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rv_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rv_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    Seg(slc).endo_epi_rv_rvi_resized_training_mask_norm = Matlab_gt_CMR_normal_orientation(flipdim(flipdim(Seg(slc).endo_epi_rv_rvi_resized_training_mask, 1), 2), slc_dir, read_dir, phase_dir);
    [I, J] = find(Seg(slc).rvi_resized_training_mask_norm>0);
    Seg(slc).rvi_resized_training_norm = [mean(J) mean(I)];

    if(numel(Seg(slc).roi)==4)
        Seg(slc).roi(5) = aif_s;
        Seg(slc).roi(6) = aif_e;
    end
    
    Seg(slc).roi_norm = [ps_x pe_x ps_y pe_y aif_s aif_e] - 1;
end

h_masks = figure; imagescn(endo_epi_rv_rvi_resized_mask_norm, [], [1 SLC], [12]);  
hold on
plot(Seg(slc).roi_norm(3), Seg(slc).roi_norm(1), 'r+');
plot(Seg(slc).roi_norm(4), Seg(slc).roi_norm(1), 'yo');
plot(Seg(slc).roi_norm(3), Seg(slc).roi_norm(2), 'yo');
plot(Seg(slc).roi_norm(4), Seg(slc).roi_norm(2), 'r+');

save(fullfile(trainingCaseDir, 'Seg_norm'), 'Seg');

saveas(h_masks, fullfile(trainingCaseDir, 'endo_epi_rv_rvi_resized_mask_norm'), 'fig');

tic
command = [' scp ' fullfile(trainingCaseDir, '*_norm.npy') ' xueh2@168.61.47.145:' fullfile('~/mrprogs/gadgetron_CMR_ML-source/TrainingData/Perf_SAX_SEG', siteName, caseName)];
disp(command)
dos(command, '-echo');
disp(['copy npy norm - ' num2str(toc)]);

tic
command = [' scp ' fullfile(trainingCaseDir, 'Seg_norm.mat') ' xueh2@168.61.47.145:' fullfile('~/mrprogs/gadgetron_CMR_ML-source/TrainingData/Perf_SAX_SEG', siteName, caseName)];
disp(command)
dos(command, '-echo');
disp(['copy Seg_norm - ' num2str(toc)]);