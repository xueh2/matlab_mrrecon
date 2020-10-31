
% ====================================================== %
% load posterior
% [post_csf, header] = LoadAnalyze('posterior/post_csf.hdr', 'Grey' );
% [post_wm1, header] = LoadAnalyze('posterior/post_wm1.hdr', 'Grey' );
% [post_wm2, header] = LoadAnalyze('posterior/post_wm2.hdr', 'Grey' );
% [post_gm, header] = LoadAnalyze('posterior/post_gm.hdr', 'Grey' );
% [post_outlier, header] = LoadAnalyze('posterior/post_outlier.hdr', 'Grey' );
% 
% ====================================================== %
% load previous global results

if ( MRF_flag_Global )
    filename = [Global_SegResult GlobalPrefix '_segResult_5classes_MRF_uint8.hdr'];
else
    filename = [Global_SegResult GlobalPrefix '_segResult_5classes.hdr'];
end
[global_SegResult, header] = LoadAnalyze(filename, 'Grey' );

csf = 1;
wm1 = 3;
wm2 = 4;
cortex =  2;
pv = 0;
nonbrain = 5;

% ====================================================== %

% create current wm & gm mask
cortex_mask = GetEachRegion(global_SegResult, cortex);
wm_mask = GetEachRegion(global_SegResult, [wm1 wm2]);
% ====================================================== %

% compute suspected voxel points
Low_Thres_Uint = Low_Thres*2048;

[suspected_volume, maxPosterior]  = EstimateSuspectedVoxels(brainmask, post_csf, post_wm1+post_wm2, post_gm, post_outlier, Low_Thres_Uint);

% Low_Thres_Uint = 1;
% [suspected_volume, entropy] = EstimateSuspectedVoxels_Entropy(brainmask, post_csf, post_wm1+post_wm2, post_gm, post_outlier, Low_Thres_Uint);

% clear post_csf post_wm1 post_wm2 post_gm post_outlier

filename = [Local_SegResult prefix '_suspected_volume.hdr'];
SaveAnalyze(uint32(suspected_volume), header, filename, 'Grey' );

filename = [Local_SegResult prefix '_maxPosterior.hdr'];
SaveAnalyze(uint32(maxPosterior), header, filename, 'Grey' );

% filename = [Local_SegResult prefix '_entropy.hdr'];
% SaveAnalyze(uint32(2048*entropy), header, filename, 'Grey' );

% x = maxPosterior(:);
% pp = x(find(x>0));
% hist(double(pp(:)), 2048);

% x = entropy(:);
% pp = x(find(x>0));
% hist(double(pp(:)), 2048);

clear maxPosterior

suspected_volume = suspected_volume & (cortex_mask | wm_mask);

clear cortex_mask wm_mask
% remove single points (outlier ... ?)
suspected_volume = RemoveSinglePoints(suspected_volume, volumeThres);
filename = [Local_SegResult prefix '_suspected_volume3.hdr'];
SaveAnalyze(uint32(suspected_volume), header, filename, 'Grey' );

% ====================================================== %

% ====================================================== %
% clustering the suspected points

% centres: x, y, z (col, row, depth)
centres = ClusteringSuspectedPoints(suspected_volume, partsNumber, TryNumber, initialCentres );
clear suspected_volume
% ====================================================== %

% ====================================================== %
% split the brain
SplittedBrain = SplitBrain_Voronoi(brainmask, centres);

filename = [Local_SegResult prefix '_SplittedBrain.hdr'];
SaveAnalyze(uint32(SplittedBrain), header, filename, 'Grey' );
% ====================================================== %

% ====================================================== %
% local segmentation

csf = 1;
wm1 = 3;
wm2 = 4;
cortex =  2;
pv = 0;
nonbrain = 5;

% tissuelabels_seg = [csf cortex wm1 wm2 nonbrain];
tissuelabels_seg = [cortex wm1 wm2];

saveflag = 1;
% ======================================== %
if ( exist('svmparameters') == 0 )
    svmparameters = struct('maxN', 2000, 'foldnum', 5, 'C_low', -5, 'C_high', 15, 'C_step', 2, ...
        'gamma_low', -15, 'gamma_high', 3, 'gamma_step', 2, 'fine_range', 2, 'fine_step', 0.25, 'finesearch_flag', 1, 'amplifier', 5);
end
local_SegResults = LocalSegmentation_Kmeans_EM_3classes_SVM(imagedata, header, ...
    brainmask, global_SegResult, SplittedBrain, partsNumber, saveflag, tissuelabels_seg, ...
    post_gm, post_wm1, post_wm2, svmparameters, LabelPVs_Image);

clear SplittedBrain
clear post_csf post_wm1 post_wm2 post_gm post_outlier
% ====================================================== %
% combine local segmentation into the previous global results
newGlobal_SegResult = Combine_Lodal_Global2(global_SegResult, tissuelabels_seg, ...
    partsNumber, local_SegResults, LabelPVs_Image);
% ====================================================== %
filename = [Local_SegResult prefix '_newGlobal_SegResult.hdr'];
SaveAnalyze(uint32(newGlobal_SegResult), header, filename,'Grey');


% detect PVs
csflabel = 1;
wmlabel = [3 4];
cortexlabel =  2;
pvlabel = 0;

nonbrainlabel = 5;
neighborwidth = 3;

LabeledSeg = LabelPVs_Seg_slow(double(newGlobal_SegResult), header, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel, neighborwidth);
filename = [Local_SegResult prefix '_newGlobal_SegResult_PVs.hdr'];
SaveAnalyze(uint32(LabeledSeg), header, filename,'Grey');

volumeThreshold = 200;
label = [2];

[label3D, largestComponent, LabeledSeg] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);

filename = [Local_SegResult prefix '_cortex_seg.hdr'];
SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');


volumeThreshold = 200;
label = [3 4];

[label3D, largestComponent, LabeledSeg] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);

filename = [Local_SegResult prefix '_wm_seg.hdr'];
SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');
