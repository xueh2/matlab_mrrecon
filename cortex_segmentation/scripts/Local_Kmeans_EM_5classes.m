
% ====================================================== %
% load posterior
% [post_csf, header] = LoadAnalyze('posterior/post_csf.hdr', 'Grey' );
% [post_wm1, header] = LoadAnalyze('posterior/post_wm1.hdr', 'Grey' );
% [post_wm2, header] = LoadAnalyze('posterior/post_wm2.hdr', 'Grey' );
% [post_gm, header] = LoadAnalyze('posterior/post_gm.hdr', 'Grey' );
% [post_outlier, header] = LoadAnalyze('posterior/post_outlier.hdr', 'Grey' );

% ====================================================== %
% load previous global results
if ( MRF_flag_Global )
    filename = [Global_SegResult prefix '_segResult_5classes_MRF_uint8.hdr'];
else
    filename = [Global_SegResult prefix '_segResult_5classes.hdr'];
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

clear post_csf post_wm1 post_wm2 post_gm post_outlier

filename = [Local_SegResult prefix '_suspected_volume.hdr'];
SaveAnalyze(uint32(suspected_volume), header, filename, 'Grey' );

filename = [Local_SegResult prefix '_maxPosterior.hdr'];
SaveAnalyze(uint32(maxPosterior), header, filename, 'Grey' );

% x = maxPosterior(:);
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

tissuelabels_seg = [cortex wm1 wm2];
saveflag = 1;
% ======================================== %
k = length(tissuelabels_seg);
replicate = TryNumber;
if ( isempty(Kmeans_InitialCentres) == 0 )
    M = zeros(k,1);
    % initial intensity centroid for clustering
    M(:,1, 1) = [Kmeans_InitialCentres(2); Kmeans_InitialCentres(3); Kmeans_InitialCentres(4)]; 

    for m = 2:replicate
        values = rand(k, 1)-0.5;
        values = values .* 150;
        M(:,1, m) = M(:,1, 1) + values; 
    end
else
    M = [];
end
% ======================================== %

local_SegResults = LocalSegmentation_Kmeans_EM_3classes(imagedata, header, ...
    brainmask, global_SegResult, SplittedBrain, partsNumber, tissuelabels_seg, ...
    replicate, M, saveflag, trynumber, LabelPVs_Image);

clear SplittedBrain
% ====================================================== %
% 
% for i = 1:partsNumber
%     filename = ['LocalSeg/' 'local_SegResult_class_' num2str(i) '.hdr' ];
%     SaveAnalyze(uint32(local_SegResults{i}), header, filename, 'Grey' );
% end
After_local_step_5classes
% % ====================================================== %
% % combine local segmentation into the previous global results
% newGlobal_SegResult = Combine_Lodal_Global(global_SegResult, tissuelabels_seg, partsNumber, local_SegResults);
% % ====================================================== %
% filename = [Local_SegResult prefix '_newGlobal_SegResult.hdr'];
% SaveAnalyze(uint32(newGlobal_SegResult), header, filename,'Grey');
% 
% 
% % detect PVs
% csflabel = 1;
% wmlabel = [3 4];
% cortexlabel =  2;
% pvlabel = 0;
% 
% nonbrainlabel = 5;
% neighborwidth = 3;
% 
% LabeledSeg = LabelPVs_Seg_slow(double(newGlobal_SegResult), header, ...
%                                     csflabel, wmlabel, cortexlabel, pvlabel,...
%                                     nonbrainlabel, neighborwidth);
% filename = [Local_SegResult prefix '_newGlobal_SegResult_PVs.hdr'];
% SaveAnalyze(uint32(LabeledSeg), header, filename,'Grey');
% 
% volumeThreshold = 200;
% label = [2];
% 
% [label3D, largestComponent, LabeledSeg] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);
% 
% filename = [Local_SegResult prefix '_cortex_seg.hdr'];
% SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');
% 
% 
% volumeThreshold = 200;
% label = [3 4];
% 
% [label3D, largestComponent, LabeledSeg] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);
% 
% filename = [Local_SegResult prefix '_wm_seg.hdr'];
% SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');
