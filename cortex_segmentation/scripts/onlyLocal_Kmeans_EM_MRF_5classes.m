
% ====================================================== %
% split the brain
filename = [Local_SegResult prefix '_SplittedBrain.hdr'];
[SplittedBrain, header] = LoadAnalyze(filename, 'Grey');
% ====================================================== %

% ====================================================== %
% local segmentation

csf = 1;
wm1 = 3;
wm2 = 4;
cortex =  2;
pv = 0;
nonbrain = 5;

tissuelabels_seg = [csf cortex wm1 wm2 nonbrain];

saveflag = 1;
% ======================================== %
k = length(tissuelabels_seg);
replicate = TryNumber;
if ( isempty(Kmeans_InitialCentres) == 0 )
    M = zeros(k,1);
    % initial intensity centroid for clustering
    M(:,1, 1) = Kmeans_InitialCentres; 

    for m = 2:replicate
        values = rand(k, 1)-0.5;
        values = values .* 150;
        M(:,1, m) = M(:,1, 1) + values; 
    end
else
    M = [];
end
% ======================================== %

local_SegResults = onlyLocalSegmentation_Kmeans_EM_MRF_5classes(imagedata, header, ...
    brainmask, SplittedBrain, partsNumber, tissuelabels_seg, replicate, M, saveflag, trynumber);

clear SplittedBrain
% ====================================================== %
% 
% for i = 1:partsNumber
%     filename = ['LocalSeg/' 'local_SegResult_class_' num2str(i) '.hdr' ];
%     SaveAnalyze(uint32(local_SegResults{i}), header, filename, 'Grey' );
% end

% ====================================================== %
% combine local segmentation into the previous global results
global_SegResult = zeros(size(SplittedBrain), 'uint8');
newGlobal_SegResult = Combine_Lodal_Global(global_SegResult, tissuelabels_seg, partsNumber, local_SegResults);
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

LabeledSeg = LabelPVs_Seg(double(newGlobal_SegResult), header, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel, neighborwidth);
filename = [Local_SegResult prefix '_newGlobal_SegResult_PVs.hdr'];
SaveAnalyze(uint32(LabeledSeg), header, filename,'Grey');

volumeThreshold = 200;
label = [2];

[label3D, largestComponent] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);

filename = [Local_SegResult prefix '_cortex_seg.hdr'];
SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');


volumeThreshold = 200;
label = [3 4];

[label3D, largestComponent] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);

filename = [Local_SegResult prefix '_wm_seg.hdr'];
SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');
