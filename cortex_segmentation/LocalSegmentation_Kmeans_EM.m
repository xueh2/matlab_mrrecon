
function local_SegResults = LocalSegmentation_Kmeans_EM(imagedata, header, ...
    brainmask, global_SegResult, SplittedBrain, partsNumber, tissuelabels_seg, ...
    replicate, M, saveflag)
% perform the local segmentation to refine the results

% allocate memory
local_SegResults = cell(partsNumber, 1);
for i = 1:partsNumber
    local_SegResults{i} = zeros(size(imagedata), 'uint8');
end

numberTissues = length(tissuelabels_seg);

% perform segmentation for each parts
for i = 1:partsNumber

    disp(['Segmenting ' num2str(i) ' part ...']);
      
    parts_mask = zeros(size(imagedata), 'uint8');
    % get image data as a feature vector
    for k = 1:numberTissues
        tissueLabel = tissuelabels_seg(k);
        parts_mask(find( (global_SegResult==tissueLabel) & (SplittedBrain==i) ) ) = 1;
    end
    % run kmeans to intially segmented the data
    [x, indexes] = kmean_init(imagedata, header, parts_mask);
    [IDX,C,sumd,D] = kmeans(x, numberTissues, 'start', M, 'distance', 'cityblock', 'display', 'iter');
    
    kmeans_label = zeros(size(imagedata), 'uint8');
    [ndata, ndim] = size(indexes);
    for tt = 1:ndata
        label = IDX(tt);
        kmeans_label(indexes(tt, 1), indexes(tt, 2), indexes(tt, 3)) = label;
    end
    
    if ( saveflag )
        SaveAnalyze(uint32(kmeans_label), header, 'LocalSeg/kmeansResult_parts.hdr', 'Grey' );
    end
    
    % create templates
    gmlabel = 1;
    wmlabel1 = 2;
    wmlabel2 = 3;

    wmSeg1 = zeros(size(LabeledSeg), 'uint8');
    wmSeg1(find(LabeledSeg==wmlabel1)) = 1;

    wmSeg2 = zeros(size(LabeledSeg), 'uint8');
    wmSeg2(find(LabeledSeg==wmlabel2)) = 1;

    gmSeg = zeros(size(LabeledSeg), 'uint8');
    gmSeg(find(LabeledSeg==gmlabel)) = 1;

    sigmaWM1 = 2*header.xvoxelsize;
    sigmaWM2 = 2*header.xvoxelsize;
    sigmaCortex = 2*header.xvoxelsize;

    halfwidthWm1 = 3;
    halfwidthWm2 = 2;
    halfwidthCortex = 3;

    [csfT, wmT1, wmT2, cortexT, outlierT] = CreateTemplates_5classes_Gaussian(header, csfSeg, wmSeg1, wmSeg2, gmSeg, outlierSeg, ...
                        sigmaCSF, sigmaWM1, sigmaWM2, sigmaCortex, sigmaOutlier, ...
                        halfwidthCsf, halfwidthWm1, halfwidthWm2, halfwidthCortex, halfwidthOutlier);

       
    % run EM to refine results

    % restore results
    local_SegResults{i} = kmeans_label;
end

return