
function local_SegResults = LocalSegmentation_Kmeans_EM_MRF_4classes(imagedata, header, ...
    brainmask, global_SegResult, SplittedBrain, partsNumber, tissuelabels_seg, ...
    replicate, M, saveflag, trynumber)

% 4 classes: use all data, csf, gm, wm, outlier
% perform the local segmentation to refine the results

% allocate memory
local_SegResults = cell(partsNumber, 1);
for i = 1:partsNumber
    local_SegResults{i} = zeros(size(imagedata), 'uint8');
end

numberTissues = length(tissuelabels_seg);
% parts_mask = zeros(size(imagedata), 'uint8');

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
    if ( isempty(M) == 0 )
        [IDX,C,sumd,D] = kmeans(x, numberTissues, 'start', M, 'distance', 'cityblock', 'display', 'iter', 'EmptyAction', 'drop');
    else
        [IDX,C,sumd,D] = kmeans(x, numberTissues, 'replicates', trynumber, 'distance', 'cityblock', 'display', 'iter', 'EmptyAction', 'drop');
    end
    
    kmeans_label = zeros(size(imagedata), 'uint8');
    [ndata, ndim] = size(indexes);
    for tt = 1:ndata
        label = IDX(tt);
        kmeans_label(indexes(tt, 1), indexes(tt, 2), indexes(tt, 3)) = label;
    end
    
    if ( saveflag )
        SaveAnalyze(uint32(kmeans_label), header, 'LocalSeg/kmeansResult_parts.hdr', 'Grey' );
    end
    
    [csflabel, gmlabel, wmlabel, outlierlabel] = GetClassLabel_4classes(imagedata, kmeans_label);
    
    % create templates
%     gmlabel = 1;
%     wmlabel = 2;
%     csflabel = 4;
%     outlierlabel = 3;

    wmSeg = zeros(size(kmeans_label), 'uint8');
    wmSeg(find(kmeans_label==wmlabel)) = 1;

    gmSeg = zeros(size(kmeans_label), 'uint8');
    gmSeg(find(kmeans_label==gmlabel)) = 1;

    csfSeg = zeros(size(kmeans_label), 'uint8');
    csfSeg(find(kmeans_label==csflabel)) = 1;

    outlierSeg = zeros(size(kmeans_label), 'uint8');
    outlierSeg(find(kmeans_label==outlierlabel)) = 1;
     
    sigmaCSF = header.xvoxelsize;
    sigmaWM = 2*header.xvoxelsize;
    sigmaCortex = 2*header.xvoxelsize;
    sigmaOutlier = header.xvoxelsize;


    halfwidthCsf = 2;
    halfwidthWm = 3;
    halfwidthCortex = 3;
    halfwidthOutlier = 2;

    [csfT, wmT, cortexT, outlierT] = CreateTemplates_4classes_Gaussian(header, csfSeg, wmSeg, gmSeg, outlierSeg, ...
                        sigmaCSF, sigmaWM, sigmaCortex, sigmaOutlier, halfwidthCsf, halfwidthWm, halfwidthCortex, halfwidthOutlier);
                  
    % run EM to refine results
    ncentres = 4;
    input_dim = 1;
    mix = gmm(input_dim, ncentres, 'full');

    clear wmSeg gmSeg csfSeg outlierSeg kmeans_label
    % brainmask = mask;
    options = [];
    initType = 'mcd';
    initParameters.minimalP = 0.6;
    initParameters.eps = 0.02;
    [mix, x, indexes] = gmminit_4classes_image(mix, imagedata, header,...
                            csfT, cortexT, wmT, outlierT, parts_mask, ...
                            options, initType, initParameters);

    clear csfT cortexT wmT outlierT

    % set up options
    options = zeros(1,18);
    options(1) = 1;
    options(14) = 100;
    options(3) = 1e-3;
    options(5) = 1;
    [mix, options, errlog] = gmmem_image(mix, x, options);

    % compute the posterior probablities
    [post, a] = gmmpost_image(mix, x);

    % restore results
    [ndata, ndim] = size(mix.indexes);
    for tt = 1:ndata
        label = find(post(tt,:) == max(post(tt,:)));
        if ( label == 1 ) % csf
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(1);
        end
        
        if ( label == 2 ) % cortex
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(2);
        end
        
        if ( label == 3 ) % wm
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(3);
        end
        
        if ( label == 4 ) % outlier
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(4);
        end
    end
    
    if ( saveflag )
        filename = ['LocalSeg/' 'localSegResult_part' num2str(i) '.hdr'];
        SaveAnalyze(uint32(local_SegResults{i}), header, filename, 'Grey' );
    end

end

return