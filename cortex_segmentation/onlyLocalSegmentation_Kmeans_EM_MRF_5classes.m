
function local_SegResults = onlyLocalSegmentation_Kmeans_EM_MRF_5classes(imagedata, header, ...
    brainmask, SplittedBrain, partsNumber, tissuelabels_seg, ...
    replicate, M, saveflag, trynumber)

% 5 classes: use all data, csf, gm, wm1, wm2, outlier
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
        parts_mask(find( SplittedBrain==i ) ) = 1;
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
    
    [csflabel, gmlabel, wmlabel1, wmlabel2, outlierlabel] = GetClassLabel_5classes(imagedata, kmeans_label);
    
    % create templates
%     csflabel = 5;
%     wmlabel1 = 3;
%     wmlabel2 = 4;
%     gmlabel = 2;
%     outlierlabel = 1;

    % create simulated atlas
    csfSeg = zeros(size(kmeans_label), 'uint8');
    csfSeg(find(kmeans_label==csflabel)) = 1;

    wmSeg1 = zeros(size(kmeans_label), 'uint8');
    wmSeg1(find(kmeans_label==wmlabel1)) = 1;

    wmSeg2 = zeros(size(kmeans_label), 'uint8');
    wmSeg2(find(kmeans_label==wmlabel2)) = 1;

    cortexSeg = zeros(size(kmeans_label), 'uint8');
    cortexSeg(find(kmeans_label==gmlabel)) = 1;

    outlierSeg = zeros(size(kmeans_label), 'uint8');
    outlierSeg(find(kmeans_label==outlierlabel)) = 1;
     
    sigmaCSF = 2*header.xvoxelsize;
    sigmaWM1 = 2*header.xvoxelsize;
    sigmaWM2 = 2*header.xvoxelsize;

    sigmaCortex = 2*header.xvoxelsize;
    sigmaOutlier = header.xvoxelsize;


    halfwidthCsf = 4;
    halfwidthWm1 = 4;
    halfwidthWm2 = 4;
    halfwidthCortex = 4;
    halfwidthOutlier = 2;

    [csfT, wmT1, wmT2, cortexT, outlierT] = CreateTemplates_5classes_Gaussian(header, csfSeg, wmSeg1, wmSeg2, cortexSeg, outlierSeg, ...
                        sigmaCSF, sigmaWM1, sigmaWM2, sigmaCortex, sigmaOutlier, ...
                        halfwidthCsf, halfwidthWm1, halfwidthWm2, halfwidthCortex, halfwidthOutlier);

    % run EM to refine results
    ncentres = 5;
    input_dim = 1;
    mix = gmm(input_dim, ncentres, 'full');

    clear wmSeg1 wmSeg2 gmSeg csfSeg outlierSeg kmeans_label
    % brainmask = mask;
    options = [];
    initType = 'mcd';
    initParameters.minimalP = 0.6;
    initParameters.eps = 0.02;
    [mix, x, indexes] = gmminit_5classes_image_uint8(mix, imagedata, header,...
                            csfT, cortexT, wmT1, wmT2, outlierT, parts_mask, ...
                            options, initType, initParameters);

    clear csfT cortexT wmT1 wmT2 outlierT
    
    samplefactor = 2;
    xsize = header.xsize;
    ysize = header.ysize;
    zsize = header.zsize;
    DIM = [ysize xsize zsize];
    ind1 = [1:samplefactor:DIM(1)]';    
    ind2 = [1:samplefactor:DIM(2)]';
    ind3 = [1:samplefactor:DIM(3)]';
    [ind1,ind2,ind3] = ndgrid(ind1,ind2,ind3);
    mix.sampleInd = ind1 + (ind2-1)*DIM(1) + (ind3-1)*DIM(1)*DIM(2);
    insideROI = find(mix.playing(mix.sampleInd(:)));
    mix.sampleInd = mix.sampleInd(insideROI);

    clear ind1 ind2 ind3 insideROI
    
    % set up options
    options = zeros(1,18);
    options(1) = 1;
    options(14) = 45;
    options(3) = 0.5;
    options(5) = 1;
    [mix, options, errlog] = gmmem_MRF_image_uint8_5classes(mix, x, options);

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
            
        if ( label == 4 ) % wm
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(4);
        end
    
        if ( label == 5 ) % outlier
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(5);
        end
    end
    
    if ( saveflag )
        filename = ['LocalSeg/' 'localSegResult_part' num2str(i) '.hdr'];
        SaveAnalyze(uint32(local_SegResults{i}), header, filename, 'Grey' );
    end

end

return