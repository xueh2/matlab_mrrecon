
function local_SegResults = LocalSegmentation_Kmeans_EM_MRF_2classes(imagedata, header, ...
    brainmask, global_SegResult, SplittedBrain, partsNumber, tissuelabels_seg, ...
    replicate, M, saveflag, trynumber)
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

    [gmlabel, wmlabel] = GetClassLabel_2classes(imagedata, kmeans_label);

    % create templates
%     gmlabel = 1;
%     wmlabel = 2;

    wmSeg = zeros(size(kmeans_label), 'uint8');
    wmSeg(find(kmeans_label==wmlabel)) = 1;

    gmSeg = zeros(size(kmeans_label), 'uint8');
    gmSeg(find(kmeans_label==gmlabel)) = 1;

    sigmaWM = 2*header.xvoxelsize;
    sigmaCortex = 3*header.xvoxelsize;

    halfwidthWm = 3;
    halfwidthCortex = 3;

    [wmT, cortexT] = CreateTemplates_2classes_Gaussian(header, wmSeg, gmSeg, ...
                                sigmaWM, sigmaCortex, halfwidthWm, halfwidthCortex);
   
    % run EM to refine results
    ncentres = 2;
    input_dim = 1;
    mix = gmm(input_dim, ncentres, 'full');

    clear wmSeg gmSeg kmeans_label
    % brainmask = mask;
    options = [];
    initType = 'mcd';
    initParameters.minimalP = 0.6;
    initParameters.eps = 0.02;
    [mix, x, indexes] = gmminit_2classes_image_uint8(mix, imagedata, header,...
                                cortexT, wmT, parts_mask, ...
                                options, initType, initParameters);

    clear cortexT wmT
    
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
    options(14) = 35;
    options(3) = 1e-3;
    options(5) = 1;
    [mix, options, errlog] = gmmem_MRF_image_uint8_2classes(mix, x, options);

    % compute the posterior probablities
    [post, a] = gmmpost_image(mix, x);

    % restore results
    [ndata, ndim] = size(mix.indexes);
    for tt = 1:ndata
        label = find(post(tt,:) == max(post(tt,:)));
        if ( label == 1 ) % cortex
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(1);
        end
        
        if ( label == 2 ) % wm
            local_SegResults{i}(mix.indexes(tt, 1), mix.indexes(tt, 2), mix.indexes(tt, 3)) = tissuelabels_seg(2);
        end
        
    end
    
    if ( saveflag )
        filename = ['LocalSeg/' 'localSegResult_part' num2str(i) '.hdr'];
        SaveAnalyze(uint32(local_SegResults{i}), header, filename, 'Grey' );
    end

end

return