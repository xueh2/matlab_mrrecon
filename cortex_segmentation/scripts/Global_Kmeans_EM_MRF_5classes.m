
% [x, indexes] = kmean_init(imagedata, header, brainmask);
% 
% k = classnumber;
% if ( isempty(Kmeans_InitialCentres) == 0)
% %     k = classnumber;
%     M = zeros(classnumber, 1);
%     M(:,1, 1) = Kmeans_InitialCentres; % initial intensity centroid for clustering
% 
%     for m = 2:TryNumber
%         values = rand(k, 1)-0.5;
%         values = values .* 150;
%         M(:,1, m) = M(:,1, 1) + values; 
%     end
%     [IDX,C,sumd,D] = kmeans(x, k, 'start', M, 'distance', 'cityblock', 'display', 'iter', 'EmptyAction', 'drop');
% else
%     [IDX,C,sumd,D] = kmeans(x, k, 'distance', 'cityblock', 'display', 'iter', 'Replicates', TryNumber, 'EmptyAction', 'drop');
% end
% % save results
% kmeans_label = zeros(size(imagedata), 'uint32');
% [ndata, ndim] = size(indexes);
% for i = 1:ndata
%     label = IDX(i);
%     kmeans_label(indexes(i, 1), indexes(i, 2), indexes(i, 3)) = label;
% end
% clear IDX C sumd D
% filename = [Global_SegResult prefix '_kmeansResult.hdr'];
% SaveAnalyze(kmeans_label, header, filename, 'Grey' );

Perform_Kmeans

[csflabel, cortexlabel, wmlabel1, wmlabel2, outlierlabel] = GetClassLabel_5classes(imagedata, kmeans_label);

wmlabel = [wmlabel1 wmlabel2];
pvlabel = 0;

neighborwidth = 3;

% LabeledSeg = LabelPVs_Seg(double(kmeans_label), header, ...
%                         csflabel, wmlabel, cortexlabel, pvlabel,...
%                         nonbrainlabel, neighborwidth);
%                     
% filename = [Global_SegResult prefix '_kmeansResult_PVs.hdr'];
% SaveAnalyze(uint32(LabeledSeg), header, filename, 'Grey' );.
LabeledSeg = kmeans_label;
% create simulated atlas
clear kmeans_label
% csflabel = 5;
% wmlabel1 = 3;
% wmlabel2 = 4;
% gmlabel = 1;
% outlierlabel = 2;

csfSeg = zeros(size(LabeledSeg), 'uint8');
csfSeg(find(LabeledSeg==csflabel)) = 1;

wmSeg1 = zeros(size(LabeledSeg), 'uint8');
wmSeg1(find(LabeledSeg==wmlabel1)) = 1;

wmSeg2 = zeros(size(LabeledSeg), 'uint8');
wmSeg2(find(LabeledSeg==wmlabel2)) = 1;

gmSeg = zeros(size(LabeledSeg), 'uint8');
gmSeg(find(LabeledSeg==cortexlabel)) = 1;

outlierSeg = zeros(size(LabeledSeg), 'uint8');
outlierSeg(find(LabeledSeg==outlierlabel)) = 1;

sigmaCSF = 2*header.xvoxelsize;
sigmaWM1 = 2*header.xvoxelsize;
sigmaWM2 = 2*header.xvoxelsize;

sigmaCortex = 2*header.xvoxelsize;
sigmaOutlier = header.xvoxelsize;


halfwidthCsf = 3;
halfwidthWm1 = 3;
halfwidthWm2 = 2;
halfwidthCortex = 3;
halfwidthOutlier = 2;

clear LabeledSeg

[csfT, wmT1, wmT2, cortexT, outlierT] = CreateTemplates_5classes_Gaussian(header, csfSeg, wmSeg1, wmSeg2, gmSeg, outlierSeg, ...
                    sigmaCSF, sigmaWM1, sigmaWM2, sigmaCortex, sigmaOutlier, ...
                    halfwidthCsf, halfwidthWm1, halfwidthWm2, halfwidthCortex, halfwidthOutlier);


% EM segmentation

% set up mixture model
ncentres = 5;
input_dim = 1;
mix = gmm(input_dim, ncentres, 'full');

clear csfSeg wmSeg1 wmSeg2 gmSeg outlierSeg
% brainmask = mask;
options = [];
initType = 'mcd';
initParameters.minimalP = 0.6;
initParameters.eps = 0.02;
[mix, x, indexes] = gmminit_5classes_image_uint8(mix, imagedata, header,...
                            csfT, cortexT, wmT1, wmT2, outlierT, brainmask, ...
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
options(14) = 35;
options(3) = 1.0e-3;
options(5) = 1;
[mix, options, errlog] = gmmem_MRF_image_uint8_5classes(mix, x, options);

% compute the posterior probablities
[post, a] = gmmpost_image(mix, x);

% save results
% results = zeros(size(imagedata), 'uint32');
% [ndata, ndim] = size(mix.indexes);
% for i = 1:ndata
%     label = find(post(i,:) == max(post(i,:)));
%     results(mix.indexes(i, 1), mix.indexes(i, 2), mix.indexes(i, 3)) = label(1);
% end

filename = [Global_SegResult prefix '_segResult_5classes_MRF_uint8.hdr'];

% SaveAnalyze(results, header, filename, 'Grey' );

% save posterior as data file
% [post_csf, post_wm1, post_wm2, post_gm, post_outlier] = GetPostImage_5classes(mix, header, post);
% 
% SaveAnalyze(uint32(post_csf), header, 'posterior/post_csf.hdr', 'Grey' );
% SaveAnalyze(uint32(post_wm1), header, 'posterior/post_wm1.hdr', 'Grey' );
% SaveAnalyze(uint32(post_wm2), header, 'posterior/post_wm2.hdr', 'Grey' );
% SaveAnalyze(uint32(post_gm), header, 'posterior/post_gm.hdr', 'Grey' );
% SaveAnalyze(uint32(post_outlier), header, 'posterior/post_outlier.hdr', 'Grey' );
After_gmmem_step_5classes
% ====================================================== %
% detect PVs
csflabel = 1;
wmlabel = [3 4];
cortexlabel =  2;
pvlabel = 0;

nonbrainlabel = 5;
neighborwidth = 3;

LabeledSeg = LabelPVs_Seg(mix, double(results), header, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel, neighborwidth);

filename = [Global_SegResult prefix '_segResult_5classes_MRF_uint8_PVs.hdr'];
SaveAnalyze(uint32(LabeledSeg), header, filename, 'Grey');


volumeThreshold = 200;
label = [2];

[label3D, largestComponent, LabeledSeg] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);

filename = [Global_SegResult prefix '_cortex_seg.hdr'];
SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');


volumeThreshold = 200;
label = [3 4];

[label3D, largestComponent, LabeledSeg] = RegionVolumeFilter_cortex(LabeledSeg, header, volumeThreshold, label);

filename = [Global_SegResult prefix '_wm_seg.hdr'];
SaveAnalyze(uint32(largestComponent), header, filename, 'Grey');


label3D = zeros(size(LabeledSeg), 'uint32');
label3D(LabeledSeg == csflabel) = 1;
filename = [Global_SegResult prefix '_csf_seg.hdr'];
SaveAnalyze(label3D, header, filename, 'Grey');

clear LabeledSeg post label3D largestComponent