
function [LabeledSeg, mix] = LabelPVs_Seg_gmm_5classes_atlas_MEX(mix, segResult, ...
                                    csflabel, wmlabel, cortexlabel, pvlabel,...
                                    nonbrainlabel, partsmask, globalSeg)
% according to the prior knowledge to label the pvs from the initial
% segmentation
% the pvs arel labelled as pvlabel

disp('Label partial volumes ...');

neighborwidth = mix.neighborwidth;
lamda = mix.lamda;

num = length(wmlabel);
if (num>1)
    
    for i = 2:num
        segResult(find(segResult==wmlabel(i))) = wmlabel(1);
    end
    
end

neighborNum = mix.neighborNum;
lamda = mix.lamda;
xsize = mix.header.xsize;
ysize = mix.header.ysize;
zsize = mix.header.zsize;

[LabeledSeg, Priors, Atlas] = LabelPVs_Seg_gmm_atlas_Local_MEX(single(mix.priors), mix.indexes, segResult, neighborNum, lamda,...
    xsize, ysize, zsize, csflabel, wmlabel, cortexlabel, pvlabel, nonbrainlabel, mix.atlas, partsmask, globalSeg);

mix.priors = Priors;
mix.atlas = Atlas;
index = find(LabeledSeg==pvlabel);
mix.LabelPVs = [mix.LabelPVs; index];
return
