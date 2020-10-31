% return

% input: home, labeled_slices
% the unsegmented orientation is marked by -1

cd(home);

classnumber = 4;
% classlabel = 2;

labeled_trans = Labeled_slices{i,2}(1);
labeled_sag = Labeled_slices{i,2}(2);
labeled_cor = Labeled_slices{i,2}(3);

% gm

% transversefile = 'transverse_gm.hdr';
% sagittalfile = 'sagittal_gm.hdr';
% coronaryfile = 'coronal_gm.hdr';

transverseLabel = [];
sagittalLabel = [];
coronaryLabel = [];

if ( labeled_trans ~= -1 )
%     transversefile = fullfile(manualSeg, transversefile);
    [transverseLabel, header] = LoadAnalyze(transversefile, 'Grey');
end

if ( labeled_sag ~= -1 )
%     sagittalfile = fullfile(manualSeg, sagittalfile);
    [sagittalLabel, header] = LoadAnalyze(sagittalfile, 'Grey');
end

if ( labeled_cor ~= -1 )
%     coronaryfile = fullfile(manualSeg, coronaryfile);
    [coronaryLabel, header] = LoadAnalyze(coronaryfile, 'Grey');
end

if ( Tlabel == 1 )
    templateDir = 'warped_templates';
    templateFile = 'warped_brainmask_nostem.hdr';
    tFileName = fullfile(home, templateDir, templateFile);
    [T, header] = LoadAnalyze(tFileName, 'Grey');
end

if ( Corpus_Callosum == 1 )
    templateDir = 'warped_templates';
    templateFile = 'warped_corpus_callosum.hdr';
    tFileName = fullfile(home, templateDir, templateFile);
    [Corpus, header] = LoadAnalyze(tFileName, 'Grey');
end

% Global_dsc
[Global_dsc_resultfile, Global_PVs_dsc_resultfile, Global_Local_PVs_dsc_resultfile] = ...
    CreateFileName_AssesSeg(home, Labeled_slices(i, :), classnumber, GlobalResultDir, LocalResultDir);
Global_PVs_dsc_resultfile

if ( isempty(dir(Global_PVs_dsc_resultfile)) == 0 )
    
    [Global_PVs_dsc_result, header] = LoadAnalyze(Global_PVs_dsc_resultfile, 'Grey');

    if ( Tlabel == 1 )
        Global_PVs_dsc_result(find(T==0)) = 0;
    end

    if ( Corpus_Callosum == 1 )
        Global_PVs_dsc_result(find(Corpus>0)) = 0;
        transverseLabel(find(Corpus>0)) = 0;
        sagittalLabel(find(Corpus>0)) = 0;
        coronaryLabel(find(Corpus>0)) = 0;
    end
    
    Global_PVs_dsc_4classes = computeDSC(Global_PVs_dsc_result, classlabel, ...
                transverseLabel, labeled_trans, ...
                sagittalLabel, labeled_sag, ...
                coronaryLabel, labeled_cor, ...
                header);

    voxelVolume = header.xvoxelsize * header.yvoxelsize * header.zvoxelsize;
    numLabel = length(classlabel);
    totalNum = 0;
    for tt = 1:numLabel
        index = find(Global_PVs_dsc_result==classlabel(tt));
        totalNum = totalNum + length(index);
    end
    tissueVolume = totalNum * voxelVolume;

    savedfilename = fullfile(AssessDir, 'Global_PVs_dsc_4classes.mat');
    save(savedfilename, 'Global_PVs_dsc_4classes', 'tissueVolume');
else
    error('can not find the segmenation file !');
end