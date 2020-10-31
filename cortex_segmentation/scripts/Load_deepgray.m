
% compute DSC for deep gray matter

cd(home);

classlabel = 1;

labeled_trans = Labeled_slices{i,2}(1);
labeled_sag = Labeled_slices{i,2}(2);
labeled_cor = Labeled_slices{i,2}(3);

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


deepgray_dsc_resultfile = fullfile(manualSeg, deepgrayfile);

[deepgray_dsc_result, header] = LoadAnalyze(deepgray_dsc_resultfile, 'Grey');
deepgray_dsc = computeDSC(deepgray_dsc_result, classlabel, ...
            transverseLabel, labeled_trans, ...
            sagittalLabel, labeled_sag, ...
            coronaryLabel, labeled_cor, ...
            header);

savedfilename = fullfile(AssessDir, 'deepgray_dsc.mat');
save(savedfilename, 'deepgray_dsc');