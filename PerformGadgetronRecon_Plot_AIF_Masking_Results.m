
function PerformGadgetronRecon_Plot_AIF_Masking_Results(perf_cases, rest_cases, resDir, maskDir)
% PerformGadgetronRecon_Plot_AIF_Masking_Results(perf_cases, rest_cases, resDir, maskDir)

mkdir(maskDir);

case_list = [];

startN = 1
endN = size(perf_cases, 1)
for ii=startN:endN
    case_list = [case_list; {perf_cases{ii, 2}}; {perf_cases{ii, 3}}];
end

endN = size(rest_cases, 1)
for ii=startN:endN
    case_list = [case_list; {rest_cases{ii, 2}}];
end

startN = 1
endN = size(case_list, 1)
for ii=startN:endN
    disp([num2str(ii-startN+1) ' out of ' num2str(endN-startN+1) ' - ' case_list{ii,1}]);
       
    caseName = case_list{ii,1};

    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(caseName);
    caseDir = fullfile(resDir, study_dates, caseName)
    cd(caseDir)

%     aif_mask = readGTPlusExportImageSeries_Squeeze(119);

    [names, num] = findFILE('.', '*119.hdr');
    aif_mask = analyze75read(names{1});
    
    h2 = figure; imagescn(aif_mask, [], [], [], 40);
    saveas(h2, fullfile(maskDir, [caseName '_AIFMASK_GT']), 'fig');
    ff = getframe(h2);
    [X, map] = frame2im(ff);
    imwrite(X, fullfile(maskDir, [caseName '_AIFMASK_GT.bmp']), 'BMP');

    reply = input('Do you want to review aif series? Y/N [N]:','s');
    if isempty(reply)
      reply = 'N';
    end

    if(reply~='N')
        aif = readGTPlusExportImageSeries_Squeeze(101);
        h = figure; imagescn(aif(:,:,1,6:end), [], [4 15], 20);
        saveas(h, fullfile(maskDir, [caseName '_AIF_GT']), 'fig');
        ff = getframe(h);
        [X, map] = frame2im(ff);
        imwrite(X, fullfile(maskDir, [caseName '_AIF_GT.bmp']), 'BMP');       
    end

    closeall
    
    disp(['====================================================================================================']);
end
