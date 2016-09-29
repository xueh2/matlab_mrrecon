
%% GT UT all tests

% setenv('GT_HOST', 'localhost'); setenv('GT_PORT', '9002');
% setenv('GT_HOST', 'barbados'); setenv('GT_PORT', '9006');
% setenv('GT_HOST', 'barbados'); setenv('GT_PORT', '9002');
% setenv('GT_HOST', '137.187.135.57'); setenv('GT_PORT', '9002');
% setenv('GT_HOST', '52.1.38.252'); setenv('GT_PORT', '9002');
setenv('GT_HOST', 'denmark'); setenv('GT_PORT', '9002');

UTCollectRef = 0;

setenv('GTPLUS_UT_DIR', 'E:\gtuser\gt_windows_setup\ut')
setenv('GTPLUS_UT_DIR', 'F:\data');
setenv('GTPLUS_UT_DIR', 'J:\ut');
setenv('GTPLUS_UT_DIR', '\\137.187.134.135\share\ut');


timeUsed_gtprep_Win = performUTValidation(set_up_UT_cases_ScannerUpdate, UTCollectRef, 0, 'localhost', '9002')
performUTValidation(set_up_UT_cases_ScannerUpdate, UTCollectRef, 0, 'barbados.nhlbi.nih.gov', '9006')

performUTValidation(set_up_UT_cases_Cloud, UTCollectRef, 0, 'localhost', '9002')

% amazon test
performUTValidation(set_up_UT_cases_Cloud, UTCollectRef, 0, 'denmark', '9222', 4)
performUTValidation(set_up_UT_cases_Cloud, UTCollectRef, 0, 'denmark', '9222', 31)
performUTValidation(set_up_UT_cases_Cloud, UTCollectRef, 0, 'denmark', '9222', 32)


timeUsed = performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'localhost', '9002')

timeUsed = performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'localhost', '9002', 6)

timeUsed = performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'denmark', '9002', 6)
performUTValidation(set_up_UT_cases_GT_UT, 0, 0, 'barbados', '9002', 6, 6, 1, 0);

performUTValidation(set_up_UT_cases_GT_UT_GPU, UTCollectRef, 0, 'barbados', '9002')


performUTValidation(set_up_UT_cases_Test, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_Cine, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_Binning, UTCollectRef, 0, 'localhost', '9002')


timeUsed_Win = performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'localhost', '9002');

timeUsed_Win = performUTValidation(set_up_UT_cases_SNRUnitRecon, UTCollectRef, 0, 'localhost', '9002');

performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'localhost', '9002', 5)

performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'denmark', '9002', 3)

performUTValidation(set_up_UT_cases_ExternalSites, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_LGE, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_Perfusion, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_3D, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_CVSparse, UTCollectRef, 0, 'localhost', '9002')

performUTValidation(set_up_UT_cases_EPI, UTCollectRef, 0, 'localhost', '9002')


performUTValidation(set_up_UT_cases_FOVTest, 1, 0, 'localhost', '9002')


performUTValidationICE(set_up_UT_cases_ICE, UTCollectRef)










plotTimming_GT_UT(timeUsed_Win, timeUsed_nauru)
plotTimming_GT_UT(timeUsed_Win, timeUsed_nauru_mkl)

timeUsed_Win = performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'localhost', '9002');
timeUsed_Win_baseline = timeUsed_Win;

timeUsed_Win = performUTValidation(set_up_UT_cases_GT_UT, UTCollectRef, 0, 'localhost', '9002');
plotTimming_GT_UT(timeUsed_Win_baseline, timeUsed_Win)

timeUsed_gtprep_Win = performUTValidation(set_up_UT_cases_ScannerUpdate, UTCollectRef, 0, 'localhost', '9002')
plotTimming_GT_UT(timeUsed_gtprep_Win_baseline, timeUsed_gtprep_Win)

data = analyze75read('PD_row0.hdr');

