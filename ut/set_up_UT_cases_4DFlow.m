function UTCases = set_up_UT_cases_4DFlow;
% run the gt recon

set_UT_Dir('D');

UTCases = {

    '4DFlow',  'meas_MID00061_FID27110_4DFlow_retro_ePAT2_64sl',  'VD',   'GTPrep_3DT_RetroGated_4DFlow.xml',             'grappa_res',              'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
