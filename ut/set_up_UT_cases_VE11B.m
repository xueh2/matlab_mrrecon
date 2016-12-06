function UTCases = set_up_UT_cases_VE11B;
% run the gt recon

set_UT_Dir('D');

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {

    'VE11B',                 'meas_MID00073_FID00975_MiniFLASH_pk',              'VB',   'Generic_Cartesian_Grappa_SNR_CoilQA.xml',             'grappa_res',              'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
