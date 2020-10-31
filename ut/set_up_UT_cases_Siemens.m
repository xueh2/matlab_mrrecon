function UTCases = set_up_UT_cases_Siemens;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'T1T2',               'sasha_hc_t1t2_bh_invivo',                         'VD',   'SASHA-HC_grappa_moco_dstore.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl'   ; ...    
    'T1T2',               'sasha_hc_t1t2_fb_invivo',                         'VD',   'SASHA-HC_grappa_moco_dstore.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl'   ; ...    
    }

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
