function UTCases = set_up_UT_cases_XA31;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {    
    'XA31',           'CMR_Localizer_ECG_Training_183206_30000022021114365275500000003_30000022021114365275500000003_91_20220214-143851',  'NX',   'Generic_Cartesian_Grappa_CoilQA_ECG.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_NX.xsl'   ; ...
    'XA31',           'RT_Cine_LIN_183206_30000022021114365275500000003_30000022021114365275500000003_92_20220214-144204',  'NX',   'Generic_RTCine_ECG.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_NX.xsl'   ; ...
    
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end


