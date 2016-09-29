function UTCases = set_up_UT_cases_Test;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    't1_dixon',                  '20160318_12h30m20s_6812',                                     'VD',   'Generic_Cartesian_Grappa.xml',                                  'grappa_res',              'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
