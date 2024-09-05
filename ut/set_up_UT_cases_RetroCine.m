function UTCases = set_up_UT_cases_RetroCine;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {

    '/data2/raw_data/data_local/retrocine',  'Retro_Lin_Cine_2DT_LAX_GLS_000000_4146679_4146688_231_00000000-000000',  'VE',   'GT_RetroCine_gmap_augmentation.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...
    
  };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
