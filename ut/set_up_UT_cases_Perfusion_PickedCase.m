function UTCases = set_up_UT_cases_Perfusion_PickedCase;
% run the gt recon

% set_UT_Dir('E')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {   

    %'/home/xueh/raw_data/data_local/perf/',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_46190_0000024082200535330500000018_0000024082200535330500000018_624_20240822-223054',  'VE',   'GT_QPerf_AI_STCNNT_on_MOCO_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            

    '/data/raw_data/qperf_rawfile_testcases/',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_182758_211151816.29669997_211151816.29669997_372_20210310-144940',  'NX20',   'GT_QPerf_AI_STCNNT_on_MOCO_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            
                         

    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
