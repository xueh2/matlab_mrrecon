function perf_xml = find_perf_xml(configNamePreset, case_name)
% perf_xml = find_perf_xml(configNamePreset, case_name)

if ( isempty(strfind(case_name, 'R3')))
    perf_xml = configNamePreset{1};
else
    perf_xml = configNamePreset{2};
end