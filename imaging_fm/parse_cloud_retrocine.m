function [patient_id, study_id, meas_id] = parse_cloud_retrocine(res_dir, str_key)
% [patient_id, study_id, meas_id] = parse_cloud_retrocine(res_dir, str_key)

if nargin < 2
    str_key = '*gmap.mrd';
end

[names, num] = findFILE(res_dir, str_key);

patient_id = [];
study_id = [];
meas_id = [];

for n=1:num

    [fpath, fname, ext] = fileparts(names{n});
    
    ind = strfind(fname, '_');

    if str2num(fname(ind(5)+1:ind(6)-1)) == 0
        patient_id = [patient_id; {fname(ind(6)+1:ind(7)-1)}];
        study_id = [study_id; {fname(ind(7)+1:ind(8)-1)}];
        meas_id = [meas_id; {fname(ind(8)+1:ind(9)-1)}];
    else
        patient_id = [patient_id; {fname(ind(5)+1:ind(6)-1)}];
        study_id = [study_id; {fname(ind(6)+1:ind(7)-1)}];
        meas_id = [meas_id; {fname(ind(7)+1:ind(8)-1)}];
    end
end
