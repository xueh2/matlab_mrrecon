function [data, headers] = find_and_load_cine(cine_dir, check_subdirs)
% [data, headers] = find_and_load_cine(cine_dir, check_subdirs)
% find and load cine folder

if(check_subdirs)
    [series, num] = FindSubDirs(cine_dir);
else
    series{1} = cine_dir;
    num = 1;
end

dicom_files = [];
dicom_headers = [];

for n=1:num

    if(check_subdirs)
        [dcm_files, N] = findFILE(fullfile(cine_dir, series{n}), '*.dcm');
        if(N==0)
            [dcm_files, N] = findFILE(fullfile(cine_dir, series{n}, '*.ima'));
        end
    else
        [dcm_files, N] = findFILE(cine_dir, '*.dcm');
        if(N==0)
            [dcm_files, N] = findFILE(cine_dir, '*.ima');
        end
    end
    
    for k=1:N
        header = dicominfo(dcm_files{k});
        dicom_headers = [dicom_headers; {header}];
        dicom_files = [dicom_files; dcm_files(k)];
    end
end

SliceLocations = [];
for k=1:numel(dicom_headers)
    SliceLocations(k) = dicom_headers{k}.SliceLocation;
end

SLoc = unique(SliceLocations);

SLC = numel(SLoc);
PHS = ceil(numel(dicom_headers)/SLC);
if(SLC*PHS<numel(dicom_headers))
    error('SLC*PHS~=numel(dicom_headers)');
end

SLoc_sorted = sort(SLoc, 'descend');

data_filling_record = ones(SLC, 1);
data = [];
headers = [];

for k=1:size(dicom_headers,1)
    aSLC = dicom_headers{k}.SliceLocation;
    
    slc = find(SLoc_sorted==aSLC);
    
    im = dicomread(dicom_files{k});
    RO = size(im,1);
    E1 = size(im,2);
    
    if(isempty(data))
        data = zeros(RO, E1, PHS, SLC);
        trigger_time = zeros(PHS, SLC);
        headers = cell(PHS, SLC);
    end
    
    ind = data_filling_record(slc);
    data(:,:,ind,slc) = im;
    trigger_time(ind,slc) = dicom_headers{k}.TriggerTime;
    headers{ind,slc} = dicom_headers{k};
    data_filling_record(slc) = data_filling_record(slc) + 1;
end

data_final = zeros(RO, E1, PHS, SLC);
trigger_time_final = zeros(PHS, SLC);
headers_final = cell(PHS, SLC);

for slc=1:SLC
    tt = trigger_time(:,slc);
    [tt_sorted, ind] = sort(tt);
    for phs=1:PHS
        data_final(:,:,ind(phs), slc) = data(:,:,phs,slc);
        headers_final{ind(phs), slc} = headers{phs, slc};
    end
end
   
data = data_final;
headers = headers_final;