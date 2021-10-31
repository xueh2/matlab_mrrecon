
function [data, header, attribs, acq_time, physio_time] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep)
% read in the gtplus create images
% data = readGTPlusExportImageSeries_h5(folderName, seriesNum);
% [data, header, acq_time, physio_time] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep)

if nargin < 3
    withTime = 0;
end

if nargin < 4
    numAsRep = 0;
end


[name, num] = findFILE(folderName, 'ref*.h5');

h5_ind = 1;
if(num>1)
    % get the latest one    
    fnum = 0;
    for n=1:num
        [fpath, fname, ext] = fileparts(name{n});
        fdate = datenum(fname(4:end),'yyyymmdd');
        if(fdate>fnum)
            fnum = fdate;
            h5_ind = n;
        end
    end
end

info = h5info(name{h5_ind});

N = numel(info.Groups.Groups);

for n=1:N
    gname = info.Groups.Groups(n).Name;
    
    ind = find(gname=='_');
    s_num = str2num(gname(ind(end)+1:end));
    
    if(s_num == seriesNum)
        
        try
            data = h5read(name{h5_ind}, [gname '/data']);
            header = h5read(name{h5_ind}, [gname '/header']);            
            attrib = h5read(name{h5_ind}, [gname '/attributes']);
            
            s = size(squeeze(data));
            s_h = s(3:end);
            
            acq_time = header.acquisition_time_stamp;
            acq_time = double(acq_time(:))*2.5;
            if(~isempty(s_h) & numel(s_h)>1)
                acq_time = reshape(acq_time, s_h);
            end
            
            physio_time = header(:).physiology_time_stamp(1,:);
            physio_time = double(physio_time(:))*2.5;
            if(~isempty(s_h) & numel(s_h)>1)
                physio_time = reshape(physio_time, s_h);
            end
            
            maxSLC = max(header.slice);
            maxAVE = max(header.average);
            maxPHS = max(header.phase);
            maxCON = max(header.contrast);
            maxREP = max(header.repetition);
            maxSET = max(header.set);
 
            if(isreal(data))
                RO = size(data, 1);
                E1 = size(data, 2);
                max_num = size(data, ndims(data));
                if(ndims(data)==2)
                    max_num = 1;
                end
                
                im = zeros([RO E1 maxSLC+1 size(data, 3) maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1]);
            else
                RO = size(data.real, 1);
                E1 = size(data.real, 2);
                max_num = size(data.real, ndims(data.real));
                if(ndims(data.real)==2)
                    max_num = 1;
                end
                
                im = zeros([RO E1 maxSLC+1 size(data, 3) maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1]);
                im = complex(im, im);
            end
            headers = cell(maxSLC+1, maxCON+1, maxPHS+1, maxREP+1, maxSET+1, maxAVE+1);
            
            headers = struct('PatientPosition', [0 0 0], 'FOV', [0 0 0], 'phase_dir', [0 0 0], 'read_dir', [0 0 0], 'slice_dir', [0 0 0], 'window_center', -1, 'window_width', -1, 'TI', 0, 'TE', 0, 'TS', 0, 'slice_location', -1, ...
                'version', 0, 'data_type', 0, 'flags', 0, 'measurement_uid', 0, 'matrix_size', [1 1 1], 'field_of_view', [0 0 0], 'channels', 0, 'position', [0 0 0], 'patient_table_position', [0 0 0], ...
                'average', 0, 'slice', 0, 'contrast', 0, 'phase', 0, 'repetition', 0, 'set', 0, 'acquisition_time_stamp', 0, 'physiology_time_stamp', [0 0 0], 'image_type', 0, ...
                'image_index', 0, 'image_series_index', 0, 'user_int', zeros(1,8), 'user_float', zeros(1,8), 'attribute_string_len', 0);
            
            attribs = cell(maxSLC+1, maxCON+1, maxPHS+1, maxREP+1, maxSET+1, maxAVE+1);
            
            for n=1:max_num
                if(isreal(data))
                    im(:,:, header.slice(n)+1, :, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1) = data(:,:,:,n);
                else
                    im(:,:, header.slice(n)+1, :, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1) = complex(data.real(:,:,:,n), data.imag(:,:,:,n));
                end
                h = extract_header(header, n, max_num);
                headers(header.slice(n)+1, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1) = h;
                
                attribs{header.slice(n)+1, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1} = attrib{n};
            end
            
            data = im;
            header = headers;
            
        catch
            disp(['error happened in read in ' gname]);
            data = [];
            header = [];
            acq_time = [];
            physio_time = [];
            attribs = [];
        end
        
    end
end

end

function h = extract_header(headers, n, max_num)
    fields = fieldnames(headers);
    
    h  = [];
    for k=1:numel(fields)
        v = getfield(headers, fields{k});
        v = squeeze(v);
        if(size(v, 1)==max_num)
            h = setfield(h, fields{k}, v(n));
        else
            h = setfield(h, fields{k}, v(:, n));
        end
    end
    
    h.PatientPosition = h.position;
    h.FOV = h.field_of_view;
    h.window_center = -1;
    h.window_width = -1;
    h.TI = 0;    
    h.TE = 0;
    h.TS = 0;
    h.slice_location = dot(h.PatientPosition, h.slice_dir);
end
