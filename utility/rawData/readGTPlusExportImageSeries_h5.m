
function [data, header, attribs, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep)
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

endo_pt = [];
epi_pt = [];

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
            
%             headers = struct('PatientPosition', [0 0 0], 'FOV', [0 0 0], 'phase_dir', [0 0 0], 'read_dir', [0 0 0], 'slice_dir', [0 0 0], 'window_center', -1, 'window_width', -1, 'TI', 0, 'TE', 0, 'TS', 0, 'slice_location', -1, ...
%                 'version', 0, 'data_type', 0, 'flags', 0, 'measurement_uid', 0, 'matrix_size', [1 1 1], 'field_of_view', [0 0 0], 'channels', 0, 'position', [0 0 0], 'patient_table_position', [0 0 0], ...
%                 'average', 0, 'slice', 0, 'contrast', 0, 'phase', 0, 'repetition', 0, 'set', 0, 'acquisition_time_stamp', 0, 'physiology_time_stamp', [0 0 0], 'image_type', 0, ...
%                 'image_index', 0, 'image_series_index', 0, 'user_int', zeros(1,8), 'user_float', zeros(1,8), 'attribute_string_len', 0);
            
            attribs = cell(maxSLC+1, maxCON+1, maxPHS+1, maxREP+1, maxSET+1, maxAVE+1);
            
            for n=1:max_num
                if(isreal(data))
                    im(:,:, header.slice(n)+1, :, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1) = data(:,:,:,n);
                else
                    im(:,:, header.slice(n)+1, :, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1) = complex(data.real(:,:,:,n), data.imag(:,:,:,n));
                end
                h = extract_header(header, n, max_num);
                headers{header.slice(n)+1, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1} = h;
                
                attribs{header.slice(n)+1, header.contrast(n)+1, header.phase(n)+1, header.repetition(n)+1, header.set(n)+1, header.average(n)+1} = attrib{n};
            end
            
            % set the header
            for tt=1:max_num
                endo = 0;
                epi = 0;
                PatientPosition = 0;
                FOV = 0;
                recon_FOV = 0;
                phase_dir = 0;
                read_dir = 0;
                slice_dir = 0;
                window_center = -1;
                window_width = -1;
                TI = -1;
                TE = -1;
                TS = -1;

                xmlContent = gt_xml_load_from_string(attrib{tt});
                xmlContent = xmlContent.ismrmrdMeta.meta;

                N = numel(xmlContent);
                for n=1:N
                    if ( strcmp(xmlContent{n}.name.Text, 'GT_acquisition_time_stamp') == 1 )
                        acq_time_index = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'acquisition_time_stamp') == 1 )
                        acq_time_index = n;
                    end

                    if ( strcmp(xmlContent{n}.name.Text, 'GT_physiology_time_stamp') == 1 )
                        physio_time_index = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'physiology_time_stamp') == 1 )
                        physio_time_index = n;
                    end

                    if ( strcmp(xmlContent{n}.name.Text, 'ENDO') == 1 )
                        endo = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'EPI') == 1 )
                        epi = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'PatientPosition') == 1 )
                        PatientPosition = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'FOV') == 1 )
                        FOV = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'recon_FOV') == 1 )
                        recon_FOV = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'phase_dir') == 1 )
                        phase_dir = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'read_dir') == 1 )
                        read_dir = n;
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'slice_dir') == 1 )
                        slice_dir = n;
                    end

                    for kk=0:7
                        user_int_str = ['user_int_' num2str(kk)];
                        if ( strcmp(xmlContent{n}.name.Text, user_int_str) == 1 )
                            try
                                user_int(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1, kk+1) = str2double(xmlContent{n}.value{1}.Text);  
                            catch
                            end
                        end
                    end

                    if ( strcmp(xmlContent{n}.name.Text, 'GADGETRON_WindowCenter') == 1 )
                        if(numel(xmlContent{n}.value)>1)
                            window_center = str2double(xmlContent{n}.value{1}.Text);
                        else
                            window_center = str2double(xmlContent{n}.value.Text);
                        end
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'GADGETRON_WindowWidth') == 1 )
                        if(numel(xmlContent{n}.value)>1)
                            window_width = str2double(xmlContent{n}.value{1}.Text);
                        else
                            window_width = str2double(xmlContent{n}.value.Text);
                        end
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'GADGETRON_TI') == 1 )
                        if(numel(xmlContent{n}.value)>1)
                            TI = str2double(xmlContent{n}.value{1}.Text);
                        else
                            TI = str2double(xmlContent{n}.value.Text);
                        end
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'GADGETRON_TE') == 1 )
                        if(numel(xmlContent{n}.value)>1)
                            TE = str2double(xmlContent{n}.value{1}.Text);
                        else
                            TE = str2double(xmlContent{n}.value.Text);
                        end
                    end
                    if ( strcmp(xmlContent{n}.name.Text, 'GADGETRON_TS') == 1 )
                        if(numel(xmlContent{n}.value)>1)
                            TS = str2double(xmlContent{n}.value{1}.Text);
                        else
                            TS = str2double(xmlContent{n}.value.Text);
                        end
                    end
                end

                try
                    if(iscell(xmlContent{acq_time_index}.value))
                        acq_time(header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1) = str2double(xmlContent{acq_time_index}.value.Text);
                    else
                        acq_time(header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1) = str2double(xmlContent{acq_time_index}.value.Text);          
                    end
                    physio_time(header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1) = str2double(xmlContent{physio_time_index}.value.Text);
                catch
                end

                curr_endo_pt = [];
                if(endo>0)
                    num_pt = numel(xmlContent{endo}.value);
                    for n=1:2:num_pt
                        curr_endo_pt = [curr_endo_pt; str2double(xmlContent{endo}.value{n}.Text) str2double(xmlContent{endo}.value{n+1}.Text)];
                    end
                end

                curr_epi_pt = [];
                if(epi>0)
                    num_pt = numel(xmlContent{epi}.value);
                    for n=1:2:num_pt
                        curr_epi_pt = [curr_epi_pt; str2double(xmlContent{epi}.value{n}.Text) str2double(xmlContent{epi}.value{n+1}.Text)];
                    end
                end

                endo_pt = [endo_pt; {header.phase(tt)+1 header.slice(tt)+1 curr_endo_pt}];
                epi_pt = [epi_pt; {header.phase(tt)+1 header.slice(tt)+1 curr_epi_pt}];

                curr_header = struct('PatientPosition', [0 0 0], 'FOV', [0 0 0], 'phase_dir', [0 0 0], 'read_dir', [0 0 0], 'slice_dir', [0 0 0], 'window_center', -1, 'window_width', -1, 'TI', 0, 'TE', 0, 'TS', 0, 'slice_location', -1);

                if(FOV==0)
                    FOV=recon_FOV;
                end

                if(PatientPosition>0)
                    curr_header.PatientPosition = [str2double(xmlContent{PatientPosition}.value{1}.Text) str2double(xmlContent{PatientPosition}.value{2}.Text) str2double(xmlContent{PatientPosition}.value{3}.Text)];
                end
                if(FOV>0)
                    curr_header.FOV = [str2double(xmlContent{FOV}.value{1}.Text) str2double(xmlContent{FOV}.value{2}.Text) str2double(xmlContent{FOV}.value{3}.Text)];
                end
                if(phase_dir>0)
                    curr_header.phase_dir = [str2double(xmlContent{phase_dir}.value{1}.Text) str2double(xmlContent{phase_dir}.value{2}.Text) str2double(xmlContent{phase_dir}.value{3}.Text)];
                end
                if(read_dir>0)
                    curr_header.read_dir = [str2double(xmlContent{read_dir}.value{1}.Text) str2double(xmlContent{read_dir}.value{2}.Text) str2double(xmlContent{read_dir}.value{3}.Text)];
                end
                if(slice_dir>0)
                    curr_header.slice_dir = [str2double(xmlContent{slice_dir}.value{1}.Text) str2double(xmlContent{slice_dir}.value{2}.Text) str2double(xmlContent{slice_dir}.value{3}.Text)];

                    if(PatientPosition>0)
                        curr_header.slice_location = dot(curr_header.PatientPosition, curr_header.slice_dir);
                    end
                end
                if(window_center>0)
                    curr_header.window_center = window_center;
                end
                if(window_width>0)
                    curr_header.window_width = window_width;
                end
                if(TI>0)
                    curr_header.TI = TI;
                end
                if(TE>0)
                    curr_header.TE = TE;
                end
                if(TS>0)
                    curr_header.TS = TS;
                end

                headers{header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1}.window_center = curr_header.window_center;
                headers{header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1}.window_width = curr_header.window_width;
                headers{header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1}.TI = curr_header.TI;
                headers{header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1}.TE = curr_header.TE;
                headers{header.slice(tt)+1, header.contrast(tt)+1, header.phase(tt)+1, header.repetition(tt)+1, header.set(tt)+1, header.average(tt)+1}.TS = curr_header.TS;
            end
            
            % ------------------------------------------
            
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
