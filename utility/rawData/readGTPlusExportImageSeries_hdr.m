
function [data, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_hdr(folderName, seriesNum, withTime, numAsRep)
% read in the gtplus create images
% data = readGTPlusExportImageSeries_hdr(folderName, seriesNum);
% [data, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_hdr(folderName, seriesNum, withTime, numAsRep);

if nargin < 3
    withTime = 0;
end

if nargin < 4
    numAsRep = 0;
end

[names, num] = findFILE(folderName, ['gadgetron_SLC' '*.hdr']);
% [names, num] = findFILE(folderName, [dataRole '*.hdr']);

if(num==0)
    [names, num] = findFILE(folderName, ['Generic_Cartesian' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_2DT_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_3DT_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GT_2DT' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GT_3DT' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['Generic_Spiral.xml' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['CMR_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['grappa_device_cpu' '*.hdr']);
end
if(num==0)
    [names, num] = findFILE(folderName, ['grappa_cpu' '*.hdr']);
end
if(num==0)
    [names, num] = findFILE(folderName, ['epi' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_2DT_Perfusion_AIF_TwoEchoe' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_2DT_Perf_AIF_2E' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_2DT_Perf_AIFR3_2E' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['BART' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['SASHA-HC' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['Generic_2DT_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['GTPrep_2DT_RetroGated_Cine_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['results_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['Generic_RTCine' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['Siemens_Gadgetron' '*.hdr']);
end

maxSLC = 0;
maxE2 = 0;
maxCON = 0;
maxPHS = 0;
maxREP = 0;
maxSET = 0;
maxAVE = 0;
maxRUN = 0;
maxImageNum = 0;

endo_pt = [];
epi_pt = [];

header = struct('PatientPosition', [0 0 0], 'FOV', [0 0 0], 'phase_dir', [0 0 0], 'read_dir', [0 0 0], 'slice_dir', [0 0 0], 'window_center', -1, 'window_width', -1, 'TI', 0, 'slice_location', -1);

hasImageSize = 0;
for ii=1:num
    name = names{ii};
    [pathstr, filename, ext] = fileparts(name);
    
    [slc, e2, con, phs, rep, set, ave, cha, run, imagenum, series] = getImageISMRMRDDim(filename);
    
    if ( series == seriesNum )
        
        if ( slc > maxSLC ) maxSLC = slc; end
        if ( e2 > maxE2 ) maxE2 = e2; end
        if ( con > maxCON ) maxCON = con; end
        if ( phs > maxPHS ) maxPHS = phs; end
        if ( rep > maxREP ) maxREP = rep; end
        if ( set > maxSET ) maxSET = set; end
        if ( ave > maxAVE ) maxAVE = ave; end
        if ( run > maxRUN ) maxRUN = run; end   
        if ( imagenum > maxImageNum ) maxImageNum = imagenum; end   
        
        if ( ~hasImageSize )                        
            real_name = fullfile(pathstr, filename);
            data2D = analyze75read([real_name '.hdr']);

            RO = size(data2D, 1);
            E1 = size(data2D, 2);
            
            hasImageSize = 1;
        end
    end
end

maxImageNum = maxImageNum -1;

if (maxSLC+maxE2+maxCON+maxPHS+maxREP+maxSET+maxAVE==0)
    numAsRep = 1;
    maxREP = maxImageNum;
end

if(numAsRep)
    maxREP = maxImageNum;
end

data = zeros([RO E1 maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1]);
acq_time = zeros([maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1]);
physio_time = zeros([maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1]);

user_int = zeros([maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1 8]);

acq_time_index = 0;
physio_time_index = 0;

for ii=1:num
    name = names{ii};
    [pathstr, filename, ext] = fileparts(name);
    
    [slc, e2, con, phs, rep, set, ave, cha, run, imagenum, series] = getImageISMRMRDDim(filename);
    
    if (numAsRep)
        rep = imagenum-1;
    end
    
    if ( series == seriesNum )

        if(withTime)
            xmlContent = xml_load(fullfile(pathstr, [filename '.attrib']));

            endo = 0;
            epi = 0;
            PatientPosition = 0;
            FOV = 0;
            phase_dir = 0;
            read_dir = 0;
            slice_dir = 0;
            window_center = -1;
            window_width = -1;
            TI = -1;
            
%             if(acq_time_index==0)
                N = numel(xmlContent);
                for n=1:N
                    if ( strcmp(xmlContent(n).meta(1).name, 'GT_acquisition_time_stamp') == 1 )
                        acq_time_index = n;
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'acquisition_time_stamp') == 1 )
                        acq_time_index = n;
                    end
                    
                    if ( strcmp(xmlContent(n).meta(1).name, 'GT_physiology_time_stamp') == 1 )
                        physio_time_index = n;
                    end
                    
                    if ( strcmp(xmlContent(n).meta(1).name, 'ENDO') == 1 )
                        endo = n;
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'EPI') == 1 )
                        epi = n;
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'PatientPosition') == 1 )
                        PatientPosition = n;
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'FOV') == 1 )
                        FOV = n;
                    end
%                     if ( strcmp(xmlContent(n).meta(1).name, 'recon_FOV') == 1 )
%                         FOV = n;
%                     end
                    if ( strcmp(xmlContent(n).meta(1).name, 'phase_dir') == 1 )
                        phase_dir = n;
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'read_dir') == 1 )
                        read_dir = n;
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'slice_dir') == 1 )
                        slice_dir = n;
                    end
                    
                    if ( ~isempty(strfind(xmlContent(n).meta(1).name, 'GT_ROI')) )
                        sector_no = str2double(xmlContent(n).meta(1).name(end));                        
                    end
                    
                    for kk=0:7
                        user_int_str = ['user_int_' num2str(kk)];
                        if ( strcmp(xmlContent(n).meta(1).name, user_int_str) == 1 )
                            try
                                user_int(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1, kk+1) = str2double(xmlContent(n).meta.value);  
                            catch
                            end
                        end
                    end
                    
                    if ( strcmp(xmlContent(n).meta(1).name, 'GADGETRON_WindowCenter') == 1 )
                        window_center = str2double(xmlContent(n).meta.value);
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'GADGETRON_WindowWidth') == 1 )
                        window_width = str2double(xmlContent(n).meta.value);
                    end
                    if ( strcmp(xmlContent(n).meta(1).name, 'GADGETRON_TI') == 1 )
                        TI = str2double(xmlContent(n).meta.value);
                    end
                end
%             end
            
            try
                acq_time(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = str2double(xmlContent(acq_time_index).meta.value);          
                physio_time(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = str2double(xmlContent(physio_time_index).meta.value);
            catch
            end
            
            curr_endo_pt = [];
            if(endo>0)
                num_pt = numel(xmlContent(endo).meta);
                for n=1:2:num_pt
                    curr_endo_pt = [curr_endo_pt; str2double(xmlContent(endo).meta(n).value) str2double(xmlContent(endo).meta(n+1).value)];
                end
            end
            
            curr_epi_pt = [];
            if(epi>0)
                num_pt = numel(xmlContent(epi).meta);
                for n=1:2:num_pt
                    curr_epi_pt = [curr_epi_pt; str2double(xmlContent(epi).meta(n).value) str2double(xmlContent(epi).meta(n+1).value)];
                end
            end
            
            endo_pt = [endo_pt; {phs slc curr_endo_pt}];
            epi_pt = [epi_pt; {phs slc curr_epi_pt}];
            
            curr_header = struct('PatientPosition', [0 0 0], 'FOV', [0 0 0], 'phase_dir', [0 0 0], 'read_dir', [0 0 0], 'slice_dir', [0 0 0], 'window_center', -1, 'window_width', -1, 'TI', 0, 'slice_location', -1);
            
            if(PatientPosition>0)
                num_pt = numel(xmlContent(PatientPosition).meta);
                curr_header.PatientPosition = [str2double(xmlContent(PatientPosition).meta(1).value) str2double(xmlContent(PatientPosition).meta(2).value) str2double(xmlContent(PatientPosition).meta(3).value)];
            end
            if(FOV>0)
                num_pt = numel(xmlContent(FOV).meta);
                curr_header.FOV = [str2double(xmlContent(FOV).meta(1).value) str2double(xmlContent(FOV).meta(2).value) str2double(xmlContent(FOV).meta(3).value)];
            end
            if(phase_dir>0)
                num_pt = numel(xmlContent(phase_dir).meta);
                curr_header.phase_dir = [str2double(xmlContent(phase_dir).meta(1).value) str2double(xmlContent(phase_dir).meta(2).value) str2double(xmlContent(phase_dir).meta(3).value)];
            end
            if(read_dir>0)
                num_pt = numel(xmlContent(read_dir).meta);
                curr_header.read_dir = [str2double(xmlContent(read_dir).meta(1).value) str2double(xmlContent(read_dir).meta(2).value) str2double(xmlContent(read_dir).meta(3).value)];
            end
            if(slice_dir>0)
                num_pt = numel(xmlContent(slice_dir).meta);
                curr_header.slice_dir = [str2double(xmlContent(slice_dir).meta(1).value) str2double(xmlContent(slice_dir).meta(2).value) str2double(xmlContent(slice_dir).meta(3).value)];
                
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
            
            header(slc+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = curr_header;
        end
        
        real_name = [filename];
        real2D = analyze75read([fullfile(pathstr, real_name) '.img']);

        try
            data(:,:,slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = real2D;
        catch
            data(:,:,slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = permute(real2D, [2 1]);
        end
    end
end

if(~isempty(endo_pt))
    N = size(endo_pt, 1);
    ind = zeros(N, 1);
    for n=1:N
        ind(n) = endo_pt{n, 1};
    end
    
    [ind_sorted, ind_order] = sort(ind, 'ascend');
    
    endo_pt_sorted = endo_pt;
    
    for n=1:N
        endo_pt_sorted{n, 1} = endo_pt{ind_order(n), 1};
        endo_pt_sorted{n, 2} = endo_pt{ind_order(n), 2};
        endo_pt_sorted{n, 3} = endo_pt{ind_order(n), 3};
    end
    endo_pt = endo_pt_sorted;
end

if(~isempty(epi_pt))
    N = size(epi_pt, 1);
    ind = zeros(N, 1);
    for n=1:N
        ind(n) = epi_pt{n, 1};
    end
    
    [ind_sorted, ind_order] = sort(ind, 'ascend');
    
    epi_pt_sorted = epi_pt;
    
    for n=1:N
        epi_pt_sorted{n, 1} = epi_pt{ind_order(n), 1};
        epi_pt_sorted{n, 2} = epi_pt{ind_order(n), 2};
        epi_pt_sorted{n, 3} = epi_pt{ind_order(n), 3};
    end
    epi_pt = epi_pt_sorted;
end
