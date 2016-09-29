
function [data, acq_time, physio_time] = readGTPlusExportImages(folderName, dataRole)
% read in the gtplus create images

[names, num] = findFILE(folderName, [dataRole '*.attrib']);
% [names, num] = findFILE(folderName, [dataRole '*.hdr']);

maxSLC = 0;
maxE2 = 0;
maxCON = 0;
maxPHS = 0;
maxREP = 0;
maxSET = 0;
maxAVE = 0;
maxRUN = 0;

for ii=1:num
    name = names{ii};
    
    [pathstr, filename, ext] = fileparts(name);
    [slc, e2, con, phs, rep, set, ave, run] = getImageISMRMRDDim(filename);
    
    if ( slc > maxSLC ) maxSLC = slc; end
    if ( e2 > maxE2 ) maxE2 = e2; end
    if ( con > maxCON ) maxCON = con; end
    if ( phs > maxPHS ) maxPHS = phs; end
    if ( rep > maxREP ) maxREP = rep; end
    if ( set > maxSET ) maxSET = set; end
    if ( ave > maxAVE ) maxAVE = ave; end
    if ( run > maxRUN ) maxRUN = run; end   
end

[pathstr, filename, ext] = fileparts(names{1});

% real_name = fullfile(pathstr, [filename(1:end-7) '.hdr']);
% data2D = analyze75read(real_name);

real_name = fullfile(pathstr, filename);
try
    [data2D, header] = Matlab_gt_read_nifti(real_name);
catch
    % [data2D, header] = Matlab_gt_read_analyze(real_name);
    data2D = analyze75read([real_name '.img']);
end

RO = size(data2D, 1);
E1 = size(data2D, 2);
data = zeros([RO E1 maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1]);
acq_time = zeros([maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1]);
physio_time = zeros([maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1 maxRUN+1]);

for ii=1:num
    name = names{ii};
    
    xmlContent = xml_load(name);
   
    [pathstr, filename, ext] = fileparts(name);
    [slc, e2, con, phs, rep, set, ave, run] = getImageISMRMRDDim(filename);   

    if ( strcmp(xmlContent(11).meta(1).name, 'AcquisitionTime') == 1 )
        acq_time(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = str2double(xmlContent(11).meta.value);
        physio_time(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = str2double(xmlContent(12).meta.value);
    end
    
    [pathstr, img_name, ext] = fileparts(name);   
    
    % real_name = [img_name(1:end-7) '.hdr'];
    % real2D = analyze75read(fullfile(pathstr, real_name));   
    
    real_name = [img_name];
    try
        [real2D, header] = Matlab_gt_read_nifti(fullfile(pathstr, real_name));
    catch
        % [real2D, header] = Matlab_gt_read_analyze(fullfile(pathstr, real_name));
        real2D = analyze75read([fullfile(pathstr, real_name) '.img']);
    end
    
    data(:,:,slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = real2D;
end



