
function data = readGTPlusExportCxImages(folderName, dataRole, dataRole2)
% read in the gtplus create images
% data = readGTPlusExportCxImages(folderName, dataRole, dataRole2)

if nargin == 2 
    [names, num] = findFILE(folderName, [dataRole '_*_attrib.xml']);
end

if nargin == 3 
    [names, num] = findFILE(folderName, [dataRole '*' dataRole2 '_*_attrib.xml']);
end

maxSLC = 0;
maxE2 = 0;
maxCON = 0;
maxPHS = 0;
maxREP = 0;
maxSET = 0;
maxAVE = 0;
maxCHA = 0;
maxRUN = 0;

for ii=1:num
    name = names{ii};

    [pathstr, filename, ext] = fileparts(name);   
    [slc, e2, con, phs, rep, set, ave, cha, run] = getImageISMRMRDDim(filename);
    
    if ( slc > maxSLC ) maxSLC = slc; end
    if ( e2 > maxE2 ) maxE2 = e2; end
    if ( con > maxCON ) maxCON = con; end
    if ( phs > maxPHS ) maxPHS = phs; end
    if ( rep > maxREP ) maxREP = rep; end
    if ( set > maxSET ) maxSET = set; end   
    if ( ave > maxAVE ) maxAVE = ave; end   
    if ( cha > maxCHA ) maxCHA = cha; end   
    if ( run > maxRUN ) maxRUN = run; end   
end

[pathstr, filename, ext] = fileparts(names{1});
real_name = fullfile(pathstr, [filename(1:end-6) 'REAL.hdr']);
try
    data2D = analyze75read(real_name);
catch
    real_name = fullfile(pathstr, [filename(1:end-6) 'REAL']);
    [data2D, header] = Matlab_gt_read_nifti(real_name);
end
RO = size(data2D, 1);
E1 = size(data2D, 2);
data = zeros([RO E1 maxSLC+1 maxE2+1 maxCON+1 maxPHS+1 maxREP+1 maxSET+1 maxAVE+1  maxCHA+1 maxRUN+1]);

for ii=1:num
    name = names{ii};
    
    [pathstr, filename, ext] = fileparts(name);
    
    disp(['--> load ' filename]);
    
    [slc, e2, con, phs, rep, set, ave, cha, run] = getImageISMRMRDDim(filename);   
      
    try
        real_name = [filename(1:end-6) 'REAL.hdr'];
        imag_name = [filename(1:end-6) 'IMAG.hdr'];

        real2D = analyze75read(fullfile(pathstr, real_name));   
        imag2D = analyze75read(fullfile(pathstr, imag_name));   
    catch
        real_name = [filename(1:end-6) 'REAL'];
        imag_name = [filename(1:end-6) 'IMAG'];

        [real2D, header] = Matlab_gt_read_nifti(fullfile(pathstr, real_name));   
        [imag2D, header] = Matlab_gt_read_nifti(fullfile(pathstr, imag_name));   
    end
    
    data(:,:,slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, cha+1, run+1) = complex(real2D, imag2D);    
end



