
function [data, acq_time, physio_time] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep)
% read in the gtplus create images
% data = readGTPlusExportImageSeries(folderName, seriesNum);
% [data, acq_time, physio_time] = readGTPlusExportImageSeries(folderName, seriesNum, withTime, numAsRep);

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
    [names, num] = findFILE(folderName, ['GTPrep_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['CMR_' '*.hdr']);
end

if(num==0)
    [names, num] = findFILE(folderName, ['grappa_device_cpu' '*.hdr']);
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

            if(acq_time_index==0)
                N = numel(xmlContent);
                for n=1:N
                    if ( strcmp(xmlContent(n).meta(1).name, 'GT_acquisition_time_stamp') == 1 )
                        acq_time_index = n;
                    end
                    
                    if ( strcmp(xmlContent(n).meta(1).name, 'GT_physiology_time_stamp') == 1 )
                        physio_time_index = n;
                    end
                end
            end
            
            acq_time(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = str2double(xmlContent(acq_time_index).meta.value);          
            physio_time(slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = str2double(xmlContent(physio_time_index).meta.value);
        end
        
        real_name = [filename];
        real2D = analyze75read([fullfile(pathstr, real_name) '.img']);

        data(:,:,slc+1, e2+1, con+1, phs+1, rep+1, set+1, ave+1, run+1) = real2D;
    end
end



