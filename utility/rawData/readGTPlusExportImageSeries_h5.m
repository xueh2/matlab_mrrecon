
function [data, header, acq_time, physio_time] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep)
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

info = h5info(name{1});

g1 = info.Groups.Name;

N = numel(info.Groups.Groups);

for n=1:N
    gname = info.Groups.Groups(n).Name;
    
    ind = find(gname=='_');
    s_num = str2num(gname(ind(end)+1:end));
    
    if(s_num == seriesNum)
        
        try
            data = h5read(name{1}, [gname '/data']);
            header = h5read(name{1}, [gname '/header']);            
            
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
            
        catch
            disp(['error happened in read in ' gname]);
            data = [];
            header = [];
            acq_time = [];
            physio_time = [];
        end
        
    end
end

