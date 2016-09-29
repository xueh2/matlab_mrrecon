
function hdr2avi(hdrfile, avifile, sizeRatio, centre, width, delay)
% save the hdr as an avi file

[data, header] = Matlab_LoadAnalyze(hdrfile);

if ( max(data(:)) > 4095 )
    data = data ./ (max(data(:))/4095);
end
        
xsize = header.sizeX;
ysize = header.sizeY;
zsize = header.sizeZ;

minD = min(data(:));
maxD = max(data(:));

% autonormalize
windowCentre = centre;
windowWidth = width;
lowR = windowCentre - windowWidth;
highR = windowCentre + windowWidth;

data(find(data<lowR)) = lowR;
data(find(data>highR)) = highR;

data = normalizeImage2Range(data, 0, 255);

for kk=1:zsize
    
    slice = data(:, :, kk);

    %slice = imresize(slice, sizeRatio, 'nearest');
    
    slice = imresize(slice, sizeRatio, 'bicubic');
    
    M(kk) = im2frame(uint8(slice), gray(256));
    
end    

FramesRate = round(1/delay);
global frameRate
if ( frameRate ~= -1 )
    FramesRate = frameRate;
end

movie2avi(M, avifile, 'Compression', 'None', 'fps', FramesRate);

% [pathstr, name, ext] = fileparts(avifile);
% gifname = fullfile(pathstr, [name '.gif']);
% 
% height = size(M(1).cdata, 1);
% width = size(M(1).cdata, 2);
% gifdata = zeros([height width 1 zsize], 'uint8');
% for kk = 1:zsize
%     gifdata(:, :, :, kk) = M(kk).cdata;
% end
% imwrite(gifdata, gray(256), gifname, 'gif', 'DelayTime', delay, 'LoopCount', Inf);
