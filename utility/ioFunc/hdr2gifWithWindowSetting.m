
function hdr2gifWithWindowSetting(hdrfile, giffile, sizeRatio, centre, width, delay)
% save the hdr as an avi file

% [data, header] = LoadAnalyze(hdrfile, GreyOrReal);
% [data, header] = Matlab_LoadAnalyze(hdrfile);
data = analyze75read(hdrfile);
header = CreateFtkHeaderInfo(data, [1 1 1 1]);
% data = single(data(:,:,5:end));
% header.sizeZ = size(data, 3);

xsize = header.sizeX;
ysize = header.sizeY;
zsize = header.sizeZ;

data = flipdim(data, 1);

minD = min(data(:));
maxD = max(data(:));
%data = normalizeImage(double(data)) .* 255;

% autonormalize
windowCentre = centre;
windowWidth = width;
lowR = windowCentre - windowWidth;
highR = windowCentre + windowWidth;

data(find(data<lowR)) = lowR;
data(find(data>highR)) = highR;

data = normalizeImage2Range(data, 0, 255);

% lowR = 0;
% highR = 255;
% mag = data; 
% mag = normalizeImage2Range(mag, lowR, highR);
% data = normalizeWindowSetting(mag, centre, width);

for kk=1:zsize
    
    slice = data(:, :, kk);

    slice = imresize(slice, sizeRatio, 'nearest');
    
    M(kk) = im2frame(uint8(slice), gray(256));
    
end    

% movie2avi(M, avifile, 'Compression', 'None', 'fps', 8);

[pathstr, name, ext] = fileparts(giffile);
gifname = fullfile(pathstr, [name '.gif']);

height = size(M(1).cdata, 1);
width = size(M(1).cdata, 2);
gifdata = zeros([height width 1 zsize], 'uint8');
for kk = 1:zsize
    gifdata(:, :, :, kk) = M(kk).cdata;
end
imwrite(gifdata, gray(256), gifname, 'gif', 'DelayTime', delay, 'LoopCount', Inf);
