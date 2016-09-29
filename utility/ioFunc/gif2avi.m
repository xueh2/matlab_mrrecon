
function gif2avi(giffile, avifile, sizeRatio)
% save the hdr as an avi file

[data, info] = imread(giffile);

FramesRate = 8;
global frameRate
if ( frameRate ~= -1 )
    FramesRate = frameRate;
end

M = [];
for f=1:size(data, 4)    
    M(f).cdata = data(:,:,1,f);
    M(f).colormap = info;   
end

movie2avi(M, avifile, 'Compression', 'None', 'fps', FramesRate);