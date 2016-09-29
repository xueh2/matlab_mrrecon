
function [sampledValues, pos, h_out] = Measure_Sharpness(data, header, stepSize, centre, width, pos, h, deleteFrame)
% This is a function to measure the sharpness of a myocardium
  
lowR = 0;
highR = 255;
mag = data; 
mag = normalizeImage2Range(mag, lowR, highR);
mag = normalizeWindowSetting(mag, centre, width);

if ( isempty(pos) )
    if ( ishandle(h) )
        figure(h); imshow(mag, []);
        h_out = h;
    else
        h_out = figure; imshow(mag, []);
    end
    L = imline;
    pos = L.getPosition;
    if ( deleteFrame ) close(h); end
else
    h_out = h;
end

stepX = stepSize;
stepY = stepX*(pos(2,2)-pos(1,2))/abs(pos(2,1)-pos(1,1));

x = pos(1,1):stepX:pos(2,1);
y = pos(1,2):stepY:pos(2,2);

sampledValues = interp2(data, x,y, 'linear');


