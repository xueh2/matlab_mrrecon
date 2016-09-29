
function dataScaled = normalizeWindowSetting(data, centre, width)

data = single(data);

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

dataScaled = data;