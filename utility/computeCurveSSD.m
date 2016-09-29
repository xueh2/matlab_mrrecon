function [ssd, m, s] = computeCurveSSD(gtIm, testIm, pts)
% compute the ssd cost of ground truth images and tested images

numOfPt = size(pts, 1);

ssd = size(numOfPt, 1);

for p=1:numOfPt
    currPt = pts(p, :);
    d = gtIm(currPt(2), currPt(1), :) - testIm(currPt(2), currPt(1), :);
    ssd(p) = (d(:)'*d(:))/numOfPt;
end

m = mean(ssd);
s = std(ssd);