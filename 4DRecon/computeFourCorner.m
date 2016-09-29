function [wc1, wc2, wc3, wc4] = computeFourCorner(data, header)
% compute the world coordinates of four corners
% wc1 ---- wc2
% |        | 
% |        |
% |        |
% wc4 ---- wc3

row = size(data, 1);
col = size(data, 2);

% compute the four corner of input images
c1 = [0 0 0];
wc1 = c1;
[wc1(1), wc1(2), wc1(3)] = Image2WorldMrFtk(header, c1(1), c1(2), c1(3));

c2 = [col-1 0 0];
wc2 = c2;
[wc2(1), wc2(2), wc2(3)] = Image2WorldMrFtk(header, c2(1), c2(2), c2(3));

c3 = [col-1 row-1 0];
wc3 = c3;
[wc3(1), wc3(2), wc3(3)] = Image2WorldMrFtk(header, c3(1), c3(2), c3(3));

c4 = [0 row-1 0];
wc4 = c4;
[wc4(1), wc4(2), wc4(3)] = Image2WorldMrFtk(header, c4(1), c4(2), c4(3));
