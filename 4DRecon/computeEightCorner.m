function [wc1, wc2, wc3, wc4, wc5, wc6, wc7, wc8] = computeEightCorner(data, header)
% compute the world coordinates of eight corners for a 3D volume
% first slice (z=0)
% wc1 ---- wc2
% |        | 
% |        |
% |        |
% wc4 ---- wc3
%
% last slice (z=N-1)
% wc5 ---- wc6
% |        | 
% |        |
% |        |
% wc8 ---- wc7


row = size(data, 1);
col = size(data, 2);
depth = size(data, 3);

% compute the four corner of first slice
[wc1(1), wc1(2), wc1(3)] = Image2WorldMrFtk(header, 0, 0, 0);

[wc2(1), wc2(2), wc2(3)] = Image2WorldMrFtk(header, col-1, 0, 0);

[wc3(1), wc3(2), wc3(3)] = Image2WorldMrFtk(header, col-1, row-1, 0);

[wc4(1), wc4(2), wc4(3)] = Image2WorldMrFtk(header, 0, row-1, 0);

% compute the four corner of last slice
[wc5(1), wc5(2), wc5(3)] = Image2WorldMrFtk(header, 0, 0, depth-1);

[wc6(1), wc6(2), wc6(3)] = Image2WorldMrFtk(header, col-1, 0, depth-1);

[wc7(1), wc7(2), wc7(3)] = Image2WorldMrFtk(header, col-1, row-1, depth-1);

[wc8(1), wc8(2), wc8(3)] = Image2WorldMrFtk(header, 0, row-1, depth-1);
