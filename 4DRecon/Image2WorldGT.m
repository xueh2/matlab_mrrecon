function [wx, wy, wz ] = Image2WorldGT(position, read_dir, phase_dir, slice_dir, pixel_spacing_ro, pixel_spacing_e1, pixel_spacing_z, RO, E1, x, y, z)
% image coordinates to world 
% [wx, wy, wz ] = Image2WorldGT(position, read_dir, phase_dir, slice_dir, pixel_spacing_ro, pixel_spacing_e1, pixel_spacing_z, x, y, z)

c = ceil([RO/2, E1/2]);

w = (x-c(1)) * pixel_spacing_ro * read_dir + (y-c(2)) * pixel_spacing_e1 * phase_dir + z * pixel_spacing_z * slice_dir + position;
wx = w(1);
wy = w(2);
wz = w(3);