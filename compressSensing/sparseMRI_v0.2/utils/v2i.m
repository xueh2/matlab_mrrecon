function im = v2i(z, ind, N)
% vector to image for points in ind
im = zeros(N);
im(ind) = z;