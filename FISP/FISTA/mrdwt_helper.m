function out = mrdwt_helper(in, h, L)

scalefactor = 1;

nt = size(in, 3);
nc = size(in, 4);
for ct = 1:nt
imt = [];
for cc = 1:nc

inim = in(:,:,ct,cc);
[wlr whr nlr] = mrdwt(real(inim), h, L - 1);
[wli whi nli] = mrdwt(imag(inim), h, L - 1);
wl = complex(wlr, wli);
wh = complex(whr, whi);
wl = wl * scalefactor^(-(L-1));
% n1 = size(wh, 1);
% n2 = size(wh, 2);
% n = min(n1, n2);
n=size(wl,2);  % change here

for ll = 1:L-1
     wh(:,(ll-1)*n*3+1:ll*n*3) = ...
     wh(:,(ll-1)*n*3+1:ll*n*3)*scalefactor^(-ll);
end
wlh = [wl wh];
out(:,:,ct,cc) = wlh;
end


end

