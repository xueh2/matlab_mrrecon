function out = mirdwt_helper(wlh, h, L)
scalefactor = 1;
n1 = size(wlh, 1);
n2 = size(wlh, 2);
%n = min(n1, n2);
n=n2/((L-1)*3+1);  % change here

nt = size(wlh, 3);
nc = size(wlh, 4);
for ct = 1:nt
for cc = 1:nc

% FIXME: only square images
wl = wlh(:, 1:n, ct, cc);
wh = wlh(:, (n + 1):max(n1, n2), ct, cc);

% normalization
wl = wl*scalefactor^((L-1));
for ll = 1:L-1
    wh(:,(ll-1)*n*3+1:ll*n*3) = ...
    wh(:,(ll-1)*n*3+1:ll*n*3)*scalefactor^(ll);
end

out(:,:,ct,cc) = complex( ...
 mirdwt(real(wl), real(wh), h, L - 1) ...
,mirdwt(imag(wl), imag(wh), h, L - 1) ...
);

end
end
