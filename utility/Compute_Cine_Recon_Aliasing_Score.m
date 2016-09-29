
function score = Compute_Cine_Recon_Aliasing_Score(Grappa_Recon, R)

a_0 = Grappa_Recon;
s = size(a_0);
if length(s)==3, s(4) = 1;end
if length(s)==2, s(3) = 1; s(4) = 1;end
mask = zeros(s(1), s(2), s(3));

score = zeros([s(3) s(4)]);

for j = 1:s(4)
    j = j;
    for i=1:s(3), %loop frames
        temp = a_0(:,:,i,j); % temp(:,1:50) = 0;
        I_0 = sort(temp(:), 'descend');
        
        % mask(:,:,i) = temp > I_0(round(s(1)*s(2)*0.05)) ; % cine
        mask(:,:,i) = temp > I_0(round(s(1)*s(2)*0.25)) ; % perfusion
        
        % mask(:,:,i) = temp > 10;
        % mask(:,:,i) = temp > 2*min(Grappa_Recon(:));
    end
    score(1:s(3), j) = Aliasing_Score_Mask(mask.*Grappa_Recon(:, :, :, j), Grappa_Recon(:, :, :, j), R);
end






























