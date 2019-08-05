function [new_data, new_cine, new_mask, data] = fixCentroid(data, cine, mask, add_centroid)
imsz = size(data);
PHS = imsz(3);
centroidE = zeros(PHS, 1);
centroidR = zeros(PHS, 1);
for i = 1:(PHS)
    [Rmesh, Emesh] = meshgrid(1:imsz(2), 1:imsz(1));
    % find centroid of mask
    centroidR(i) = sum(sum(Rmesh(mask(:, :, i) == 1)))/sum(sum(mask(:, :, i) == 1));
    centroidE(i) = sum(sum(Emesh(mask(:, :, i) == 1)))/sum(sum(mask(:, :, i) == 1));
    centroid = [centroidE(i), centroidR(i)];
end

dE = centroidE - centroidE(10);
dR = centroidR - centroidR(10);
maxE = floor(imsz(1) - (centroidE(10) + max(dE)))-1;
minE = floor(centroidE(10) + min(dE))-1;
maxR = floor(imsz(2) - (centroidR(10) + max(dR)))-1;
minR = floor(centroidR(10) + min(dR))-1;


new_sz = [maxE+minE, maxR+minR, PHS];

new_data = zeros(new_sz);
new_mask = zeros(new_sz);
new_cine = zeros(new_sz);
centroidE = round(centroidE);
centroidR = round(centroidR);

for i = 1:PHS
    [centroidE(i)-minE,maxE+centroidE(i), centroidR(i)-minR,maxR+centroidR(i)]
    new_data(:, :, i) = data(centroidE(i)-minE:maxE+centroidE(i)-1, centroidR(i)-minR:maxR+centroidR(i)-1, i);
    new_mask(:, :, i) = mask(centroidE(i)-minE:maxE+centroidE(i)-1, centroidR(i)-minR:maxR+centroidR(i)-1, i);
    new_cine(:, :, i) = cine(centroidE(i)-minE:maxE+centroidE(i)-1, centroidR(i)-minR:maxR+centroidR(i)-1, i);
end

if add_centroid
    for i = 1:PHS
        data(centroidE(i), centroidR(i), i) = 2000;
        new_data(minE, minR, i) = 2000;
    end
end

end