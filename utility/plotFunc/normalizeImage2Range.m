function output = normalizeImage2Range(input, lowR, highR)

minI = min(input(:));
maxI = max(input(:));

output = input;
output = lowR + (highR-lowR) .* (input-minI)./(maxI - minI);
