function AIF_signal = perfusion_AIF_signal(AIF, mask, numPD, percentageAIF)
% get AIF signal

ind = find(mask(:)>0);

RO = size(AIF, 1);
E1 = size(AIF, 2);
N = size(AIF, 3);

AIF_signal = zeros(N, 1);

for n=1:N
    d = AIF(:,:,n);
    if(n<=numPD)
        AIF_signal(n) = mean(d(ind(:)));
    else
        sp = d(ind(:));
        if(percentageAIF>0)
            sp = sort(sp(:));
            m = numel(sp);
            sp = sp( floor(m*percentageAIF/100):end);
        end
        
        AIF_signal(n) = mean(sp(:));
    end
end