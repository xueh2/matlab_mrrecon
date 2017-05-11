function res = PerformSpatioTemporalFiltering(a, thres_temporal, thres_spatial, wave_name, level)
% perform the joint spatio-temporal filtering
% a: [RO E1 N], N is the temporal dimension
% wave_name : db1, db2, db3, db4, db5
% level: number of wavelet transform levels
% thres_temporal: ratio of energy loss allowed along the temporal dimension, compared to the total mean, default 0.001
% thres_spatial: ratio of spatial smoothing allowed, default 1.0
% res = PerformSpatioTemporalFiltering(input, thres_temporal, thres_spatial, wave_name, level)

RO = size(a, 1);
E1 = size(a, 2);
N = size(a, 3);
s = size(a);

% perform KL transform
% [eigenImages, V, D] = KL_Eigenimage(a);

b = (reshape(a, [s(1)*s(2),s(3)] ));
t_mean = mean(b,1);
A = (b'*b - s(1)*s(2)*t_mean'*t_mean)/(s(1)*s(2)-1);
[V,D] = eig((A + A')/2);
eigenImages = b*V;
eigenImages = reshape(eigenImages, s);

D = diag(D);

for ii=1:N-1
    e = sum(D(1:ii));
    % if(e/D(end)>=thres_temporal)
    if(sqrt(e/sum(D))>=thres_temporal)
        break;
    end
end

for ii2=N-1:-1:1
    if(sqrt(D(ii2)/sum(D))<=thres_temporal)
        break;
    end
end

ii2 = N-1 - ii2;

ii = floor((ii+ii2)/2);

if(ii>3*N/4)
    ii = floor(N/2);
end

disp(['discard ' num2str(ii) ' modes along the temporal dimension ... ']);
sigma2_noise = mean(D(1:ii)); % noise level

% perform wavelet transform
wav = Matlab_gt_redundant_wavelet_transform(eigenImages(:,:,ii+1:N), wave_name, 2, 1, level);

% perform filtering on the wav coeff
wN = size(wav, 4);
wW = size(wav, 3);
sigma2_signal = zeros(wW, wN);
thres_signal = zeros(wW, wN);
for n=1:wN
    for p=1:wW % for eveyr high freq wavelet coefficient
        
        pp = wav(:,:,p,n);
        sigma2_signal(p, n) = var(pp(:),1);
        t = sigma2_signal(p, n) - sigma2_noise;
        if(t<0) 
            t = sigma2_signal(p, n); 
        end
        
        thres_signal(p, n) = thres_spatial * sigma2_noise / (sqrt(t)+eps);
    end
end

wav_f = wav;

% perform the soft thresholding
for n=1:wN-1
    for p=1:wW % for eveyr high freq wavelet coefficient
        t = thres_signal(p, n);
        w = wav(:,:,p,n);
        wav_f(:,:,p,n) = softThresholding(w, t);
    end
end

% perform inverse wavelet transform
eigenImages2 = Matlab_gt_redundant_wavelet_transform(wav_f, wave_name, 2, -1, level);
eigenImages(:,:,ii+1:N) = eigenImages2;

% inverse KL transform
V(:,1:ii) = 0;
V_in = V';

tt = reshape(double(eigenImages), [RO*E1 N]);
res = (reshape(tt*V_in, [RO,E1,N] ) ) ;

