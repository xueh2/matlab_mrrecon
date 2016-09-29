function [thres, filteredWaveCoeff, filteredImg] = estimateBayesWaveletThreshold(im, N, wname, ratioForThres)
% Ref to paper "Adaptive wavelet thresholding for image denoising and compression"
% given a input 2D image im, estimate the wavelet thresholds
% N : number of wavelet transform levels
% wname : wavelet name, following the wavelet toolbox
% ratioForThres : ratio to scale the estiamted thresholds, realThresUsed = thres*ratioForThres

% first, perform the non-decimated wavelet transform

if ( isreal(im) )
    % perform estimation for real image
    WT_R = ndwt2(im, N, wname);
    
    % estimate noise variance from HH subband
    HH = WT_R.dec{end}; % always use the HH for the level 1
    
    sigmaNoise = median(abs(HH(:)))/0.6745;
    
%     options.alpha = 0.99;
%     options.lts = 1;
%     [res,raw] = fastmcd(abs(HH(:)),options);
%     sigmaNoise = res.center;
    
    numOfLevels = numel(WT_R.dec);
    thres = zeros(numOfLevels, 1);
    
    % estimate signal variance
    for n=1:numOfLevels        
        coeff = WT_R.dec{n}; % imshow(coeff, [])
       
        sigmaYSquare = coeff(:)'*coeff(:) / numel(coeff);
        
        % sigmaYSquare = std(coeff(:))^2;
        sigmaX = sqrt(max(0, sigmaYSquare-sigmaNoise*sigmaNoise));
        thres(n) = sigmaNoise*sigmaNoise/sigmaX;
    end
    
    % perform soft-thresholding
    filteredWaveCoeff = WT_R.dec;
    for n=1:numOfLevels
        if ( thres(n) == Inf )
            filteredWaveCoeff{n} = zeros(size(filteredWaveCoeff{n}));
        else
            coeff = WT_R.dec{n};
            filteredWaveCoeff{n} = wthresh(coeff, 's', thres(n)*ratioForThres);
            
%             res = abs(coeff)-thres(n)*ratioForThres;
%             res = (res + abs(res))/2;
%             filteredWaveCoeff{n} = sign(coeff).*res;    
        end
    end
    
    % get the filtered image
    WT_R.dec = filteredWaveCoeff;
    filteredImg = indwt2(WT_R);
else
    % perform estimation for complex image
    WT_R = ndwt2(real(im), N, wname);
    WT_I = ndwt2(imag(im), N, wname);
       
    % estimate noise variance from HH subband of level 1
    HH = WT_R.dec{end} + i* WT_I.dec{end};
    % sigmaNoise = median(abs(HH(:)))/0.6745;
    % sigmaNoise = median(HH(:))/0.6745; <-- oversmoothing
    
    % sigmaNoise = mean(abs(HH(:))); <-- looks ok
    
    options.alpha = 0.99;
    options.lts = 1;
    [res,raw] = fastmcd(abs(HH(:)),options);
    sigmaNoise = res.center;
    
    numOfLevels = numel(WT_R.dec);
    thres = zeros(numOfLevels, 1);
    
    % estimate signal variance
    for n=1:numOfLevels        
        coeff = WT_R.dec{n} + i*WT_I.dec{n};
        
        sigmaYSquare = coeff(:)'*coeff(:) / numel(coeff);
        % sigmaYSquare = std(coeff(:));
        sigmaX = sqrt(max(0, sigmaYSquare-sigmaNoise*conj(sigmaNoise)));
        % sigmaX = sigmaYSquare;
        thres(n) = sigmaNoise*conj(sigmaNoise)/sigmaX;
    end
    
    % perform soft-thresholding
    filteredWaveCoeff = WT_R.dec;
    filteredWaveCoeff_R = WT_R.dec;
    filteredWaveCoeff_I = WT_I.dec;
    for n=1:numOfLevels
        if ( thres(n) == Inf )
            filteredWaveCoeff_R{n} = zeros(size(filteredWaveCoeff_R{n}));
            filteredWaveCoeff_I{n} = zeros(size(filteredWaveCoeff_I{n}));
        else
            coeff = WT_R.dec{n} + i*WT_I.dec{n}; 
            y = wthresh(coeff, 's', thres(n)*ratioForThres);
            filteredWaveCoeff_R{n} = real(y);
            filteredWaveCoeff_I{n} = imag(y);
            filteredWaveCoeff{n} = y;

%             res = abs(coeff)-thres(n)*ratioForThres;
%             res = (res + abs(res))/2;
%             unity = coeff./(abs(coeff)+eps);
%             filteredWaveCoeff_R{n} = real(unity.*res);
%             filteredWaveCoeff_I{n} = imag(unity.*res);
%             filteredWaveCoeff{n} = filteredWaveCoeff_R{n} + i * filteredWaveCoeff_I{n};
        end
    end
    
%     % do not filtering the LL
%     filteredWaveCoeff_R{1} = WT_R.dec{1};
%     filteredWaveCoeff_I{1} = WT_I.dec{1};
%     filteredWaveCoeff{1} = filteredWaveCoeff_R{1} + i*filteredWaveCoeff_I{1}; 
    
    % get the filtered image
    WT_R.dec = filteredWaveCoeff_R;
    WT_I.dec = filteredWaveCoeff_I;
    filteredImg = indwt2(WT_R) + i*indwt2(WT_I);

% wname='db5';
% 
% [swa,swh,swv,swd]=swt2(im,3,wname);
% tmp=swd(:,:,1);
% Nvar=median(abs(tmp(:)))/0.6745;
% sorh='s';
% 
% % For subband - LL
%     for j=1:3
%         beta1=sqrt(log(length(swa(:,:,j))/3));
%         Ssig1=std(swa(:));
%         thr1=(beta1*Nvar^2)/Ssig1;
%         swa(:,:,j)=wthresh(swa(:,:,j),sorh,thr1);
%     end    
% % For subband - Horizontal Details
%     for j=1:3
%         beta1=sqrt(log(length(swh(:,:,j))/3));
%         Ssig1=std(swh(:));
%         thr1=(beta1*Nvar^2)/Ssig1;
%         swh(:,:,j)=wthresh(swh(:,:,j),sorh,thr1);
%     end
% % For subband - Vertical Details
%     for k=1:3
%         beta2=sqrt(log(length(swv(:,:,k))/3));
%         Ssig2=std(swv(:));
%         thr2=(beta2*Nvar^2)/Ssig2;
%         swv(:,:,k)=wthresh(swv(:,:,k),sorh,thr2);
%     end
% % For subband - Diagonal Details
%     for l=1:3
%         beta3=sqrt(log(length(swd(:,:,l))/3));
%         Ssig3=std(swd(:));
%         thr3=(beta3*Nvar^2)/Ssig3;
%         swd(:,:,l)=wthresh(swd(:,:,l),sorh,thr3);
%     end
%     clean_eigen1=iswt2(swa,swh,swv,swd,wname);
%     filteredImg = clean_eigen1;
%     thres = [];
%     filteredWaveCoeff = [];
end