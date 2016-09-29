function  noisePrewhiteningMatrix = calculateNoisePrewhitener2(noise, rx_dwelltime_data, rx_dwelltime_prescan_noise);
% function  noisePrewhiteningMatrix = calculateNoisePrewhitener(noise, rx_dwelltime_data, readout_length);
% 
% inputs:
%     noise is the prescan noise: CHA (channels), COL (readout), LIN (lines)
%     rx_dwelltime_data is the receiver dwell time (ADC sample period in nsec) used for acquisition of imaging data
% 
% output:
%     noisePrewhiteningMatrix is the complex noise prewhitening matrix (CHA x CHA)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

s=size(noise);
readout_length=s(2); % # readout samples including oversampling
noise = reshape(noise, [s(1) prod(s(2:end))]); % CHA x samples
M=size(noise,2); % number of noise samples

% calculate noise scale factor assuming prescan noise acquired using
% 130 Hz/pixel (using SBBNoiseAdj)
if ( rx_dwelltime_prescan_noise < 0 )
    rx_dwelltime_prescan_noise=10^9/readout_length/130; % ADC sample period in nsec for noise measurement
end

rcvr_noise_bandwidth =0.79; % noise equivalent bandwidth of digital receiver filter
noise_scale_factor=rx_dwelltime_prescan_noise/rx_dwelltime_data/rcvr_noise_bandwidth;

% calculate noise correlation matrix (Hermitian)
Rn=(1/(M-1))*(noise*noise');
% Rn=Rn*noise_scale_factor;
Rn=(Rn+Rn')/2; % guarantee positive definite matrix
% calculate Cholesky factorization
L=chol(Rn,'lower'); % L is lower triangular matrix such that Rn = L*L'
% calculate pre-whitenening matrix
noisePrewhiteningMatrix=inv(L)*sqrt(2)*sqrt(1/noise_scale_factor); % sqrt(2) used so that std(real) & std(imag) both equal 1

return






