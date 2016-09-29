% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load input data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load phantom_cine_PAT2_inplace_32ch
% load volunteer_sax_cine_PAT2_inplace_12ch
% load volunteer_sax_realtime_PAT4_tsense_12ch
load phantom_cine_PAT2_inplace_32ch_20081121; % 20081121_phantom_StudyLoid_4\meas_MID246_CV_fisp_cine_pat2_FID693


% contains Matlab variables:
%   headers           1x1               1018104  struct              
%   idir              1x49                   98  char                
%   output            1x1             305464888  struct              
%   protocol          1x1                 74144  struct              

ConfigHeader=DatHeader2struct(headers.Config);


NFirstRefLin = ConfigHeader.NFirstRefLin;
NRefLin = ConfigHeader.NRefLin;
if isfield(ConfigHeader,'AccelFactPE')
    AccelFactPE = ConfigHeader.AccelFactPE;
else
    AccelFactPE = protocol.sPat.lAccelFactPE;
end
if isfield(ConfigHeader,'nLinMeasOrig')
    nLinMeasOrig = ConfigHeader.nLinMeasOrig;
else
    nLinMeasOrig = ConfigHeader.NoOfFourierLines;
end

data=squeeze(output.data);
for i=size(data,1)+1:nLinMeasOrig
    data(i,:,:,:)=0; % fill out data matrix with zeros that were not acquired
end
acs=data(NFirstRefLin+1:NFirstRefLin+NRefLin,:,:,:); % (LIN, COL, CHA, PHA)
tmp=data;
data=zeros(size(tmp),'single');
data(1:AccelFactPE:end,:,:,:)=tmp(1:AccelFactPE:end,:,:,:); % (LIN, COL, CHA, PHA)
noise=output.noise;    % (LIN, COL, CHA)
clear tmp output

nLINmeas = ConfigHeader.NLinMeas;
nCOLmeas = ConfigHeader.NColMeas;

nRoFTLen = ConfigHeader.NRoFtlen;
nPeFTLen = ConfigHeader.NPeftLen;
% nPaFTLen = ConfigHeader.NPaFTLen;

nImageCols = ConfigHeader.NImageCols; % nImaLin
nImageLins = ConfigHeader.NImageLins; % nImaCol
% nImagePars = ConfigHeader.NImagePars; % nImaParts

rx_dwelltime_data = protocol.sRXSPEC.alDwellTime{1}; % from MeasYaps ascii header
% % input data
% data     k-space data
% acs      auto-calibration signal (k-space) data
% noise    noise data
%     
% % input parameters
% AccelFactPE          acceleration rate
% rx_dwelltime_noise   noise measurement dwell time
% rx_dwelltime_data    signal measurement dwell time
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recon configuration parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.recon.PAT.AccelFactPE = AccelFactPE;
params.recon.PAT.PAT_algorithm = 'sense';
params.recon.PAT.SENSE.regularization_parameter= 1.0; %0.5;


params.recon.raw_filter.data.type = 'gaussian';
params.recon.raw_filter.data.parameter = 1;
params.recon.raw_filter.acs.type  = 'hamming';
params.recon.raw_filter.acs.parameter = [];

params.recon.b1map.spatial_smoothing = 1;
params.recon.b1map.coilmap_smoothing_blocksize = 5;

params.recon.PAT_subtract_temporal_average_flag = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise pre-whitening (data and acs)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_dwelltime_data = protocol.sRXSPEC.alDwellTime{1};
noise=permute(noise, [3 2 1]); % set dimensions to (CHA, COL, LIN)
data=permute(data, [3 1 2 4]); % set dimensions to (CHA, LIN, COL, PHA)
acs=permute(acs, [3 1 2 4]); % set dimensions to (CHA, LIN, COL, PHA)
noisePrewhiteningMatrix = calculateNoisePrewhitener(noise, rx_dwelltime_data);
data = applyNoisePrewhitener(data, noisePrewhiteningMatrix);
acs = applyNoisePrewhitener(acs, noisePrewhiteningMatrix);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw filtering
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RawFilterCoeffsACS = calculateRawFilter (NRefLin,nCOLmeas, params.recon.raw_filter.acs);
RawFilterCoeffsData = calculateRawFilter (nLINmeas,nCOLmeas, params.recon.raw_filter.data);
acs_filtered = applyRawFilter (acs,RawFilterCoeffsACS);
data_filtered = applyRawFilter (data,RawFilterCoeffsData);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=fftshift(fft(data_filtered,nRoFTLen,3),3); % Readout FFT
tmp=tmp(:,:,end/4+1:3*end/4,:); % extract central FOV in RO dimension (for 2x RO oversampling)
tmp=fftshift(fft(tmp,nPeFTLen,2),2); % phase encode FFT
fftscale = 1/sqrt(nLINmeas*nCOLmeas/AccelFactPE); % scale FFT output for unity noise gain
imagedata_aliased = tmp*fftscale; % scale FFT output for unity noise gain

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT acs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=fftshift(fft(acs_filtered,nRoFTLen,3),3); % Readout FFT
tmp=tmp(:,:,end/4+1:3*end/4,:); % extract central FOV in RO dimension (for 2x RO oversampling)
tmp=fftshift(fft(tmp,nPeFTLen,2),2); % phase encode FFT
% tmp=tmp*(1/sqrt(NRefLin*nCOLmeas)); % scale FFT output for unity noise gain
imagedata_acs = tmp*fftscale; % scale FFT output with same scale factor as alias image data
% imagedata_acs = tmp;


% permute order of dimensions to LIN, COL, CHA, PHA
imagedata_aliased = permute(imagedata_aliased, [ 2 3 1 4]);
imagedata_acs = permute(imagedata_acs, [ 2 3 1 4]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate SENSE unmixing coeffs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHA=1; % calculate for 1st phase
[b1map]=calculateB1map(imagedata_acs(:,:,:,PHA), params.recon);
% mag=abs(sum(conj(b1map).*imagedata_acs(:,:,:,PHA),3));
mag=abs(sum(conj(b1map).*mean(imagedata_acs,4),3));
% mag = rss(imagedata_acs(:,:,:,PHA),3);
% b1map = imagedata_acs(:,:,:,PHA)./repmat(mag,[1 1 size(imagedata_acs,3)]);
[UnmixingCoeffs] = calculateSENSEUnmixingCoeffs(b1map, mag, params.recon);


imagedata_sense = zeros(size(imagedata_aliased,1),size(imagedata_aliased,2),size(imagedata_aliased,4),'single'); % initialize array
for PHA = 1:size(imagedata_aliased,4); % loop over phases
    disp(['phase: ',num2str(PHA),' out of ',num2str(size(imagedata_aliased,4))])
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply SENSE unmixing coeffs
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imagedata_sense(:,:,PHA) = sum(UnmixingCoeffs.*imagedata_aliased(:,:,:,PHA),3);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute g-factor map
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gmap=rss(UnmixingCoeffs,3);
figure; imagescn(gmap,[1 1.5]); colormap('jet')

figure; imagescn(abs(imagedata_sense),[],[],[],3)


