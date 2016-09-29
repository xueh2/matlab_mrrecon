
function perf_upsampled = PerformQuantitativePerfusion_OneImageSeries_Upsample(perf, dstAcqTimes, perf_acqTime)
% upsample AIF time stamps to regular grid

interp_method = 'linear';

SLC = size(perf, 4);

REP_perf = size(perf, 3);

AIF_time_stamp = dstAcqTimes;
N = numel(AIF_time_stamp);

% enforce the border-value boundary condition

RO = size(perf, 1);
E1 = size(perf, 2);

perf_upsampled = zeros(RO, E1, N, SLC);

for s=1:SLC
       
    perfAcq = perf_acqTime(:,s);
    perfAcq = [AIF_time_stamp(1)-0.01; perfAcq; AIF_time_stamp(end)+0.01;];
    
    perf_slc = zeros(RO, E1, REP_perf+2);
    perf_slc(:,:,1) = perf(:,:,1,s);
    perf_slc(:,:,2:end-1) = perf(:,:,:,s);
    perf_slc(:,:,end) = perf(:,:,end, s);
    
    perf_upsampled_slc = interp1(perfAcq, permute(perf_slc, [3 1 2]), AIF_time_stamp, interp_method);
    perf_upsampled_slc = permute(perf_upsampled_slc, [2 3 1]);
    perf_upsampled(:,:,:,s) = perf_upsampled_slc;
end

ind = find(isnan(perf_upsampled(:))==1);
perf_upsampled(ind(:)) = 0;
