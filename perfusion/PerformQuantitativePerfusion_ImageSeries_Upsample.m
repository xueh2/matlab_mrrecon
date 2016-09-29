
function [AIF_upsampled, perf_upsampled, deltaT] = PerformQuantitativePerfusion_ImageSeries_Upsample(AIF, perf, AIF_acqTime, perf_acqTime, upsamplingRatio, deltaT)
% upsample AIF time stamps to regular grid

interp_method = 'linear';

REP = size(AIF, 3);
SLC = size(perf, 4);

REP_perf = size(perf, 3);

t_start = AIF_acqTime(1);
t_end = AIF_acqTime(end);

if(deltaT==0)
    N = floor(REP*upsamplingRatio);
    deltaT = (t_end - t_start) / (N-1);
end

AIF_time_stamp = t_start:deltaT:t_end;
N = numel(AIF_time_stamp);

AIF_upsampled = interp1(AIF_acqTime, permute(AIF, [3 1 2]), AIF_time_stamp, interp_method);
AIF_upsampled = permute(AIF_upsampled, [2 3 1]);

% enforce the border-value boundary condition

RO = size(perf, 1);
E1 = size(perf, 2);

RO_AIF = size(AIF, 1);
E1_AIF = size(AIF, 2);

perf_upsampled = zeros(RO, E1, N, SLC);

for s=1:SLC
       
    perfAcq = perf_acqTime(:,s);
    perfAcq = [AIF_time_stamp(1)-0.01; perfAcq];
    
    perf_slc = zeros(RO, E1, REP_perf+1);
    perf_slc(:,:,1) = perf(:,:,1,s);
    perf_slc(:,:,2:end) = perf(:,:,:,s);
    
    perf_upsampled_slc = interp1(perfAcq, permute(perf_slc, [3 1 2]), AIF_time_stamp, interp_method);
    perf_upsampled_slc = permute(perf_upsampled_slc, [2 3 1]);
    perf_upsampled(:,:,:,s) = perf_upsampled_slc;
end

ind = find(isnan(perf_upsampled(:))==1);
perf_upsampled(ind(:)) = 0;
