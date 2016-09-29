
function perf_Gd = PerformQuantitativePerfusion_GetGdFromLUT(perf_SRNorm, Gd, LUT)
% perf_Gd = PerformQuantitativePerfusion_GetGdFromLUT(perf_SRNorm, Gd, LUT)
% compute Gd from SR/PD

interp_method = 'linear';

minLUT = LUT(1);
maxLUT = LUT(end);

ind_min = find(perf_SRNorm(:)<=minLUT);
ind_max = find(perf_SRNorm(:)>maxLUT);

perf_SRNorm_Cut = perf_SRNorm;
perf_SRNorm_Cut(ind_max) = maxLUT;
% perf_SRNorm_Cut(ind_min) = 2*minLUT - perf_SRNorm(ind_min);
perf_SRNorm_Cut(ind_min) = minLUT;
   
perf_Gd = interp1(LUT, Gd, perf_SRNorm_Cut, interp_method);
% perf_Gd(ind_min) = -perf_Gd(ind_min);

ind = find(isnan(perf_Gd(:))==1);
perf_Gd(ind(:)) = 0;
