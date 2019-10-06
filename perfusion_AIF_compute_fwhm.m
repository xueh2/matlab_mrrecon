function FWHM = perfusion_AIF_compute_fwhm(AIF, foot, peak, valley)
% compute FWHM for aif signal
% FWHM = perfusion_AIF_compute_fwhm(AIF, foot, peak, valley)

N = size(AIF, 1);

dt = 0.01;

t = 0:dt:N-1;

AIF_sampled = interp1(0:N-1, AIF, t);

if(valley<foot)
    FWHM = 0;
    return;
end

if(valley<peak)
    valley = peak + (peak-foot);
    if(valley>N-1)
        valley = N-1;
    end
end

foot_sampled = ceil(foot/dt);
valley_sampled = floor(valley/dt);
peak_sampled = round(peak/dt);

maxA = max(AIF_sampled(foot_sampled:valley_sampled));

for f=foot_sampled:peak_sampled
    v = AIF_sampled(f+1);
    if(v>=0.5*maxA)
        break;
    end
end

for f2=peak_sampled:-1:foot_sampled
    v = AIF_sampled(f2+1);
    if(v<=0.5*maxA)
        break;
    end
end

fwhm_s = (f+f2)/2;

for f=peak_sampled:valley_sampled
    v = AIF_sampled(f+1);
    if(v<=0.5*maxA)
        break;
    end
end

for f2=valley_sampled:-1:peak_sampled
    v = AIF_sampled(f2+1);
    if(v>=0.5*maxA)
        break;
    end
end

fwhm_e = (f+f2)/2;

FWHM = (fwhm_e - fwhm_s) * dt;
if(FWHM<0)
    FWHM = (valley - foot)/3;
end