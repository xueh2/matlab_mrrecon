function [params]=calculate_window_level(mag1,mag2,f1dnorm,f1d,Iref,params);

% function to automatically calculate values for image window-level and window-width

% notes
%   1) the 97 percentile calculation may be done more efficiently
%   2) the cleaning of stray pixels in binary masks using binary
%   morphological operations may not be necessary and/or may be done better

switch params.autowindowlevel
    case 'off'
        params.windowlevel.mag1_wl=0;
        params.windowlevel.mag1_ww=0;
        params.windowlevel.mag2_wl=0;
        params.windowlevel.mag2_ww=0;
		params.windowlevel.f1dnormreal_ww=0;
		params.windowlevel.f1dnormreal_wl=0;
		params.windowlevel.f1dnormreal_ww=0;
		params.windowlevel.f1dnormreal_wl=0;
        params.windowlevel.f1dreal_ww=0;
        params.windowlevel.f1dreal_wl=0;
    case 'on'
        % threshold is applied to mag2
        % mag2 is normalized such  that the noise has std dev = 1
        % therefore threshold_value is in units of std dev
        if params.ncoils==1
            threshold_value=2;
        else
			threshold_value=5;
        end		
        % Iref is median smoothed version of mag2
		% Iref=medfilt2(mag2,[Ns Ns],'symmetric'); % Ns=7 computed in
		% intensity_normalize.m
		mask=Iref>threshold_value;
        params.windowlevel.Iref=Iref;
        params.windowlevel.mask=mask;
        
		% mask represents regions of blood, muscle, fat, all above
		% threshold
    
        % compute mask2, a cleaned up version of mask using
        % binary morphological functions to delete stray pixels
        % by eroding and dilating
		% mask2=bwmorph(mask,'clean');
        
% 		mask2=bwmorph(mask,'erode',2);
% 		mask2=bwmorph(mask2,'dilate',2);
        mask2=mask;
        
        % tmp is a vector of values of f1dnorm for pixels in mask =1
		tmp=f1dnorm(mask);
        params.windowlevel.masked_values=real(tmp);
        % real(tmp)  are the values of phase sensitive IR image for pixels
        % which contain any tissue T1 (blood, fat, muscle, water)
        % Since f1norm is surface coil intensity corrected, the
        % distribution of values of real(tmp) takes on a series of distinct
        % levels for each value of tissue T1. typically, the blood and fat
        % are high values, and the muscle and any possible water (effusion)
        % are distinctly lower values. The median generally falls between
        % these distributions, therefore, the binary mask created by
        % f1dnorm<median(real(tmp(:))) is generally muscle and possible
        % water due to any effusion, and, of course, noise in regions with
        % no tissue. the value threshold, combined this with mask2 to
        % remove the noise. threshold is therefore mostly muscle region.
% 		threshold=mask2.*(f1dnorm<median(real(tmp(:))));
        threshold=mask2.*(real(f1dnorm)<median(real(tmp(:)))); % changed 7/2/03
        % threshold2 is slightly cleaned up
% 		threshold2=bwmorph(threshold,'erode',1);
        threshold2=threshold;
        params.windowlevel.threshold=threshold;
        
		% 		threshold2=bwmorph(threshold2,'dilate',2);
		% normal_level is the level of muscle (or normal myocardium)
        normal_level=mean(real(f1dnorm(threshold2==1)));
		params.windowlevel.normal_level=normal_level;
        params.windowlevel.normal_values=real(f1dnorm(threshold2==1));
        % the window level will be determined as follows:
        % find the level which corresponds to the 97 percentile in regions
        % which have tissue (any type) using mask2. This percentile can be
        % adjusted. Use this value as the upper value displayed. Calculate
        % the window as 10% (this can be adjusted) greater than the range
        % calculated as (upperlevel-normal_level), i.e normal level of
        % myocardium should be approximately nulled, but increasing the
        % range slightly improves appearance.
        % next find the upperlevel using a histogram function (pdf) and
        % then a cumulative sum to calculate the cumulative distribution
        % (cdf). There is probably a much simpler way to compute the
        % percentile statistic.
        tmp=real(f1dnorm(mask2));
        if ~isempty(tmp)
			[pdf,val]=hist(tmp,64);
			cdf=cumsum(pdf);
			cdf=cdf/cdf(end);
            params.windowlevel.histogram=pdf;
            params.windowlevel.histogram_bins=val;
			f1dnormreal_val97=val(find(cdf>=0.97));
			f1dnormreal_val97=f1dnormreal_val97(1);
            params.windowlevel.f1dnormreal_val97=f1dnormreal_val97;
%             percentile = prctile(tmp,97);
%             disp(['f1dnormreal_val97: ',num2str(f1dnormreal_val97)])
%             disp(['percentile: ',num2str(percentile)])
%             
			range=1.1*f1dnormreal_val97-normal_level;
			min_level=f1dnormreal_val97-range;
            params.windowlevel.normal_level=normal_level;
            params.windowlevel.range=range;
            params.windowlevel.min_level=min_level;
            
            params.windowlevel.f1dnormreal_wl=min_level+range/2;
			params.windowlevel.f1dnormreal_ww=range;
        else
            params.windowlevel.f1dnormreal_wl=0;
            params.windowlevel.f1dnormreal_ww=0
        end
        
        tmp=real(f1d(mask2));
        if ~isempty(tmp)
			[pdf,val]=hist(tmp,64);
			cdf=cumsum(pdf);
			cdf=cdf/cdf(end);
			% f1dreal_val05=val(find(cdf<=0.05));
			% f1dreal_val05=f1dreal_val05(end);
			f1dreal_val97=val(find(cdf>=0.97));
			f1dreal_val97=f1dreal_val97(1);
			normal_level=mean(real(f1d(threshold2==1)));
			range=1.1*f1dreal_val97-normal_level;
			min_level=f1dreal_val97-range;
			params.windowlevel.f1dreal_wl=min_level+range/2;
			params.windowlevel.f1dreal_ww=range;
        else
            params.windowlevel.f1dreal_wl=0;
            params.windowlevel.f1dreal_ww=0
        end   
		
		[pdf,val]=hist(mag1(:),64);
		cdf=cumsum(pdf);
		cdf=cdf/cdf(end);
		mag1_val97=val(find(cdf>=0.97));
		mag1_val97=mag1_val97(1);        
		params.windowlevel.mag1_wl=(mag1_val97+0)/2;
		params.windowlevel.mag1_ww=(mag1_val97-0);
		
		[pdf,val]=hist(mag2(:),64);
		cdf=cumsum(pdf);
		cdf=cdf/cdf(end);
		mag2_val97=val(find(cdf>=0.97));
		mag2_val97=mag2_val97(1);
		params.windowlevel.mag2_wl=(mag2_val97+0)/2;
		params.windowlevel.mag2_ww=(mag2_val97-0);
end

return

%         const double dRangePSIRNorm    = 1.1 * dCutoffPSIRNorm-dMeanPSIRNorm;
%         const double dMinLevelPSIRNorm = dCutoffPSIRNorm - dRangePSIRNorm;



