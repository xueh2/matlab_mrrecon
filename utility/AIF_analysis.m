function AIF_analysis(idir, odir);
% function AIF_analysis(idir, odir);

% idir = 'T:\ReconResults\Barts_new\20191117\Perfusion_AIFR3_2E_Interleaved_42110_760681510_760681519_39_20191117-092558\';
% idir = 'T:\ReconResults\Barts_new\20190831\Perfusion_AIFR3_2E_Interleaved_42110_488936648_488936657_322_20190831-132944\';
% idir = 'Z:\ReconResults\Barts_new\20190831\Perfusion_AIFR3_2E_Interleaved_42110_488936486_488936495_43_20190831-083713\';


idir = [idir,filesep,'DebugOutput',filesep];
if nargin == 1; 
    odir = 'D:\work\projects\perfusion\AIF_analysis\';
end
if isempty(idir)
    return
end

idir = strrep(idir,'\\','\');
ind = strfind(idir,'\');
basefilename = idir(ind(4)+1:ind(5)-1);
h5filename = find_h5filename(basefilename);


% get aif
try
    aif1 = analyze75read([idir,'aif_moco_upsampled.img']);
    aif2 = analyze75read([idir,'aif_moco_second_echo_upsampled.img']);
    AcqTimes = analyze75read([idir,'Perf_AcqTimes_0.img']);
catch
    return
end

time = 1e-3*(AcqTimes-AcqTimes(1));
rr = diff(time);
RR = mean(rr);


% segment LV and RV using AI masks
filename =  [idir,filesep,'aif_mask_AI_probs.img'];
probs = analyze75read(filename);
%figure; imagescn(probs)
lv_mask = probs(:,:,2)>0.5;
rv_mask = probs(:,:,3)>0.5;
% check whether orientation needs permuting:
s1 = size(lv_mask);
s2 = size(aif1); s2 = s2(1:2);
% read dimension is always 128
ind = find (s2 == 128);
if isempty(ind) % old process before interp up to 128.
    ind = find(s2 ==64);
    if length(ind)>1;
        % then must be 64 x 64 -> do not permute
    elseif s1(ind) == 64;
        % then aif and mask both have same read dimension
    else
        % permute
        lv_mask = lv_mask';
        rv_mask = rv_mask';    
    end
else
    if length(ind)>1;
        % then must be 128 x 128 -> do not permute
    elseif s1(ind) == 128;
        % then aif and mask both have same read dimension
    else
        % permute
        lv_mask = lv_mask';
        rv_mask = rv_mask';    
    end
end
% pad mask to size of aif:
try
    lv_mask = zeropad(lv_mask, size(aif1,2), size(aif1,1), 'center');
    rv_mask = zeropad(rv_mask, size(aif1,2), size(aif1,1), 'center');
catch
    lv_mask = zpad(lv_mask, size(aif1,1), size(aif1,2));
    rv_mask = zpad(rv_mask, size(aif1,1), size(aif1,2));
end


% figure; imagescn(lv_mask + 2*rv_mask)
% rv_mask_index = find(rv_mask' >0);
rv_mask_index = find(rv_mask >0);
if isempty(rv_mask_index)
    return
end
            % aif = analyze75read([idir,'aif_moco.img']);
            % figure; imagescn(aif,[],[],[],3)
            % % apply RV mask
            % delete rv
            % for i = 1:size(aif,3)
            %     tmp = aif(:,:,i);
            %     rv(i) = mean(tmp(mask_index));
            % end
            % figure; plot(rv)

clear rv1 rv2
for i = 1:size(aif1,3)
    tmp = aif1(:,:,i);
    rv1(i) = mean(tmp(rv_mask_index));
    tmp = aif2(:,:,i);
    rv2(i) = mean(tmp(rv_mask_index));    
end
% figure; plot([rv1; rv2]')

% perform T2* correction
try
    pd = analyze75read([idir,'aifPD_for_TwoEcho_T2StartCorrection_0.img']);
catch
    return
end
addpath('Z:\Share\temp\gt_mex');
% ---------------------------------------------------------------------
% 8 Input paras:
% 	aif                                        : RO*E1*N, aif series, in double
% 	PD                                         : RO*E1, PD image, in double
% 	aif_second_echo                            : RO*E1*N, aif series, in double
% 	aif_LV_mask                                : RO*E1, aif LV mask, in double
% 	aif_signal_picking_percentage              : percentage to keep aif pixels (e.g. 0.15)
% 	TE0                                        : echo time for 1st echo
% 	TE1                                        : echo time for 2nd echo
% 	DebugFolder                                : if not empty, intermediate results are stored in this folder
% 10 Output para:
% 	final_mask                                 : AIF LV mask
% 	foot                                       : detected foot
% 	peak                                       : detected peak
% 	valley                                     : detected valley
% 	upslope                                    : detected upslope
% 	auc                                        : detected auc
% 	FWHM                                       : detected FWHM
% 	aif_echo0                                  : AIF signal for 1st echo
% 	aif_echo1                                  : AIF signal for 2nd echo
% 	aif_echo0_R2Star_Corrected                 : AIF signal for 1st echo, with T2* correction
% ==============================================================================================
aif = double(aif1);
aif_second_echo =  double(aif2);
aif_PD = double(pd(:,:,1));
aif_mask = double(rv_mask);
aif_signal_picking_percentage = 15; % on percent units
TE0 = 0.76;
TE1 = 1.76;
DebugFolder = [];
try
    [final_mask,rv_foot, rv_peak, rv_valley, upslope, auc, FWHM, aif_echo0, aif_echo1, aif_echo0_R2Star_Corrected] = Matlab_gt_perfusion_flow_aif_signal_from_mask(aif, aif_PD, aif_second_echo, aif_mask, aif_signal_picking_percentage, TE0, TE1, DebugFolder);
catch
    return
end
% figure; plot([aif_echo0,aif_echo1,aif_echo0_R2Star_Corrected])


% perform LUT correction conversion to [Gd]
aif_cin_LUT_Valid = analyze75read([idir,'aif_cin_LUT_Valid']);
aif_cin_Gd_Valid = analyze75read([idir,'aif_cin_Gd_Valid']);
% figure; plot(aif_cin_Gd_Valid, aif_cin_LUT_Valid);
clear lut
lut(:,1) = aif_cin_Gd_Valid';
lut(:,2) = aif_cin_LUT_Valid';
addpath('d:\work\simulation\perfusion_lut_revised');
in = aif_echo0_R2Star_Corrected;
rv_gd = apply_lut_correction(in, lut);
rv_gd = double(rv_gd);
% subtract baseline Gd
rv_gd = rv_gd - mean(rv_gd(1:10));

if max(rv_gd) < 1.5
    return
end

time = [0:.5:floor(length(rv_gd)/2)]'; % 0.5 sec steps after interpolation
time = time(1:length(rv_gd));

% hfig = figure; hold on; 
% plot(time, rv_gd)
% plot(time(round([rv_foot rv_foot])),[0 1.25*max(rv_gd)],'b:')
% plot(time(round([rv_valley rv_valley])),[0 1.25*max(rv_gd)],'b:')
% plot(time(round([rv_peak rv_peak])),[0 1.25*max(rv_gd)],'b:')
% axis([0 max(time) 0 1.25*max(rv_gd)]);
% box on

% Process LV
aif_mask = double(lv_mask);
try
    [final_mask,lv_foot, lv_peak, lv_valley, upslope, auc, FWHM, aif_echo0, aif_echo1, aif_echo0_R2Star_Corrected] = Matlab_gt_perfusion_flow_aif_signal_from_mask(aif, aif_PD, aif_second_echo, aif_mask, aif_signal_picking_percentage, TE0, TE1, DebugFolder);
catch
    return
end
% figure; plot([aif_echo0,aif_echo1,aif_echo0_R2Star_Corrected])
in = aif_echo0_R2Star_Corrected;
lv_gd = apply_lut_correction(in, lut);
lv_gd = double(lv_gd);
% subtract baseline Gd
lv_gd = lv_gd - mean(lv_gd(1:10));
lv_mask_size = sum(final_mask);

if max(lv_gd) < 1.5
    return
end


figure
hold on
plot(time, lv_gd)
plot(time(round([lv_foot lv_foot])),[0 1.25*max(rv_gd)],'r:')
plot(time(round([lv_valley lv_valley])),[0 1.25*max(rv_gd)],'r:')
plot(time(round([lv_peak lv_peak])),[0 1.25*max(rv_gd)],'r:')
axis([0 max(time) 0 1.25*max(rv_gd)]);
box on

% integrate concentration curves from foot to peak (without fitting)
rv_dose = 0.5 * sum(rv_gd(round(rv_foot):round(rv_valley)));
lv_dose = 0.5 * sum(lv_gd(round(lv_foot):round(lv_valley)));



% fit to log normal and integrate total Gd dose
addpath('d:\work\utilities')
% y = my_lognpdf(x,mu,sigma, x0, ampl)
func = @(fit,xdata)my_lognpdf(xdata,fit(1),fit(2),fit(3),fit(4));

time_interp = linspace(time(1),time(end),1000); % upsampled
try
    fit = lsqcurvefit(func,[1 1  time(round(rv_foot)) rv_dose], time(1:round(rv_valley)), rv_gd(1:round(rv_valley)));
    rv_gd_est = my_lognpdf(time_interp ,fit(1),fit(2),fit(3),fit(4));
    % figure; hold on; plot(time,rv_gd); plot(time_interp, rv_gd_est, 'r');
    fit = lsqcurvefit(func,[1 1  time(round(lv_foot)) lv_dose], time(1:round(lv_valley)), lv_gd(1:round(lv_valley)));
    lv_gd_est = my_lognpdf(time_interp ,fit(1),fit(2),fit(3),fit(4));
    % figure; hold on; plot(time,lv_gd); plot(time_interp, lv_gd_est, 'r');
    
    x = double(time(1:round(rv_valley)));
    y = double(rv_gd(1:round(rv_valley)));
    mus = [0.1:0.1:2.0];
    sigmas = [0.1:0.1:2.0];
    x0s = [time(round(rv_foot))-3:time(round(rv_foot))+5];
    ampls = [0.1*rv_dose:0.1*rv_dose:2.0*rv_dose];
    [best_params_rv, best_cost_rv, initial_params_rv, initial_cost_rv] = Matlab_gt_QPerf_PTT_fitting(x, y, mus, sigmas, x0s, ampls);
    
    x = double(time(1:round(lv_valley)));
    y = double(lv_gd(1:round(lv_valley)));
    mus = [0.1:0.1:2.0];
    sigmas = [0.1:0.1:2.0];
    x0s = [time(round(lv_foot))-3:time(round(lv_foot))+5];
    ampls = [0.1*lv_dose:0.1*lv_dose:2.0*lv_dose];
    [best_params_lv, best_cost_lv, initial_params_lv, initial_cost_lv] = Matlab_gt_QPerf_PTT_fitting(x, y, mus, sigmas, x0s, ampls);

catch
    return
end

utDir = 'D:\gtuser\mrprogs\gtprep\ut\data\PTT';
t_gd = single(time(:));
Matlab_gt_write_analyze(t_gd, CreateGtImageHeader(t_gd), fullfile(utDir, 't_gd'));

lv_gd = single(lv_gd(:));
Matlab_gt_write_analyze(lv_gd, CreateGtImageHeader(lv_gd), fullfile(utDir, 'lv_gd'));
rv_gd = single(rv_gd(:));
Matlab_gt_write_analyze(rv_gd, CreateGtImageHeader(rv_gd), fullfile(utDir, 'rv_gd'));

rv_rec = [rv_foot rv_valley rv_dose];
lv_rec = [lv_foot lv_valley lv_dose];
rv_rec = single(rv_rec(:));
Matlab_gt_write_analyze(rv_rec, CreateGtImageHeader(rv_rec), fullfile(utDir, 'rv_rec'));
lv_rec = single(lv_rec(:));
Matlab_gt_write_analyze(lv_rec, CreateGtImageHeader(lv_rec), fullfile(utDir, 'lv_rec'));

% find time of peak
lv_peak_time  = time_interp(min(find(lv_gd_est == max(lv_gd_est))));
rv_peak_time  = time_interp(min(find(rv_gd_est == max(rv_gd_est))));
gd_max = 1.15* max(max(lv_gd_est),max(rv_gd_est));
ptt = lv_peak_time - rv_peak_time;

% calculate duration (FWHM) of RV and LV
lv_duration  = time_interp(max(find(lv_gd_est >= max(lv_gd_est)/2))) - time_interp(min(find(lv_gd_est >= max(lv_gd_est)/2)));
rv_duration  = time_interp(max(find(rv_gd_est >= max(rv_gd_est)/2))) - time_interp(min(find(rv_gd_est >= max(rv_gd_est)/2)));

% calculate center of gravity (centroid)
lv_centroid = sum(time_interp .* lv_gd_est) / sum(lv_gd_est);
rv_centroid = sum(time_interp .* rv_gd_est) / sum(rv_gd_est);
ptt_centroid = lv_centroid - rv_centroid;

if ptt < 2
    return
end

% plot fits
hfig = figure; hold on
plot(time,[rv_gd],'b','linewidth',1)
plot(time_interp,[rv_gd_est],'k','linewidth',1)
plot(time,[lv_gd],'r','linewidth',1)
plot(time_interp,[lv_gd_est],'k','linewidth',1)
plot([lv_peak_time lv_peak_time],[-1 gd_max],'k:')
plot([rv_peak_time rv_peak_time],[-1 gd_max],'k:')
plot([lv_centroid lv_centroid],[-1 gd_max],'k:','linewidth',2)
plot([rv_centroid rv_centroid],[-1 gd_max],'k:','linewidth',2)
axis([0 max(time_interp) 0 gd_max])
set(hfig,'position',[194         519        1630         979]);
box on
title(['PTT_c_e_n_t_r_o_i_d  = ',num2str(ptt_centroid),' / PTT_p_k_-_p_k  = ',num2str(ptt),' / RV duration  = ',num2str(rv_duration),' / LV duration  = ',num2str(lv_duration),'  sec'],'fontsize',16)
xlabel('time (sec)','fontsize',14)
ylabel('[Gd], (mmol/L)','fontsize',14)
figurefilename = [odir, filesep, strrep(basefilename,'-','_'),'_rv_lv_fits'];

saveas(hfig, figurefilename,'fig')


hfig_rr = figure;
plot(1000*rr,'.'); axis([0 length(rr) 0 2500])
figurefilename_rr = [odir, filesep, strrep(basefilename,'-','_'),'_rr_intervals'];
xlabel('measurement #','fontsize',14)
ylabel('RR interval (ms)','fontsize',14)
saveas(hfig_rr, figurefilename_rr,'fig')
close(hfig_rr)


if 0
    % remove 1st pass from LV signal and plot residual(2nd pass)
    lv_gd_est2 = my_lognpdf(time,fit(1),fit(2),fit(3),fit(4));
    % hfig = figure; hold on
    % plot(time,[lv_gd],'b','linewidth',1)
    % plot(time,[lv_gd_est2],'k','linewidth',1)
    residual = lv_gd - lv_gd_est2;
    % plot(time,[residual],'r','linewidth',1)
    % box on
    % lv_dose_2nd_pass = 0.5 * sum(residual(round(lv_foot):end));
    fit = lsqcurvefit(func,[1 1  time(round(lv_valley)) 2*lv_dose], time(round(lv_valley):end), residual(round(lv_valley):end));
    lv_gd_2nd_pass_est = my_lognpdf(time,fit(1),fit(2),fit(3),fit(4));
    % plot(time,[lv_gd_2nd_pass_est],'m','linewidth',1); shg
    hfig = figure; hold on
    plot(time,[lv_gd],'b','linewidth',1)
    plot(time,[lv_gd_est2 + lv_gd_2nd_pass_est],'k','linewidth',1)
    box on
    title(['AIF fit'],'fontsize',16)
    xlabel('time (sec)','fontsize',14)
    ylabel('[Gd], (mmol/L)','fontsize',14)
    figurefilename = [odir, filesep, strrep(basefilename,'-','_'),'_aif_fit'];
    saveas(hfig, figurefilename,'fig')
end

out.time = time;
out.time_interp = time_interp;
out.aif = cat(4,aif1,aif2);
out.rv.rv_mask = rv_mask;
out.rv.rv1 = rv1;
out.rv.rv2 = rv2;
out.rv.rv_dose = rv_dose;
out.rv.rv_gd = rv_gd;
out.rv.rv_gd_est = rv_gd_est;
out.rv.rv_foot =  rv_foot;
out.rv.rv_peak =  rv_peak;
out.rv.rv_valley =  rv_valley;
out.rv.rv_peak_time = rv_peak_time;
out.rv.rv_duration = rv_duration;
out.rv.rv_centroid_time = rv_centroid;

out.lv.lv_mask = lv_mask;
% out.lv.lv1 = lv1;
% out.lv.lv2 = lv2;
out.lv.lv_dose = lv_dose;
out.lv.lv_gd = lv_gd;
out.lv.lv_gd_est = lv_gd_est;
out.lv.lv_foot =  lv_foot;
out.lv.lv_peak =  lv_peak;
out.lv.lv_valley =  lv_valley;
out.lv.lv_peak_time = lv_peak_time;
out.lv.lv_duration = lv_duration;
out.lv.lv_centroid_time = lv_centroid;

out.lv.lv_mask_size = lv_mask_size;
out.ptt = ptt;
out.ptt_centroid = ptt_centroid;
% out.rv.tc = sum(out.time_interp.*out.rv.rv_gd_est)/sum(out.rv.rv_gd_est);
% out.lv.tc = sum(out.time_interp.*out.lv.lv_gd_est)/sum(out.lv.lv_gd_est);
% out.ptt2 = out.lv.tc - out.rv.tc;
out.rr = 1000*rr; %(ms) store individual rr's to assess ECG quality and heart block
out.RR = RR;  % 
out.HR = 60/RR;

pmufilename = [odir,basefilename,'_pmudata.mat'];
if exist(pmufilename)
    tmp = load(pmufilename);
    rr = diff(tmp.pmudata.t*2.5);
    RR_pmu = median(rr);
    HR_pmu = 60000/RR_pmu;
    out.RR_pmu = RR_pmu;
    out.HR_pmu = HR_pmu;
end

    

% read h5file to get total Gd dose administered
% h5filename = 'Z:\RawData\Barts\20190831\Perfusion_AIFR3_2E_Interleaved_42110_488936486_488936495_43_20190831-083713.h5';
patientID = [];
protocol = [];
try
    [patientID, protocol,ismrmrd_header] = read_ismrmrd_protocol(h5filename);
    Gd_volume= [];
    if isfield(ismrmrd_header,'userParameters')
        if isfield(ismrmrd_header.userParameters,'userParameterDouble')
            if length(ismrmrd_header.userParameters.userParameterDouble)<8
                Gd_volume = [];
            else
                fieldname = ismrmrd_header.userParameters.userParameterDouble(8).name;
                if strfind(fieldname,'ContrastBolusVolume')
                    Gd_volume = ismrmrd_header.userParameters.userParameterDouble(8).value;
                end
            end
            if isfield(ismrmrd_header.userParameters, 'userParameterString')
                CA_type = ismrmrd_header.userParameters.userParameterString.value;
            else
                CA_type = 'unknown';
            end
        end
    end
    switch lower(CA_type)
        case {'gadavist', 'gadovist'}
                mmol_per_ml = 1.0;
        case {'dotarem', 'doterem'}
                mmol_per_ml = 0.5;
        otherwise
                mmol_per_ml = 0.5;
    end  
    out.CA_type = CA_type;
    out.mmol_per_ml = mmol_per_ml;
    
    if length(ismrmrd_header.userParameters.userParameterDouble) >=12
        out.age = ismrmrd_header.userParameters.userParameterDouble(12).value;
    else
        out.age = 0;
    end
    if length(ismrmrd_header.userParameters.userParameterDouble) >=12
        out.PatientHeight = ismrmrd_header.userParameters.userParameterDouble(11).value/10; % cm
    else
        out.PatientHeight = 0;
    end    
    if length(ismrmrd_header.userParameters.userParameterDouble) >=12
        out.PatientWeight = ismrmrd_header.userParameters.userParameterDouble(10).value; % kg
    else
        out.PatientWeight = 0;
    end        
    out.BSA = sqrt(out.PatientHeight * out.PatientWeight/3600); % m^2
    
    Indicator_Dose = [];
    if ~isempty(Gd_volume)
        Indicator_Dose = Gd_volume  * mmol_per_ml; % 0.5 mmol/ml x ml = mmol
    else
        Indicator_Dose = 0;
    end

    Cardiac_Output_RV_estimate = [];
    Cardiac_Output_LV_estimate = [];
%     HCT = 0.42;
    if ~isempty(Indicator_Dose)
        % Stroke Volume (SV) L/beat
        out.SV_rv = Indicator_Dose / sum(rv_gd_est);
        out.SV_lv = Indicator_Dose / sum(lv_gd_est);
        % estimate cardiac output
        dt = time_interp(2)-time_interp(1); % seconds sample interval
        % Cardiac_Output = Indicator_Dose * 60 / Integrated_estimate * (1-HCT);
        % units are Liters per minute (L/min)
        out.Cardiac_Output_RV_estimate = Indicator_Dose * 60 / (dt*sum(rv_gd_est));
        out.Cardiac_Output_LV_estimate = Indicator_Dose * 60 / (dt*sum(lv_gd_est));

        out.Pulmonary_Blood_Volume = out.Cardiac_Output_LV_estimate .* ptt * 1000/60;  %ml
        out.Pulmonary_Blood_Volume_index = out.Pulmonary_Blood_Volume ./ out.BSA; % ml/m^2
        out.Pulmonary_Blood_Volume_centroid = out.Cardiac_Output_LV_estimate .* ptt_centroid * 1000/60;  %ml
        out.Pulmonary_Blood_Volume_index_centroid = out.Pulmonary_Blood_Volume_centroid ./ out.BSA; % ml/m^2
    end
catch
    disp(['could not read: ',h5filename])
end

out.basefilename = basefilename;
out.h5filename = h5filename;
out.idir = idir;
% out.assumed_HCT = HCT;
out.Gd_volume = Gd_volume;
out.Indicator_Dose = Indicator_Dose;
% out.SV.lv = SV_lv;
% out.SV.rv = SV_rv;
% out.Cardiac_Output_RV_estimate = Cardiac_Output_RV_estimate;
% out.Cardiac_Output_LV_estimate = Cardiac_Output_LV_estimate;
% out.Pulmonary_Blood_Volume = Pulmonary_Blood_Volume;  %ml
% out.Pulmonary_Blood_Volume_index = Pulmonary_Blood_Volume_index; % ml/m^2
% out.Pulmonary_Blood_Volume_centroid = Cardiac_Output_LV_estimate .* ptt_centroid * 1000/60;  %ml
% out.Pulmonary_Blood_Volume_index_centroid = Pulmonary_Blood_Volume_centroid ./ out.BSA; % ml/m^2
out.patientID = patientID;
out.protocol = protocol;

save([odir, filesep, basefilename,'_data'], 'out')

close(hfig);

end

function h5filename = find_h5filename(basefilename);

    % basefilename = 'Perfusion_AIFR3_2E_Interleaved_42110_488936486_488936495_43_20190831-083713';

    ind = strfind(basefilename,'_');

    siteid = basefilename(ind(end-4)+1:ind(end-3)-1);
    studydate = basefilename(ind(end)+1:ind(end)+8);


    switch siteid
        case '41837'
            sitename = 'BARTS';
        case '42110'
            sitename = 'BARTS';
        case '66016'
            sitename = 'BARTS';
        case '42363'
            sitename = 'ROYALFREE';
        case '141550'
            sitename = 'ROYALBROMPTON';
        case '42629'
            sitename = 'ROYALBROMPTON';
        case '169517'
            sitename = 'ROYALBROMPTON';        
        case '42498'
            sitename = 'ROYALBROMPTON';
        case '169507'
            sitename = 'ROYALBROMPTON';
        case '66097'
            sitename = 'LEEDS';
    end
    
    h5filename = ['T:\RawData\',sitename,filesep,studydate,filesep,basefilename,'.h5'];
end

function out = apply_lut_correction(in, lut);
% function out = apply_lut_correction(in, lut);
% 
% function to correct perfusion signal using lut calculated by Bloch
% equation smulation

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

% lut(:,1) is Gd
% lut(:,2) is normalized signal S = SR/PD

    out = interp1(lut(:,2),lut(:,1), in); % linear interpolated table look-up
  
end


function [patientID, protocol , ismrmrd_header] = read_ismrmrd_protocol(filename);
% function [patientID, protocol, ismrmrd_header] = read_ismrmrd_protocol(filename);

% uses Matlab Central base64decode.m
    try
        dset = ismrmrd.Dataset(filename, 'dataset');
    catch
        patientID = [];
        protocol = [];
        return
    end
    try
        hdr = ismrmrd.xml.deserialize(dset.readxml);
    catch
       hdr.userParameters = []; 
    end
    if isfield(hdr.userParameters,'userParameterBase64')
        str = hdr.userParameters.userParameterBase64.value;
    else
        patientID = [];
        protocol = [];
        ismrmrd_header = [];
        return    
    end


    if isfield(hdr.subjectInformation,'patientID')
        patientID = hdr.subjectInformation.patientID;
    else
        patientID = [];
    end

 
% if nargout > 2
%     if isfield(hdr.userParameters,'userParameterDouble')
%         patinfo.PatientAge = hdr.userParameters.userParameterDouble(12);
%     else
%         patinfo.PatientAge = [];
%     end
%     if isfield(hdr.userParameters,'userParameterDouble')
%         patinfo.PatientWeight = hdr.userParameters.userParameterDouble(10);
%     else
%         patinfo.PatientWeight = [];
%     end
%     if isfield(hdr.userParameters,'userParameterDouble')
%         patinfo.PatientHeight = hdr.userParameters.userParameterDouble(11);
%     else
%         patinfo.PatientHeight = [];
%     end
% end
    
    if nargout > 1
        header = char(base64decode(str));
        protocol = header;
    end
    if nargout > 2
        ismrmrd_header = hdr;
    end

end

