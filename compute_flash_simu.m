function [mag_all, mag_center] = compute_flash_simu(T1, T2, B1, TR, FA, E1, accel_factor)
% T1, T2, TR: in ms

num_tubes = numel(T1);
num = numel(FA);
mag_center = zeros(num_tubes, num);

mag_all = zeros(E1, num_tubes, num);

for tt=1:num_tubes
    disp(['processing tube ' num2str(tt) ' ... ']);
    for fa=1:num
       
        params.T1 = T1(tt)/1e3;
        params.T2 = T2(tt)/1e3;
        params.TR = TR/1e3;
        params.flip_angle = FA(fa)*B1(tt);
        params.offresonance = 0;
        params.Npe_full = E1;
        params.PAT_accel_factor =  accel_factor;
        params.steadystate_prep = 'linear';
        params.N_runup = 3;
        params.rf_phase_spoiler_increment = 112;

        Npe_total = params.Npe_full/params.PAT_accel_factor;
        Npe_center = floor(Npe_total/2)+1;
        params.Npulses = E1;

        disp(['processing flip angle : ' num2str(params.flip_angle)]);
        params
        
        [PD] = flash_readout(params);
        
        
        mag = squeeze(rss(PD(1:2,:),1));
        mag_all(:, tt, fa) = mag(:);
        center = mean(mag(Npe_center-2:Npe_center+2))
        
        mag_center(tt, fa) = center;
    end
end
