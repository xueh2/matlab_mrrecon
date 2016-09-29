function R = ssfp_readout(params);
% function R = ssfp_readout(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize input paramaters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npulses = params.Npulses;
alpha_deg = params.flip_angle;
TR1 = params.TR;
steadystate_prep = params.steadystate_prep;
N_runup = params.N_runup;

T1 = params.T1;
T2 = params.T2;

off_reson_vector = params.offresonance; % Hz
num_iso = length(off_reson_vector);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize functional and operational references
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = Mops_defs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%								
alpha     = M.rad(alpha_deg);           % readout flip angle (deg to radians)

% Define other constants
% M0 = [ 0 0 1 1]';                % Equilibrium magnetization
TE1 = TR1/2;                     % Echo time, s

% Define the precession angles for TR1 and TR2
theta_per_TR1     = 2*pi*off_reson_vector * TR1;
theta_per_TE1     = theta_per_TR1/2;

% Create Operation matrices
R1_rf_odd     = M.RFHard(alpha);                             % Assume instantanoues hard RF
R1_rf_even    = M.RFHard(-alpha);                            % Assume instantanoues hard RF
R_relax1      = M.Relax(TE1, T1, T2);                        % Relaxation that takes place during one TE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:num_iso
	
    % set-up run up (catalization of steady state)
    for k=1:length(alpha)
        
        switch steadystate_prep
            case 'alpha/2'
                alpha_prep(k) = - alpha(k)/2; % alpha/2 prep;
                R1_ssfp_startup = M.RFHard(alpha_prep(k));
            case {'linear'}
                N_runup_even=mod(N_runup+1,2); %should be 1 when N_runup is even
                R1_ssfp_startup = eye (4);
                R_Precess1 = M.Precess(theta_per_TE1(i));
                R1_pr = R_Precess1 * R_relax1 ;
                for preppulse=1:N_runup
                    alpha_prep(k)=(2*preppulse-1)*alpha(k)/(2*N_runup)*((-1)^(preppulse+N_runup_even));
                    R1_temp      = M.RFHard(alpha_prep(k));
                    if preppulse < N_runup
                        R1_ssfp_startup = R1_pr * R1_pr * R1_temp * R1_ssfp_startup;
                    else
                        R1_ssfp_startup = R1_pr * R1_temp *R1_ssfp_startup;
                    end
                end

            case 'kaiser-bessel'
                N_runup_even=mod(N_runup+1,2); %should be 1 when N_runup is even
                R1_ssfp_startup = eye (4);
                R_Precess1 = M.Precess(theta_per_TE1(i));
                R1_pr = R_Precess1 * R_relax1 ;

                beta = 0;
                dFlipAngle = kb(alpha(k), N_runup, beta);
                for preppulse=1:N_runup
                    alpha_prep(k)=dFlipAngle(preppulse)*((-1)^(preppulse+N_runup_even));
                    R1_temp      = M.RFHard(alpha_prep(k));
                    if preppulse < N_runup
                        R1_ssfp_startup = R1_pr * R1_pr * R1_temp * R1_ssfp_startup;
                    else
                        R1_ssfp_startup = R1_pr * R1_temp *R1_ssfp_startup;
                    end
                end
        end
    end
    
    % set-up readout (R_readout)
    R_Precess1 = M.Precess(theta_per_TE1(i));
    R1_pr = R_Precess1 * R_relax1 ;
    R_readout_odd = R1_pr * R1_rf_odd * R1_pr;
    R_readout_even = R1_pr * R1_rf_even * R1_pr;
    R_readout = cat (3, R_readout_odd, R_readout_even);

    % set-up store
    %     alpha_prep(k) = (-1)^(Npulses)* alpha(k)/2; % alpha/2 prep;
    %     R1_ssfp_store = M.RFHard(alpha_prep(k));

    % apply operators to create output

    % RF prep
%     R = R1_ssfp_startup; % at echo time
% 
%     % readout                                      	
%     for n=1:Npulses
%         R(:,:,i) = R_readout(:,:,mod(n+1,2)+1) * R; % R stored for i-th frequency (isochromat)
%     end

    R0 = R1_ssfp_startup; % at echo time

    % readout                                      	
    for n=1:Npulses
        R0 = R_readout(:,:,mod(n+1,2)+1) * R0; % R stored for i-th frequency (isochromat)
    end
    
    R(:,:,i) = R0;
end







