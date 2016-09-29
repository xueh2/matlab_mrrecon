function [M] = flash_readout(params);
% function [M] = flash_readout(params);

%% First method in Sekihara
% Define operators


    Mops = Mops_defs;
    %% Define Constants
%     TE = params.TE;                  % Echo Time, s
    TR = params.TR;                  % Repetition Time, s
    T1 = params.T1;                  % T1,T2 for gray matter,s
    T2 = params.T2;                   
    flip_angle=params.flip_angle;    % Flip angle, deg
    offR = params.offresonance;      % Off-resonance, Hz  
    nRF = params.Npulses;            % Number of RF pulses
    rf_phase_spoiler_increment = params.rf_phase_spoiler_increment; % degrees

    nIso = 71;                    % Number of Isochromats per voxel to sim

    % Calculate other values
    theta_per_TR = 2*pi*offR*TR;              % Precession per TR due to off-resonance, rad

    rf_phase_spoiler_increment_radians = rf_phase_spoiler_increment*pi/180; % RF spoiling increment
    crush_angles = linspace(0, 2*pi, nIso+1);
    crush_angles = crush_angles(1:end-1);

    % Generate Block-Diagonal Operators (const through all sim)
    % Relaxation over 1 TR
    RelaxTR_Op   = Mops.Relax(TR,T1,T2);

    % Hard (ideal) RF Excitation
    RF_Op = Mops.RFHard(Mops.rad(-flip_angle),0);

    if ~isfield(params,'M0')
        M0 = [0 0 1 1]';
    else
        M0 = params.M0;
    end

    % Off resonance precession
    OffR_TR_Op = Mops.Precess(theta_per_TR);

    % Define magnetization matrices before (Mminus) and after (Mplus) the RF pulse
    % Also make a cell to hold M, one for each phase increment
    %                    M0 + nHB*NTRs 
    Mp = zeros(size(M0,1), 1 + nRF, nIso);
    Mm = Mp;

    % nIso Ischromats distributed over 2pi (to be summed/integrated)
    for iso = 1:nIso
        Iso_Op = RelaxTR_Op * Mops.Precess(crush_angles(iso));

        % Note: the "TR operator" is cumulative so it is applied to M0 at each
        % recorded time step
        phase_increment = 0;
        PI_Op = Mops.Precess(phase_increment); 

        % Initialize unit operators and temp mag storage 
        t = 1;
        TR_Op = eye(4);
        Mp(:,t,iso) = M0; 
        Mm(:,t,iso) = M0; t=t+1;
        % Cycle over nRF TRs to generate a single cumulative operator for each nIso
        % Increment phase operator (PI_Op)
        for n = 1:nRF
            % magnetization before and after the RF pulse 
            Mm(:,t,iso) = TR_Op * M0;
            Mp(:,t,iso) = RF_Op * TR_Op * M0; t=t+1;

            % TR
            TR_Op = OffR_TR_Op * PI_Op * Iso_Op * RF_Op * TR_Op;
            phase_increment = phase_increment + rf_phase_spoiler_increment_radians;
            PI_Op = Mops.Precess(phase_increment);
        end		
    end

    % Integrate over all isochromats (nIso) to determine voxel signal
    Mp_voxel = mean(Mp,3);

    M = Mp_voxel;

    M = M(:,2:end); % first value is initial condition

return

	


    



