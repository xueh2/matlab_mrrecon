
%function [Fat, Water, InPhase, OppPhase, B0Map, T2StarMap] = FatWaterSeparation(MultiEchoImage, header, gyroRatio, B0, TEs, params)
% Perform the Fat Water Separation
% MultiEchoImage : multi-echo image [Nfe * Npe * numOfEcho] 
% gyroRatio : gyromagnetic ration in the unit of MHz/Tesla
% B0 : B0 strength in Tesla
% TEs : echo time for every echo in second
% params : parameters including the following
% -------------------------------------------------
% RANGE_FREQ_OFFSET:    the range of field map values to consider 
% NUM_FREQ_OFFSETS:     number of discretized field map values
% RANGE_T2S:            range of T2* values, in ms
% NUM_T2S:              the number of R2* values (1 if no decay is assumed)
% NUM_T2S_FINE:         the number of R2* values for final estimation (1 if no decay is assumed)
% GET_DAMPED_SIGNALS:   obtain signal estimates at TE1
% NUM_FAT_PEAKS:        number of fat peaks (should be 6 currently, note that some rel_amps can be set to zero)
% DELTA_F:              an array of size NUM_FAT_PEAKS (currently fixed at 6) with the chemical shifts (in ppm) of water and fat
% REL_AMPS:             complex valued relative amplitudes of fat peaks
% LAMDA:                regularization parameter
% LMAP_POWER:           controls spatially-varying regularization (2 gives ~ uniform resolution, 0 gives constant regularization parameter)
% LMAP_EXTRA:           additional spatially-nonvarying regularization (with LMAP_POWER=2, this gives more smoothing in noisy regions)
% MAXIT:                maximual number of random walker iterations
% -------------------------------------------------

%% perpare the processing
S0 = size(MultiEchoImage);
Nfe = S0(1);
Npe = S0(2);
numOfEcho = S0(3);

RANGE_FREQ_OFFSET = params.RANGE_FREQ_OFFSET;
NUM_FREQ_OFFSETS = params.NUM_FREQ_OFFSETS;
freqStep = (RANGE_FREQ_OFFSET(2)-RANGE_FREQ_OFFSET(1))/NUM_FREQ_OFFSETS;
freqs = [RANGE_FREQ_OFFSET(1):freqStep:RANGE_FREQ_OFFSET(2)];
numOfFreqs = length(freqs);

RANGE_T2S = params.RANGE_T2S;
NUM_R2S = params.NUM_T2S_FINE;
R2S_step = (1/RANGE_T2S(1)-1/RANGE_T2S(2))/NUM_R2S;
R2S = [1/RANGE_T2S(2) :R2S_step: 1/RANGE_T2S(1)];
R2S = R2S * 1e-3; % ms -> second
numOfR2S = length(R2S);

fatDeltaFreq = gyroRatio*B0*params.DELTA_F; % in Hz

tstart = tic;

%% step 0
% pre-compute the projection matrix
tic

% build the phi matrix
phi = zeros(numOfEcho, 2);
phi(:,1) = 1;
for eco=1:numOfEcho
    phi(eco,2) = sum(params.REL_AMPS.*exp(i*2*pi*fatDeltaFreq*TEs(eco)));
end
phiInv = pinv(phi);

% build the lamda/phi/proj matrix, for every freq and one R2*, one matrix is pre-computed
lamdaAll = zeros(numOfEcho, numOfEcho, numOfFreqs, numOfR2S+1);
lamdaInvAll = zeros(numOfEcho, numOfEcho, numOfFreqs, numOfR2S+1);

psiAll = zeros(numOfEcho, 2, numOfFreqs, numOfR2S+1);
psiInvAll = zeros(2, numOfEcho, numOfFreqs, numOfR2S+1);
projAll = zeros(numOfEcho, numOfEcho, numOfFreqs, numOfR2S+1);

for r=1:numOfR2S
    r2s = R2S(r);
    for f=1:numOfFreqs
        fb = freqs(f);
        value = (-r2s+i*2*pi*fb).*TEs;
        lamdaAll(:,:,f,r) = diag(exp(value));
        lamdaInvAll(:,:,f,r) = diag(exp(-1*value));
        psiAll(:,:,f,r) = lamdaAll(:,:,f,r)*phi;
        psiInvAll(:,:,f,r) = phiInv*lamdaInvAll(:,:,f,r);
        projAll(:,:,f,r) = eye(numOfEcho) - psiAll(:,:,f,r)*psiInvAll(:,:,f,r);
    end
end

% last R2S, no decay
for f=1:numOfFreqs
    fb = freqs(f);
    value                           = (i*2*pi*fb).*TEs;
    lamdaAll(:,:,f,numOfR2S+1)      = diag(exp(value));
    lamdaInvAll(:,:,f,numOfR2S+1)   = diag(exp(-1*value));
    psiAll(:,:,f,numOfR2S+1)        = lamdaAll(:,:,f,numOfR2S+1)*phi;
    psiInvAll(:,:,f,numOfR2S+1)     = phiInv*lamdaInvAll(:,:,f,numOfR2S+1);
    projAll(:,:,f,numOfR2S+1)       = eye(numOfEcho) - psiAll(:,:,f,numOfR2S+1)*psiInvAll(:,:,f,numOfR2S+1);
end

clear lamdaAll lamdaInvAll psiAll
disp(['Step 0 (Pre-compute the matrix): ' num2str(toc)])

%% step 1
% Assume no T2* decay, estimate the pixel-wise field map by 1D brutal search
tic

% for every pixel, find the optimal field map
R0 = zeros(Nfe, Npe, numOfFreqs);
signal = reshape(MultiEchoImage, [Nfe*Npe numOfEcho]);
signal = signal';
for f=1:numOfFreqs
    % projF = proj(:,:,f);
    projF = projAll(:,:,f,numOfR2S+1);
    R0_F = projF*signal; % [numOfEcho Nfe*Npe]
    R0_F = sum(R0_F.*conj(R0_F), 1);
    R0(:,:,f) = reshape(R0_F, [Nfe Npe]);
end

[minR0, ind] = min(R0, [], 3);
B0Map = freqs(ind);

clear R0

% compute fat and water image
Fat = zeros(Nfe, Npe);
Water = zeros(Nfe, Npe);

for pe=1:Npe
    for fe=1:Nfe        
        psiInv_pixel = psiInvAll(:,:,ind(fe, pe),numOfR2S+1);
        rou = psiInv_pixel*signal(:, fe+(pe-1)*Nfe);
        Water(fe,pe) = rou(1);
        Fat(fe,pe) = rou(2);
    end
end

B0Map_pixel_wise = B0Map;
Water_pixel_wise = Water;
Fat_pixel_wise = Fat;

figure; imshow(abs(Water), [], 'InitialMagnification', 300);
figure; imshow(B0Map, [], 'InitialMagnification', 300);

disp(['Step 1 (Estimate the pixel-wise field map by 1D brutal search): ' num2str(toc)])

%% step 2 : pre-compute the graph and related matrixs
% compute T
T = sparse(1:Nfe*Npe, 1:Nfe*Npe, 8*pi*pi*(TEs'*TEs), Nfe*Npe, Nfe*Npe);

% compute combinatorial Laplacian matrix L
% L : Nfe*Npe by Nfe*Npe

[points, edges] = lattice(Nfe,Npe,2);
weights = ones([1 size(edges, 1)]);
L = laplacian(edges,weights);

% % compose the L matrix
% m = 8*Nfe*Npe;
% n = Nfe*Npe;
% 
% A_Inc = sparse(m, n);
% C = sparse(1:m, 1:m, 1, m, m);
% 
% % fill incidence matrix
% indE = 1;
% edgeI = zeros(m,1);
% edgeJ = zeros(m,1);
% % go through edges
% for pe2=1:Npe
%     for fe2=1:Nfe
%         i_ind = fe2+(pe2-1)*Nfe;
%         % go through the 8 neighbors
%         for pe3=pe2-1:pe2+1
%             for fe3=fe2-1:fe2+1
%                 if ( pe3==pe2 & fe3==fe2 )
%                     continue;
%                 end
%                 j_Ind = fe3+(pe3-1)*Nfe;
%                 edgeI(indE) = i_ind;
%                 edgeJ(indE) = i_ind;
%                 indE = indE + 1;
%             end
%         end
%     end
% end
% 
% for pe=1:Npe
%     for fe=1:Nfe        
%         k = fe+(pe-1)*Nfe;
%         indI = find(edgeI==k);
%         indJ = find(edgeJ==k);
%         A_Inc(indI(:), k) = 1;
%         A_Inc(indJ(:), k) = -1;
%     end
% end
% 
% L = A_Inc'*C*A_Inc;

A = T - 2*params.LAMDA*L;

% initial R2*
R2SMap = zeros(Nfe, Npe);

deltaChangeWater = zeros(params.MAXIT, 1);
deltaChangeMap = zeros(params.MAXIT, 1);

B0Map = B0Map_pixel_wise;
Water = Water_pixel_wise;
Fat = Fat_pixel_wise;

% B0Map = B0Map_Rep';
% Water = Water_Rep';
% Fat = Fat_Rep';
 
figure; imshow(abs(Water), []);drawnow;
figure; imshow(abs(Fat), []);drawnow;

for iter=1:params.MAXIT
    iter
    
    
    %% step 3 : solve random walker equation to get an updated freq map
    tic
    % 1) compute every item of S, column-wise storage
    S = zeros(Nfe, Npe);
    for pe=1:Npe
        for fe=1:Nfe               
            % s = signal(:, fe+(pe-1)*Nfe);
            s = MultiEchoImage(fe, pe, :);
            s = reshape(s, [numOfEcho 1]);
            d_q = log(s) + R2SMap(fe,pe).*TEs - log(Water(fe,pe)+Fat(fe,pe).*phi(:,2));
            S(fe, pe) = -4*pi*i*(TEs'*d_q);
        end
    end
        
    % 5) call the sparse matrix solver  
%     A = T - 2*LAMDA*L;
      x = A\S(:);
    
%     tol = 1e-4;
%     maxit = 30;
%     M1 = [];
%     M2 = [];
%     x0 = B0Map(:);
%     
%     [x,flag] = lsqr(A,S,tol,maxit,M1,M2,x0);
%     flag
   
    x0 = B0Map(:);
    deltaChangeMap(iter) = norm(x(:)-x0(:));
    deltaChangeMap(iter)
    
    % update freq map estimation
    B0Map = reshape(x, [Nfe Npe]); 
    imshow(abs(B0Map), [], 'InitialMagnification', 300);
    disp(['Step 3 (Solve random walker equation to get an updated freq map): ' num2str(toc)])
    
    %% step 4: Update the R2* estimation, given the estimated B0Map
%     tic
%     % for every pixel, find the optimal R2*
%     R0_Update = zeros(Nfe, Npe, numOfR2S);
%     for pe=1:Npe
%         tic
%         for fe=1:Nfe        
% 
%             fb = B0Map(fe, pe);
%             fb = abs(fb);
%             
%             % find the index
%             ind_Pixel = round( numOfFreqs *(fb-RANGE_FREQ_OFFSET(1))/(RANGE_FREQ_OFFSET(2)-RANGE_FREQ_OFFSET(1)) );
%             if ( ind_Pixel <= 0 )
%                 ind_Pixel = 1;
%             end
% 
%             if ( ind_Pixel >= numOfFreqs )
%                 ind_Pixel = numOfFreqs;
%             end
% 
%             for r=1:numOfR2S            
%                 proj_Pixel = projAll(:,:,ind_Pixel, r);
%                 R0_Pixel = proj_Pixel*signal(:, fe+(pe-1)*Nfe);            
%                 R0_Pixel = sum(R0_Pixel.*conj(R0_Pixel));
%                 R0_Update(fe, pe, r) = R0_Pixel;
%             end 
%         end
%         toc
%     end
% 
%     % for every pixel, find the optimal R2*
%     [minR0_Update, indR2S] = min(R0_Update, [], 3);
%     R2SMap = R2S(indR2S);
%     clear minR0_Update
    
    disp(['Step 3 (Estimate the pixel-wise R2* by 1D brutal search): ' num2str(toc)])

    %% step 4 : Update the water and fat estimation
    prevWater = Water;
    tic
    for pe=1:Npe
        for fe=1:Nfe 
            fb = B0Map(fe, pe);
            r2s = R2SMap(fe, pe);
            lamdaInv_Pixel = diag(  exp( -1*(-r2s+i*2*pi*fb).*TEs )  );
            psiInv_Pixel = phiInv*lamdaInv_Pixel;
            rou = psiInv_Pixel*signal(:, fe+(pe-1)*Nfe);
            Water(fe,pe) = rou(1);
            Fat(fe,pe) = rou(2);
        end
    end

    disp(['Step 4 (Update the water and fat estimation): ' num2str(toc)])

    figure;imshow(abs(Water), [], 'InitialMagnification', 300);
    
    deltaChangeWater(iter) = norm(Water(:) - prevWater(:));
    deltaChangeWater(iter)
    
    disp('next iteration ... ');
end

%% step 5 : compute signal estimates at the first TE, correct fat shift, etc.
InPhase = Fat + Water; 
OppPhase = Water - Fat;
T2StarMap = 1./(R2SMap+eps);

%% finish
disp(['Entire process : ' num2str(toc(tstart))])

%% plot
figure; plot(deltaChangeMap);
figure; plot(deltaChangeWater);

