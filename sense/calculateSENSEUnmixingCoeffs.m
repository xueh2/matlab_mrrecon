function [UnmixingCoeffs] = calculateSENSEUnmixingCoeffs(b1map, snr_map, params, swap_freq_phase);
% function [UnmixingCoeffs] = calculateSENSEUnmixingCoeffs(b1map, snr_map, params, swap_freq_phase);
%
% function to compute the SENSE un-mixing matrix
%   UnmixingCoeffs = inv(S'S + lambda)*S';
%     (least squares inverse solution using SNR scaled regularization)
%      assumes that the data and b1map have been noise-prewhitened
%
% inputs:
%	b1map:    complex multi-coil sensitivity profile, (LIN,COL,CHA)
%             assumes b1map are normalized by root-sum-of-squares (i.e.,
%             rss(b1map)=1 );
%	snr_map:  reference image for adaptive regularization in SNR units
%   params:  structure containing
%       params.PAT.SENSE.regularization_parameter, sets degree of regularization
%       params.sense.AccelFactPE,    PE acceleration factor
% outputs:
%	UnmixingCoeffs      SENSE un-mixing coefficients (LIN,COL,CHA) over full FOV
%

alpha = params.PAT.SENSE.regularization_parameter; % typically alpha <<1
R     = params.PAT.AccelFactPE;

if isfield(params.PAT.SENSE,'regularization_mode');
    regularization_mode = params.PAT.SENSE.regularization_mode;
else
    regularization_mode = 'square'; % Default
end

S=b1map; % complex coil sensitivity estimates

if nargin < 4; swap_freq_phase =0; end % default
if swap_freq_phase ==1 % swaps rows and columns
    S = permute(S,[2 1 3 4]);
    snr_map = permute(snr_map,[2 1]);
end


% set-up sensitivity matrix
[rows,cols,ncoils]=size(S);
for k=1:R
	s(:,:,:,k)=S(((k-1)*rows/R+1):k*rows/R,:,:);
end

sh=conj(permute(s,[1 2 4 3]));

ShS=matrix_mult_image(sh,s);
I=permute(repmat(eye(R),[1 1 rows/R cols]),[3 4 1 2]);




lambda=zeros(rows/R,cols,R,R,'single');
for i=1:R    
    switch regularization_mode
        case 'square'
            lambda(:,:,i,i)=alpha./(snr_map(((i-1)*rows/R+1):i*rows/R,:).^2+eps);
        case 'linear'
            lambda(:,:,i,i)=alpha./(snr_map(((i-1)*rows/R+1):i*rows/R,:)+eps);
    end
end

if R<9
    ShSinv=inv_nxn(ShS + lambda);
else % LOOP OVER ROWS & COLS
    tmp = ShS + lambda;
    tmp = permute(tmp,[3 4 1 2]);
    ShSinv = zeros(size(ShS),'single');
    for row=1:size(ShS,1)
        for col=1:size(ShS,2)
            ShSinv(row,col,:,:)=inv(tmp(:,:,row,col));
        end
    end
end

u=matrix_mult_image(ShSinv,sh);

UnmixingCoeffs=zeros(size(S),'single');
for k=1:R
	UnmixingCoeffs(((k-1)*rows/R+1):k*rows/R,:,:)=squeeze(u(:,:,k,:));	
end

if swap_freq_phase ==1 % swaps rows and columns
    UnmixingCoeffs = permute(UnmixingCoeffs,[2 1 3 4]);
end

return


