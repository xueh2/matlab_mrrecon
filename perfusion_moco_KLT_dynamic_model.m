function mocoRes = perfusion_moco_KLT_dynamic_model(data, Ns, sigma, numOfPD)
% perfusion moco using KLT dynamic modeling strategy

a = data(:,:,numOfPD+1:end);
N = size(a, 3)
sa = size(a);

% Ns = [N floor(N/2) floor(N/4) 16 8];
% sigma = [16 14 12 10 8];
warped = a;

warpedAll = zeros(sa(1), sa(2), sa(3), numel(Ns));
modelAll = zeros(sa(1), sa(2), sa(3), numel(Ns));

for iter=1:numel(Ns)
    widthForModel = Ns(iter);
    ratio = 0.02;
    modelSeries = KLT_dynamic_model(warped, widthForModel, ratio);

    modelAll(:,:,:,iter) = modelSeries;

    sigmas = repmat(sigma(iter), [1 3]);
    [dx, dy, warped, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series_pairwise(single(modelSeries), single(a), 'LocalCCR', 3, [16 64 100], sigmas, [2 2 2 2], 1, 0, 3, 10, 0.5, [], 0);    
    
    header = CreateGtImageHeader(dx);
    [jac, meanLogJac, maxLogJac] = Matlab_gt_deformation_field_jacobian_2D_image(dx, dy, header);
    
    ind = find(abs(maxLogJac)>0.65 | meanLogJac>0.1);
    if(~isempty(ind))
        if ( iter > 1 )
            warped(:,:,ind(:)) = warpedAll(:,:,ind(:), iter-1);
        else
            warped(:,:,ind(:)) = a(:,:,ind(:));
        end
    end
    
    warpedAll(:,:,:,iter) = warped;
end

aAll = repmat(a, [1 1 1 numel(Ns)]);
mocoRes = cat(4, aAll, modelAll, warpedAll);
