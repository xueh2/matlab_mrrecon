function senMap = CoilSensitivityMapEstimation(complexImageSen, csmMethod, spatialSmoothingKernel, kSizeEigenVectorCoilSensitivity, percentageEigenVectorCoilSensitivity, kSize)
% method: 'Walsh', 'Jun', 'Souheil'

senMap = complexImageSen;
REP = size(complexImageSen, 4);

for rep=1:REP

    if ( strcmp(csmMethod, 'Walsh') )
        params.b1map.spatial_smoothing = 1;
        params.b1map.coilmap_smoothing_blocksize = spatialSmoothingKernel;
        senMap(:,:,:,rep) = calculateB1map(complexImageSen(:,:,:,rep), params);
    end

    if ( strcmp(csmMethod, 'Jun') )
        opts.choice = 1;
        opts.percentage = percentageEigenVectorCoilSensitivity;
        tic
        [senMap(:,:,:,rep), eigD, tau]=ED_eigen_2D(fft2c(complexImageSen(:,:,:,rep)), kSizeEigenVectorCoilSensitivity, size(complexImageSen), opts);
        toc
    end

end

if ( strcmp(csmMethod, 'Souheil') )
    tic
    senMap = CoilSensitivity_Souheil_Walsh(complexImageSen, kSize, 0);
    toc
end
