
function [csfT, wmT, cortexT, outlierT] = CreateTemplates_4classes_Gaussian(header, csfSeg, wmSeg, cortexSeg, outlierSeg, ...
                    sigmaCSF, sigmaWM, sigmaCortex, sigmaOutlier, halfwidthCsf, halfwidthWm, halfwidthCortex, halfwidthOutlier)

disp('create templates ...');

% create probability templates

csfT = single(csfSeg);
cortexT = single(cortexSeg);
wmT = single(wmSeg);
outlierT = single(outlierSeg);

G_csfT = Gaussianfilter(csfT, header, sigmaCSF, sigmaCSF, sigmaCSF, halfwidthCsf, 'Real');
G_cortexT = Gaussianfilter(cortexT, header, sigmaCortex, sigmaCortex, sigmaCortex, halfwidthCortex, 'Real');
G_wmT = Gaussianfilter(wmT, header, sigmaWM, sigmaWM, sigmaWM, halfwidthWm, 'Real');
G_outlierT = Gaussianfilter(outlierT, header, sigmaOutlier, sigmaOutlier, sigmaOutlier, halfwidthOutlier, 'Real');

csfflag = 0;
wmflag = 0;
cortexflag = 0;
outlierflag = 0;

minP = 0.02;

for k = 1:header.zsize
    for j = 1:header.ysize
        for i = 1:header.xsize
            
            p_csf = G_csfT(j, i, k);
            p_wm = G_wmT(j, i, k);
            p_cortex = G_cortexT(j, i, k);
            p_outlier = G_outlierT(j, i, k);

            if ( G_csfT(j, i, k)>minP )
                csfflag = 1;
            end
            
            if ( G_wmT(j, i, k)>minP )
                wmflag = 1;
            end
            
            if ( G_cortexT(j, i, k)>minP )
                cortexflag = 1;
            end
            
            if ( G_outlierT(j, i, k)>minP )
                outlierflag = 1;
            end
            
            if ( csfflag + wmflag + cortexflag + outlierflag == 0 )
                continue;
            end
            
            csfT(j, i, k) = p_csf/(p_csf+p_wm+p_cortex+p_outlier);
            wmT(j, i, k) = p_wm/(p_csf+p_wm+p_cortex+p_outlier);
            cortexT(j, i, k) = p_cortex/(p_csf+p_wm+p_cortex+p_outlier);
            outlierT(j, i, k) = p_outlier/(p_csf+p_wm+p_cortex+p_outlier);
            
%             csfT(j, i, k) = 0.33;
%             wmT(j, i, k) = 0.33;
%             cortexT(j, i, k) = 0.33;
            
            csfflag = 0;
            wmflag = 0;
            cortexflag = 0;
            outlierflag = 0;
        end
    end
end
return;
