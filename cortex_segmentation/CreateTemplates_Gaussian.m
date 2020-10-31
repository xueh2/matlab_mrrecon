
function [csfT, wmT, cortexT] = CreateTemplates_Gaussian(header, csfSeg, wmSeg, cortexSeg,...
                    sigmaCSF, sigmaWM, sigmaCortex, halfwidthCsf, halfwidthWm, halfwidthCortex)
% create probability templates

csfT = double(csfSeg);
cortexT = double(cortexSeg);
wmT = double(wmSeg);

G_csfT = Gaussianfilter(csfT, header, sigmaCSF, sigmaCSF, sigmaCSF, halfwidthCsf, 'Real');
G_cortexT = Gaussianfilter(cortexT, header, sigmaCortex, sigmaCortex, sigmaCortex, halfwidthCortex, 'Real');
G_wmT = Gaussianfilter(wmT, header, sigmaWM, sigmaWM, sigmaWM, halfwidthWm, 'Real');

csfflag = 0;
wmflag = 0;
cortexflag = 0;

minP = 0.02;

for k = 1:header.zsize
    for j = 1:header.ysize
        for i = 1:header.xsize
            
            p_csf = G_csfT(j, i, k);
            p_wm = G_wmT(j, i, k);
            p_cortex = G_cortexT(j, i, k);
            
            if ( G_csfT(j, i, k)>minP )
                csfflag = 1;
            end
            
            if ( G_wmT(j, i, k)>minP )
                wmflag = 1;
            end
            
            if ( G_cortexT(j, i, k)>minP )
                cortexflag = 1;
            end

            if ( csfflag + wmflag + cortexflag == 0 )
                continue;
            end
            
            csfT(j, i, k) = p_csf/(p_csf+p_wm+p_cortex);
            wmT(j, i, k) = p_wm/(p_csf+p_wm+p_cortex);
            cortexT(j, i, k) = p_cortex/(p_csf+p_wm+p_cortex);
            
            csfflag = 0;
            wmflag = 0;
            cortexflag = 0;
        end
    end
end
return;
