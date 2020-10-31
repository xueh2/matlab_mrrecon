
function [wmT, cortexT] = CreateTemplates_2classes_Gaussian(header, wmSeg, cortexSeg, ...
                    sigmaWM, sigmaCortex, halfwidthWm, halfwidthCortex)

disp('create templates ...');

% create probability templates

cortexT = single(cortexSeg);
G_cortexT = Gaussianfilter(cortexT, header, sigmaCortex, sigmaCortex, sigmaCortex, halfwidthCortex, 'Real');
clear cortexSeg

wmT = single(wmSeg);
G_wmT = Gaussianfilter(wmT, header, sigmaWM, sigmaWM, sigmaWM, halfwidthWm, 'Real');
clear wmSeg

% clear csfSeg wmSeg1 wmSeg2 cortexSeg outlierSeg

wmflag = 0;
cortexflag = 0;

minP = 0.02;

for k = 1:header.zsize
    for j = 1:header.ysize
        for i = 1:header.xsize
            
            p_wm = G_wmT(j, i, k);
            p_cortex = G_cortexT(j, i, k);

            if ( G_wmT(j, i, k)>minP )
                wmflag = 1;
            end
            
            if ( G_cortexT(j, i, k)>minP )
                cortexflag = 1;
            end
                     
            if ( wmflag + cortexflag == 0 )
                continue;
            end
            
            wmT(j, i, k) = p_wm/(p_wm+p_cortex);
            cortexT(j, i, k) = p_cortex/(p_wm+p_cortex);
            
            wmflag = 0;
            cortexflag = 0;
        end
    end
end
return;
