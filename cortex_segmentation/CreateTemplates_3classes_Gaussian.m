
function [wmT1, wmT2, cortexT] = CreateTemplates_3classes_Gaussian(header, wmSeg1, wmSeg2, cortexSeg, ...
                    sigmaWM1, sigmaWM2, sigmaCortex, halfwidthWm1, halfwidthWm2, halfwidthCortex)

disp('create templates ...');


% create probability templates

cortexT = single(cortexSeg);
G_cortexT = Gaussianfilter(cortexT, header, sigmaCortex, sigmaCortex, sigmaCortex, halfwidthCortex, 'Real');
clear cortexSeg

wmT1 = single(wmSeg1);
G_wmT1 = Gaussianfilter(wmT1, header, sigmaWM1, sigmaWM1, sigmaWM1, halfwidthWm1, 'Real');
clear wmSeg1

wmT2 = single(wmSeg2);
G_wmT2 = Gaussianfilter(wmT2, header, sigmaWM2, sigmaWM2, sigmaWM2, halfwidthWm2, 'Real');
clear wmSeg2

% clear csfSeg wmSeg1 wmSeg2 cortexSeg outlierSeg

wm1flag = 0;
wm2flag = 0;
cortexflag = 0;

minP = 0.02;

for k = 1:header.zsize
    for j = 1:header.ysize
        for i = 1:header.xsize
            
            p_wm1 = G_wmT1(j, i, k);
            p_wm2 = G_wmT2(j, i, k);
            p_cortex = G_cortexT(j, i, k);

            if ( G_wmT1(j, i, k)>minP )
                wm1flag = 1;
            end
            
            if ( G_wmT2(j, i, k)>minP )
                wm2flag = 1;
            end   
            
            if ( G_cortexT(j, i, k)>minP )
                cortexflag = 1;
            end
                     
            if ( wm1flag + wm2flag + cortexflag == 0 )
                continue;
            end
            
            wmT1(j, i, k) = p_wm1/(p_wm1+p_wm2+p_cortex);
            wmT2(j, i, k) = p_wm2/(p_wm1+p_wm2+p_cortex);
            cortexT(j, i, k) = p_cortex/(p_wm1+p_wm2+p_cortex);
            
            wm1flag = 0;
            wm2flag = 0;
            cortexflag = 0;
        end
    end
end
return;
