function sampledLines = createSamplingParttern_RegularStepSizePELines(Npe, R, numOfCentralLines)
% generate the kspace sampled lines with reduction factor R and central lines

numOfLines = ceil(Npe/R);

outsideLines = numOfLines - numOfCentralLines;

N = outsideLines/2;
outsideRange = Npe/2-numOfCentralLines/2;

stepSize = outsideRange/N;

sampledLines = zeros([Npe 1]);
sampledLines(Npe/2-numOfCentralLines/2 : Npe/2+numOfCentralLines/2-1) = 1;

for k=1:N    
    n = k*stepSize;
    
    l = floor(Npe/2-numOfCentralLines/2-n);
    if ( l < 1 )
        l = 1;
    end
    sampledLines(l) = 1;
    
    l = floor(Npe/2+numOfCentralLines/2-1+n);
    if ( l > Npe )
        l = Npe;
    end
    sampledLines(l) = 1;
end
