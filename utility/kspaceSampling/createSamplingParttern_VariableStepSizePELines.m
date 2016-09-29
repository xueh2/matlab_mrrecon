function sampledLines = createSamplingParttern_VariableStepSizePELines(Npe, R, numOfCentralLines)
% generate the kspace sampled lines with reduction factor R and central lines

numOfLines = ceil(Npe/R);

outsideLines = numOfLines - numOfCentralLines;

N = outsideLines/2;
outsideRange = Npe/2-numOfCentralLines/2;

stepSize = 1.1;

while ( (1-stepSize^N)/(1-stepSize) < outsideRange )
    stepSize = stepSize + 0.1;
end

sampledLines = zeros([Npe 1]);
sampledLines(Npe/2-numOfCentralLines/2 : Npe/2+numOfCentralLines/2-1) = 1;

for k=1:N    
    n = (1-stepSize^k)/(1-stepSize);
    
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
