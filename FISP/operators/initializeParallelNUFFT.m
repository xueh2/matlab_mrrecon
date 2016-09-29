function operator=initializeParallelNUFFT(data)
%This is a parallel a parallel implementation of the NUFFT operator usually
%created with initializeMultipleOperator
%it aims at being used with the parallel toolbox but is not efficient
%enough

show=0;




[xSize,ySize,coilSize]=size(data.csm);
imSize=[xSize,ySize];
timeSize=size(data.kcoor,3);






w = 1; phase = 1; shift = [0,0]; mode=1;
NUFFT_operator=cell(timeSize,1);
for i=1:timeSize
    ik=data.kcoor(:,1,i)/(2*pi)+1i*data.kcoor(:,2,i)/(2*pi);
    if show
        plot(real(ik),imag(ik),'.')
        pause(0.01)
    end
    if (max(abs(real(ik)))>0.50001)||(max(abs(imag(ik)))>0.50001)
        fprintf(1,'Error: kcoord needs to be normalized in [-pi,pi].\n')
    end
    NUFFT_operator{i}=NUFFT(ik,w,shift,imSize); %=NUFFT(k(indices),w,phase,shift,imSize,mode); %new version, not used (COMPATIBILITY)
end






index=1;
for i=1:data.YSize(5)
    for j=1:data.YSize(3)
        NUFFTs{index}=struct(NUFFT_operator{i});
        index=index+1;
    end
end

inputSize=[data.XSize(1) data.XSize(2) data.YSize(3) data.YSize(4) data.YSize(5) data.YSize(6)];
outputSize=[data.YSize(1) data.YSize(2) data.YSize(3) data.YSize(4) data.YSize(5) data.YSize(6)];

operator=SimpleOperator(@(x) processFct(x,NUFFTs,outputSize),@(x) processFct_t(x,NUFFTs,inputSize));

end

function y=processFct(x,NUFFTs,outputSize)

    y=zeros(outputSize);
    nLoops=outputSize(3)*outputSize(5);
    parfor i=1:nLoops
        y(:,:,i) = performNUFFT(NUFFTs{i},x(:,:,i));
    end

end


function y=processFct_t(x,NUFFTs,outputSize)

    %size of input
    y=zeros(outputSize);
    nLoops=outputSize(3)*outputSize(5);
    parfor i=1:nLoops
        y(:,:,i) = performNUFFT_t(NUFFTs{i},x(:,:,i));
    end

end


function res=performNUFFT_t(a,b)

    b = b(:).*a.w(:);
	res = nufft_adj(b, a.st)/sqrt(prod(a.imSize));
	res = reshape(res, a.imSize(1), a.imSize(2));

end


function res=performNUFFT(a,b)

    b = reshape(b,a.imSize(1),a.imSize(2));
	res = nufft(b, a.st)/sqrt(prod(a.imSize)).*a.w(:);
	res = reshape(res,a.dataSize(1),a.dataSize(2));

end

