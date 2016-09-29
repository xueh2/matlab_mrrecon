function operator=initializeImageParallelMultipleOperator(operatorFct,transposeOperatorFct,inputSize,outputSize,fctIncrementationDimension)
%create an operator acting on a 6 dimensional matrix by repeating a same operation:
%[row][column][coil][slice][timepoint][repetition]
%operatorFct: direct operator handle function
%transposeOperatorFct: transpose operator handle function
%inputSize: size of the input of the direct function
%(operators on images will have inputSize=[row][column], the operator will
%repeat the process on all other dimensions)
%outputSize: size of the output of the direct function
%fctIncrementationDimension: if the function is a list of functions (for
%instance the operator is changing with time), specify on which dimension
%to iterate it (in the case of time, that would be 5)
%
%Example
%For a variable
%x=rand(30,30,5,3);
%the operator which will compute the 2D fourier transform on each slice is:
%F=initializeMultipleOperator(@(z) fftshift(fft2(z))/sqrt(numel(x)),@(z) sqrt(numel(x))*ifft2(fftshift(z)),[size(x,1),size(x,2)]);
%
%CAUTION: this version can only act on images (2D) and tries to use the
%parallel imaging toolbox (though it does not seem to improve the speed..)

%compute output size if not given
if nargin<4
    x=zeros(inputSize);
    outputSize=size(operatorFct(x));
end

%operator may vary with time or coil
if nargin<5
    %if argument is not specified, the function is the same for all
    %iterations
    fctIncrementationDimension=0;
end


%uniformise dimensions
N=length(inputSize);
if N<6
    inputSize=[inputSize,ones(1,6-N)];
end
N=length(outputSize);
if N<6
    outputSize=[outputSize,ones(1,6-N)];
end




operator=SimpleOperator(@(x) processFct(x,operatorFct,inputSize,outputSize,fctIncrementationDimension),@(x) processFct(x,transposeOperatorFct,outputSize,inputSize,fctIncrementationDimension));

end

function y=processFct(x,Fct,inputSize,outputSize,fctIncrementationDimension)

%size of input
dimX=size(x);
N=length(dimX);
if N<6
    dimX=[dimX,ones(1,6-N)];
end

[sx,sy,nLoops]=size(x);
clear sx sy;


%prepare y dimensions
dimY=outputSize;
i=(dimX>inputSize);
dimY(i)=dimX(i);
y=zeros(dimY);

if (fctIncrementationDimension~=0)
    increment=1:nLoops;
    increment=mod(floor((increment-1)/prod(dimX(3:fctIncrementationDimension-1))),dimX(fctIncrementationDimension))+1;
end

matlabpool
parfor i=1:nLoops

    %%compute and store
    if (fctIncrementationDimension==0)
        y(:,:,i) = Fct(x(:,:,i));
    else
        y(:,:,i) = Fct(increment(i),x(:,:,i));
    end


end
matlabpool close






end
