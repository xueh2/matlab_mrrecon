function operator=initializeMultipleOperator(operatorFct,transposeOperatorFct,inputSize,outputSize,fctIncrementationDimension)
%operator acting on a 6 dimensional matrix by repeating a same operation:
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

%number of loops
loops(1)=1*(inputSize(1)>1)+dimX(1)*(inputSize(1)==1);
loops(2)=1*(inputSize(2)>1)+dimX(2)*(inputSize(2)==1);
loops(3)=1*(inputSize(3)>1)+dimX(3)*(inputSize(3)==1);
loops(4)=1*(inputSize(4)>1)+dimX(4)*(inputSize(4)==1);
loops(5)=1*(inputSize(5)>1)+dimX(5)*(inputSize(5)==1);
loops(6)=1*(inputSize(6)>1)+dimX(6)*(inputSize(6)==1);


%prepare y dimensions
dimY=outputSize;
i=(dimX>inputSize);
dimY(i)=dimX(i);
y=zeros(dimY);

%order of the loop (if there is only one loop, perform it first)
[temp,o]=sort(inputSize,'descend');

for i1=1:loops(o(1))
    if loops(o(1))==1
        inputIndices{o(1)}=1:inputSize(o(1));
        outputIndices{o(1)}=1:outputSize(o(1));
    else
        inputIndices{o(1)}=i1;
        outputIndices{o(1)}=i1;
    end
for i2=1:loops(o(2))
    if loops(o(2))==1
        inputIndices{o(2)}=1:inputSize(o(2));
        outputIndices{o(2)}=1:outputSize(o(2));
    else
        inputIndices{o(2)}=i2;
        outputIndices{o(2)}=i2;
    end
for i3=1:loops(o(3))
    if loops(o(3))==1
        inputIndices{o(3)}=1:inputSize(o(3));
        outputIndices{o(3)}=1:outputSize(o(3));
    else
        inputIndices{o(3)}=i3;
        outputIndices{o(3)}=i3;
    end
for i4=1:loops(o(4))
    if loops(o(4))==1
        inputIndices{o(4)}=1:inputSize(o(4));
        outputIndices{o(4)}=1:outputSize(o(4));
    else
        inputIndices{o(4)}=i4;
        outputIndices{o(4)}=i4;
    end
for i5=1:loops(o(5))
    if loops(o(5))==1
        inputIndices{o(5)}=1:inputSize(o(5));
        outputIndices{o(5)}=1:outputSize(o(5));
    else
        inputIndices{o(5)}=i5;
        outputIndices{o(5)}=i5;
    end
for i6=1:loops(o(6))
    if loops(o(6))==1
        inputIndices{o(6)}=1:inputSize(o(6));
        outputIndices{o(6)}=1:outputSize(o(6));
    else
        inputIndices{o(6)}=i6;
        outputIndices{o(6)}=i6;
    end
    
    %dimX
    %dimY
    
%         fprintf(1,'compute\n')
%          %Fct(inputIndices{fctIncrementationDimension},x(inputIndices{1},inputIndices{2},inputIndices{3},inputIndices{4},inputIndices{5},inputIndices{6}));
%          n1=inputIndices{1}
%          n2=inputIndices{2}
%          n3=inputIndices{3}
%          n4=inputIndices{4}
%          n5=inputIndices{5}
%          n6=inputIndices{6}
%          fprintf(1,'store\n')
%          o1=outputIndices{1}
%          o2=outputIndices{2}
%          o3=outputIndices{3}
%          o4=outputIndices{4}
%          o5=outputIndices{5}
%          o6=outputIndices{6}
         
         
    
    %%compute and store
    if (fctIncrementationDimension==0)
        y(outputIndices{1},outputIndices{2},outputIndices{3},outputIndices{4},outputIndices{5},outputIndices{6})=...
        Fct(x(inputIndices{1},inputIndices{2},inputIndices{3},inputIndices{4},inputIndices{5},inputIndices{6}));
    else
        y(outputIndices{1},outputIndices{2},outputIndices{3},outputIndices{4},outputIndices{5},outputIndices{6})=...
        Fct(inputIndices{fctIncrementationDimension},x(inputIndices{1},inputIndices{2},inputIndices{3},inputIndices{4},inputIndices{5},inputIndices{6}));
    end

        
        
end
end
end
end
end
end





end
