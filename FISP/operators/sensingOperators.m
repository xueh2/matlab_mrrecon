function [nuFT,coilOpe]=sensingOperators(kcoor,csm,show)
%creates the NUFFT operator for the coordinates kcoor [x,y,t]
%also creates the csm operator with csm [x,y,coil]
%they can be directly applied on data x of type [x y coil z t repetition]
%by calling nuFT*x or coilOpe*x or nuFT*coilOpe*x
%the transpose of those operators are also defined.

if nargin<2
    show=0;
end

%sizes
[xSize,ySize,coilSize]=size(csm);
imSize=[xSize,ySize];
timeSize=size(kcoor,3);



%% Coil operator
%vidshow(csm)
coilOpe=initializeMultipleOperator(@(x) repmat(x,[1,1,coilSize]).*csm, @(x) sum(x.*conj(csm),3),imSize,[imSize,coilSize]);



%% NUFFT operators
%NUFFT parameters
w = 1; phase = 1; shift = [0,0]; mode=1;
NUFFT_operator=cell(timeSize,1);
for i=1:timeSize
    ik=kcoor(:,1,i)/(2*pi)+1i*kcoor(:,2,i)/(2*pi);
    if show
        plot(real(ik),imag(ik),'.')
        pause(0.01)
    end
    if (max(abs(real(ik)))>0.50001)||(max(abs(imag(ik)))>0.50001)
        fprintf(1,'Error: kcoord needs to be normalized in [-pi,pi].\n')
    end
    NUFFT_operator{i}=NUFFT(ik,w,shift,imSize); %=NUFFT(k(indices),w,phase,shift,imSize,mode); %new version, not used (COMPATIBILITY)
end

%% operator on all data
if exist('matlabpool') %implementation with parfor: needs the adequat toolbox
    nuFT=initializeImageParallelMultipleOperator(@(t,x) NUFFT_operator{t}*x, @(t,x) NUFFT_operator{t}'*x,imSize,size(NUFFT_operator{1}*zeros(imSize)),5);
	fprintf(1,'Please be aware that the NUFFT multiple operator was initialized using Parallel Computing toolbox.\n')
else
    nuFT=initializeMultipleOperator(@(t,x) NUFFT_operator{t}*x, @(t,x) NUFFT_operator{t}'*x,imSize,size(NUFFT_operator{1}*zeros(imSize)),5);
end

%% close window if open

    if show
        close
        pause(0.01)
    end

end
