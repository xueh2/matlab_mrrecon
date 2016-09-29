%scales the high frequency part of the redundant Haar wavelets for a given dimension (from 1 to 6)
%sizeX: size of the wavelet image/volume
%dimension: dimension on which to apply the scale
%scale: multiplication factor for high frequency components
function ope = scaleHaarWav(sizeX,dimension,scale)

		


if nargin<1
	dimension=5;
end
if nargin<2
	scale=2;
end
t=sizeX(dimension);

	if dimension==1
		ope=SimpleOperator(@(x) cat(1,x(1:(t/2),:,:,:,:,:),scale*x((t/2)+1:t,:,:,:,:,:)),@(x) cat(1,x(1:(t/2),:,:,:,:,:),(1/scale)*x((t/2)+1:t,:,:,:,:,:)));
	else if dimension==2
		ope=SimpleOperator(@(x) cat(2,x(:,1:(t/2),:,:,:,:),scale*x(:,(t/2)+1:t,:,:,:,:)),@(x) cat(2,x(:,1:(t/2),:,:,:,:),(1/scale)*x(:,(t/2)+1:t,:,:,:,:)));
	else if dimension==3
		ope=SimpleOperator(@(x) cat(3,x(:,:,1:(t/2),:,:,:),scale*x(:,:,(t/2)+1:t,:,:,:)),@(x) cat(3,x(:,:,1:(t/2),:,:,:),(1/scale)*x(:,:,(t/2)+1:t,:,:,:)));
	else if dimension==4
		ope=SimpleOperator(@(x) cat(4,x(:,:,:,1:(t/2),:,:),scale*x(:,:,:,(t/2)+1:t,:,:)),@(x) cat(4,x(:,:,:,1:(t/2),:,:),(1/scale)*x(:,:,:,(t/2)+1:t,:,:)));
	else if dimension==5
		ope=SimpleOperator(@(x) cat(5,x(:,:,:,:,1:(t/2),:),scale*x(:,:,:,:,(t/2)+1:t,:)),@(x) cat(5,x(:,:,:,:,1:(t/2),:),(1/scale)*x(:,:,:,:,(t/2)+1:t,:)));
	else if dimension==6
		ope=SimpleOperator(@(x) cat(6,x(:,:,:,:,:,1:(t/2)),scale*x(:,:,:,:,:,(t/2)+1:t)),@(x) cat(6,x(:,:,:,:,:,1:(t/2)),(1/scale)*x(:,:,:,:,:,(t/2)+1:t)));
	else
		fprintf('Wrong input dimension.\n')
		ope= 0;
    end
    end
    end
    end
    end
    end


end
