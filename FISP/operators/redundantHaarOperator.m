%creates an operator (and its transpose/pseudoinverse) which computes low frequency and high frequency haar redundant wavelets for a given dimension
%the algorithm is the one implemented for haar 2-level redundant wavelets on the Rice university library
%Example:
%to create a 2d redundant wavelet operator, type wavelet=redundantHaarOperator(1)*redundantHaarOperator(2)
%for an image X, WX=wavelet*X
%WX will have its dimensions multiplicated by 2 compared to X, first half of one dimension containing low-pass(es) on this dimension, second half high pass(es)
function ope = redundantHaarOperator(dimension)

		


if nargin<1
	dimension=5;
end



	if dimension==1
		ope=SimpleOperator(@(x) directOp1(x),@(x) inverseOp1(x));
	else if dimension==2
		ope=SimpleOperator(@(x) directOp2(x),@(x) inverseOp2(x));
	else if dimension==3
		ope=SimpleOperator(@(x) directOp3(x),@(x) inverseOp3(x));
	else if dimension==4
		ope=SimpleOperator(@(x) directOp4(x),@(x) inverseOp4(x));
	else if dimension==5
		ope=SimpleOperator(@(x) directOp5(x),@(x) inverseOp5(x));
	else if dimension==6
		ope=SimpleOperator(@(x) directOp6(x),@(x) inverseOp6(x));
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

%% dimension 1

function y=directOp1(x)
    t=size(x,1);
    m=(x+x([2:t,1],:,:,:,:,:))/2;
    y=cat(1,m,(x-m));
end

function y=inverseOp1(x)
	t=size(x,1);
	y=(x(1:(t/2),:,:,:,:,:)+x((t/2)+1:end,:,:,:,:,:)...
        +x([t/2,1:(t/2)-1],:,:,:,:,:)-x([t,(t/2)+1:t-1],:,:,:,:,:))/2;

end

%% dimension 2

function y=directOp2(x)
    t=size(x,2);
    m=(x+x(:,[2:t,1],:,:,:,:))/2;
    y=cat(2,m,(x-m));
end

function y=inverseOp2(x)
	t=size(x,2);
	y=(x(:,1:(t/2),:,:,:,:)+x(:,(t/2)+1:end,:,:,:,:)+...
            x(:,[t/2,1:(t/2)-1],:,:,:,:)-x(:,[t,(t/2)+1:t-1],:,:,:,:))/2;

end

%% dimension 3

function y=directOp3(x)
    t=size(x,3);
    m=(x+x(:,:,[2:t,1],:,:,:))/2;
    y=cat(3,m,(x-m));
end

function y=inverseOp3(x)
	t=size(x,3);
	y=(x(:,:,1:(t/2),:,:,:)+x(:,:,(t/2)+1:end,:,:,:)+...
            x(:,:,[t/2,1:(t/2)-1],:,:,:)-x(:,:,[t,(t/2)+1:t-1],:,:,:))/2;

end

%% dimension 4

function y=directOp4(x)
    t=size(x,4);
    m=(x+x(:,:,:,[2:t,1],:,:))/2;
    y=cat(4,m,(x-m));
end

function y=inverseOp4(x)
	t=size(x,4);
	y=(x(:,:,:,1:(t/2),:,:)+x(:,:,:,(t/2)+1:end,:,:)+...
            x(:,:,:,[t/2,1:(t/2)-1],:,:)-x(:,:,:,[t,(t/2)+1:t-1],:,:))/2;

end

%% dimension 5

function y=directOp5(x)
    t=size(x,5);
    m=(x+x(:,:,:,:,[2:t,1],:))/2;
    y=cat(5,m,(x-m));
end

function y=inverseOp5(x)
	t=size(x,5);
	y=(x(:,:,:,:,1:(t/2),:)+x(:,:,:,:,(t/2)+1:end,:)+...
            x(:,:,:,:,[t/2,1:(t/2)-1],:)-x(:,:,:,:,[t,(t/2)+1:t-1],:))/2;

end

%% dimension 6

function y=directOp6(x)
    t=size(x,6);
    m=(x+x(:,:,:,:,:,[2:t,1]))/2;
    y=cat(6,m,(x-m));
end

function y=inverseOp6(x)
	t=size(x,6);
	y=(x(:,:,:,:,:,1:(t/2))+x(:,:,:,:,:,(t/2)+1:end)+...
            x(:,:,:,:,:,[t/2,1:(t/2)-1])-x(:,:,:,:,:,[t,(t/2)+1:t-1]))/2;

end
