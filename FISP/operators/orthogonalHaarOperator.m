% Orthogonal Haar Wavelet Transform for the first 6 dimensions

% Qiu Wang
% 07/18/2012

function ope = orthogonalHaarOperator(dimension)

		


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
    m=(x(1:2:end-1,:,:,:,:,:)+x(2:2:end,:,:,:,:,:))/2;
    m2 = (x(1:2:end-1,:,:,:,:,:)-x(2:2:end,:,:,:,:,:))/2;
    y=cat(1,m,m2);
end

function y=inverseOp1(x)
	t=size(x,1);
    y1 = x(1:(t/2),:,:,:,:,:)+x((t/2)+1:end,:,:,:,:,:);
    y2 = x(1:(t/2),:,:,:,:,:)-x((t/2)+1:end,:,:,:,:,:);
    y = zeros(size(x));
    y(1:2:end-1,:,:,:,:,:) = y1;
    y(2:2:end,:,:,:,:,:)=y2;
end

%% dimension 2

function y=directOp2(x)
    t=size(x,1);
    m=(x(:,1:2:end-1,:,:,:,:)+x(:,2:2:end,:,:,:,:))/2;
    m2 = (x(:,1:2:end-1,:,:,:,:)-x(:,2:2:end,:,:,:,:))/2;
    y=cat(2,m,m2);
end

function y=inverseOp2(x)
	t=size(x,1);
    y1 = x(:,1:(t/2),:,:,:,:)+x(:,(t/2)+1:end,:,:,:,:);
    y2 = x(:,1:(t/2),:,:,:,:)-x(:,(t/2)+1:end,:,:,:,:);
    y = zeros(size(x));
    y(:,1:2:end-1,:,:,:,:) = y1;
    y(:,2:2:end,:,:,:,:)=y2;
end


%% dimension 3

function y=directOp3(x)
    t=size(x,1);
    m=(x(:,:,1:2:end-1,:,:,:)+x(:,:,2:2:end,:,:,:))/2;
    m2 = (x(:,:,1:2:end-1,:,:,:)-x(:,:,2:2:end,:,:,:))/2;
    y=cat(3,m,m2);
end

function y=inverseOp3(x)
	t=size(x,1);
    y1 = x(:,:,1:(t/2),:,:,:)+x(:,:,(t/2)+1:end,:,:,:);
    y2 = x(:,:,1:(t/2),:,:,:)-x(:,:,(t/2)+1:end,:,:,:);
    y = zeros(size(x));
    y(:,:,1:2:end-1,:,:,:) = y1;
    y(:,:,2:2:end,:,:,:)=y2;
end


%% dimension 4

function y=directOp4(x)
    t=size(x,1);
    m=(x(:,:,:,1:2:end-1,:,:)+x(:,:,:,2:2:end,:,:))/2;
    m2 = (x(:,:,:,1:2:end-1,:,:)-x(:,:,:,2:2:end,:,:))/2;
    y=cat(4,m,m2);
end

function y=inverseOp4(x)
	t=size(x,1);
    y1 = x(:,:,:,1:(t/2),:,:)+x(:,:,:,(t/2)+1:end,:,:);
    y2 = x(:,:,:,1:(t/2),:,:)-x(:,:,:,(t/2)+1:end,:,:);
    y = zeros(size(x));
    y(:,:,:,1:2:end-1,:,:) = y1;
    y(:,:,:,2:2:end,:,:)=y2;
end

%% dimension 5

function y=directOp5(x)
    t=size(x,1);
    m=(x(:,:,:,:,1:2:end-1,:)+x(:,:,:,:,2:2:end,:))/2;
    m2 = (x(:,:,:,:,1:2:end-1,:)-x(:,:,:,:,2:2:end,:))/2;
    y=cat(5,m,m2);
end

function y=inverseOp5(x)
	t=size(x,1);
    y1 = x(:,:,:,:,1:(t/2),:)+x(:,:,:,:,(t/2)+1:end,:);
    y2 = x(:,:,:,:,1:(t/2),:)-x(:,:,:,:,(t/2)+1:end,:);
    y = zeros(size(x));
    y(:,:,:,:,1:2:end-1,:) = y1;
    y(:,:,:,:,2:2:end,:)=y2;
end

%% dimension 6

function y=directOp6(x)
    t=size(x,1);
    m=(x(:,:,:,:,:,1:2:end-1)+x(:,:,:,:,:,2:2:end))/2;
    m2 = (x(:,:,:,:,:,1:2:end-1)-x(:,:,:,:,:,2:2:end))/2;
    y=cat(6,m,m2);
end

function y=inverseOp6(x)
	t=size(x,1);
    y1 = x(:,:,:,:,:,1:(t/2))+x(:,:,:,:,:,(t/2)+1:end);
    y2 = x(:,:,:,:,:,1:(t/2))-x(:,:,:,:,:,(t/2)+1:end);
    y = zeros(size(x));
    y(:,:,:,:,:,1:2:end-1) = y1;
    y(:,:,:,:,:,2:2:end)=y2;
end
