function opediff = buildDiffOperator(dimension)
%builds the operator (and its transpose) computing the volume of
%differences between 2 consecutive slices along a given dimension
%it uses the class SimpleOperator
if nargin<1
	dimension=5;
end

	if dimension==1
		imdiff=@(x) x(2:end,:,:,:,:,:)-x(1:(end-1),:,:,:,:,:);
		opediff=SimpleOperator(@(x) imdiff(x),@(x) imdiffT1(x));
	else if dimension==2
		imdiff=@(x) x(:,2:end,:,:,:,:)-x(:,1:(end-1),:,:,:,:);
		opediff=SimpleOperator(@(x) imdiff(x),@(x) imdiffT2(x));
	else if dimension==3
		imdiff=@(x) x(:,:,2:end,:,:,:)-x(:,:,1:(end-1),:,:,:);
		opediff=SimpleOperator(@(x) imdiff(x),@(x) imdiffT3(x));
	else if dimension==4
		imdiff=@(x) x(:,:,:,2:end,:,:)-x(:,:,:,1:(end-1),:,:);
		opediff=SimpleOperator(@(x) imdiff(x),@(x) imdiffT4(x));
	else if dimension==5
		imdiff=@(x) x(:,:,:,:,2:end,:)-x(:,:,:,:,1:(end-1),:);
		opediff=SimpleOperator(@(x) imdiff(x),@(x) imdiffT5(x));
	else if dimension==6
		imdiff=@(x) x(:,:,:,:,:,2:end)-x(:,:,:,:,:,1:(end-1));
		opediff=SimpleOperator(@(x) imdiff(x),@(x) imdiffT6(x));
	else
		fprintf('Wrong input dimension.\n')
		opediff= 0;
    end
    end
    end
    end
    end
    end




end

function y=imdiffT1(x)
	[sx,xy,c,sz,t,r]=size(x);
	y=zeros(sx+1,xy,c,sz,t,r);
	y(1:end-1,:,:,:,:,:)=-x;
	y(2:end,:,:,:,:,:)=y(2:end,:,:,:,:,:)+x;
end

function y=imdiffT2(x)
	[sx,xy,c,sz,t,r]=size(x);
	y=zeros(sx,xy+1,c,sz,t,r);
	y(:,1:end-1,:,:,:,:)=-x;
	y(:,2:end,:,:,:,:)=y(:,2:end,:,:,:,:)+x;
end

function y=imdiffT3(x)
	[sx,xy,c,sz,t,r]=size(x);
	y=zeros(sx,xy,c+1,sz,t,r);
	y(:,:,1:end-1,:,:,:)=-x;
	y(:,:,2:end,:,:,:)=y(:,:,2:end,:,:,:)+x;
end

function y=imdiffT4(x)
	[sx,xy,c,sz,t,r]=size(x);
	y=zeros(sx,xy,c,sz+1,t,r);
	y(:,:,:,1:end-1,:,:)=-x;
	y(:,:,:,2:end,:,:)=y(:,:,:,2:end,:,:)+x;
end

function y=imdiffT5(x)
	[sx,xy,c,sz,t,r]=size(x);
	y=zeros(sx,xy,c,sz,t+1,r);
	y(:,:,:,:,1:end-1,:)=-x;
	y(:,:,:,:,2:end,:)=y(:,:,:,:,2:end,:)+x;
end

function y=imdiffT6(x)
	[sx,xy,c,sz,t,r]=size(x);
	y=zeros(sx,xy,c,sz,t,r+1);
	y(:,:,:,:,:,1:end-1)=-x;
	y(:,:,:,:,:,2:end)=y(:,:,:,:,:,2:end)+x;
end
