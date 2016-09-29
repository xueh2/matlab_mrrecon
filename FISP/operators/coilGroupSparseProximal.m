function y = coilGroupSparseProximal(x,lambdas)
    %implements W'*softtresh( W*x,Lambdas) where the
    %groups are along the coil dimension
    y=zeros(size(x));
    size_w=size(x).*((size(x)>1)+1);
    size_w(3)=size(x,3);
    all_w=zeros(size_w);
    size_w(3)=1;
    group_w=zeros(size_w);
    for c=1:size(x,3)
        %we need to repeat the wavelet computation for each coil dimension,
        %because if we directly give x to the current implementation of
        %redundantHaar6D, it will compute wavelets along the coil dimension
        %(which does not make any sense)
       all_w(:,:,c,:,:,:)=redundantHaar6D(x(:,:,c,:,:,:),1);
       group_w=group_w+all_w(:,:,c,:,:,:).*conj(all_w(:,:,c,:,:,:));
    end
    group_w=sqrt(group_w);
	tresh_w=wavSoftTreshold(group_w,lambdas);
    ind=tresh_w~=0;
    tresh_w(ind)=tresh_w(ind)./group_w(ind);
    
    for c=1:size(x,3)
       y(:,:,c,:,:,:)=redundantHaar6D(all_w(:,:,c,:,:,:).*tresh_w,-1);
    end
end
