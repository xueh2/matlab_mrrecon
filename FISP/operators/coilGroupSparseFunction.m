function val = coilGroupSparseFunction(x,Lambdas)
    %implements ||Lambdas \otimes W*x||_{2,1} (group sparsity) where the
    %groups are along the coil dimension
    size_w=size(x).*((size(x)>1)+1);
    size_w(3)=1;
    w=zeros(size_w);
    for c=1:size(x,3)
        %we need to repeat the wavelet computation for each coil dimension,
        %because if we directly give x to the current implementation of
        %redundantHaar6D, it will compute wavelets along the coil dimension
        %(which does not make any sense)
       w=w+abs(redundantHaar6D(x(:,:,c,:,:,:),1)).^2;
    end
	val=wavL1norm(sqrt(w),Lambdas);
end
