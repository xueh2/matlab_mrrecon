function modelSeries = KLT_dynamic_model(a, widthForModel, ratio)
% compute KLT sliding window model series
% a : [RO E1 N]
% for evry widthForModel frames, the eigen image is computed
% this process slides over all frames

s = size(a);
N = s(3);

modelSeries = a;

if mod(widthForModel,2)==0
    halfw = floor(widthForModel/2)-1;
else
    halfw = floor(widthForModel/2);
end

for n=1:N
    frames = n-halfw : n+halfw;
    
    mf = min(frames(:));
    if ( mf <= 0 )
        frames = frames - mf + 1;
    end
    
    mf = max(frames(:));
    if ( mf > N )
        frames = frames - (mf-N);
    end
    
    im = a(:,:,frames(:));
    
%     [im_r, V, D] = KLT_Filter3D(im, ratio);
    
    [imE,V,D] = KL_Eigenimage(im);
    modelSeries(:,:,n) = imE(:,:,end);    
end