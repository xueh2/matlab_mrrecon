function I = apply_normorientation(I, Rotate90, Flip);
% function I = apply_normorientation(I, Rotate90, Flip);
% 

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

Isize=size(I);
tmp=I(:,:,:);
% rotate then flip
if ndims(tmp)>2
    for i=1:size(tmp,3)
        tmp2(:,:,i)=rot90(tmp(:,:,i), Rotate90);
    end   
    tmp=tmp2;
else
    tmp=rot90(tmp,Rotate90);
end
switch Flip
    case 0 % do nothing
    case 1
        tmp=tmp(end:-1:1,:,:);
    case 2
        tmp=tmp(:,end:-1:1,:);
end
I = reshape(tmp,[size(tmp,1) size(tmp,2) Isize(3:end)]);


return

   
