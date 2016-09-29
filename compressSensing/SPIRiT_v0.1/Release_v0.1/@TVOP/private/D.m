function res = D(image)

%
% res = D(image)
%
% image = a 2D multi-coil image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005
% modified by Hui Xue

Dx = image([2:end,end],:, :, :) - image;
Dy = image(:,[2:end,end], :, :) - image;

if ( numel(size(image)) == 3 )
    res = cat(4,Dx,Dy);
end

if ( numel(size(image)) == 4 )
    res = cat(5,Dx,Dy);
end

