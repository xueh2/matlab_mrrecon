function imagescnVolumes(v1,v2,v3,v4,v5,v6,v7)

v1 = flipdim(v1, 1);
v1 = flipdim(v1, 2);

v = zeros([size(v1) nargin]);
v(:,:,1) = v1;

v2 = flipdim(v2, 1);
v2 = flipdim(v2, 2);
v(:,:,2) = v2;

if ( nargin > 2 )
    v3 = flipdim(v3, 1);
    v3 = flipdim(v3, 2);
    v(:,:,3) = v3;
end

if ( nargin > 3 )
    v4 = flipdim(v4, 1);
    v4 = flipdim(v4, 2);
    v(:,:,4) = v4;
end

if ( nargin > 4 )
    v5 = flipdim(v5, 1);
    v5 = flipdim(v5, 2);    
    v(:,:,5) = v5;
end

if ( nargin > 5 )
    v6 = flipdim(v6, 1);
    v6 = flipdim(v6, 2);
    v(:,:,6) = v6;
end

if ( nargin > 6 )
    v7 = flipdim(v7, 1);
    v7 = flipdim(v7, 2);
    v(:,:,7) = v7;
end

imagescn(v, [], [1 nargin], []);
