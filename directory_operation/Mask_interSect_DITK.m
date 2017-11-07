
function newMask = Mask_interSect_DITK(targetMask, headerTarget, mask, headerMask)
% intersect two masks in the DITK center coordinates

[t_ysize, t_xsize, t_zsize] = size(targetMask);
t_xvoxelsize = headerTarget.xvoxelsize;
t_yvoxelsize = headerTarget.yvoxelsize;
t_zvoxelsize = headerTarget.zvoxelsize;

[ysize, xsize, zsize] = size(mask);
xvoxelsize = headerMask.xvoxelsize;
yvoxelsize = headerMask.yvoxelsize;
zvoxelsize = headerMask.zvoxelsize;

newMask = false(size(targetMask));

head = 1;
ender = floor(t_zsize/4);
oldindex = 0;
points = [];
for tt = 1:4
    index = find(targetMask(:,:,head:ender)>0);
%     index = single(index);
    [pointy, pointx, pointz] = ind2sub(size(targetMask(:,:,head:ender)), index);
    pointy = uint16(pointy);
    pointx = uint16(pointx);
    pointz = uint16(pointz);
    points = [points; pointx pointy pointz+head-1 ];
    oldindex = index(end);
    clear pointy pointx pointz index
    head = ender + 1;
    ender = floor((tt+1)*t_zsize/4);
end

[numOld, ndim] = size(points);
head = 1;
ender = floor(numOld/4);
for tt = 1:4
    tt
    pointinworld = image2world_DITK(points(head:ender, :), t_xvoxelsize, t_yvoxelsize, t_zvoxelsize,...
                                    t_xsize, t_ysize, t_zsize);

    pointinImage = world2image_DITK(pointinworld, xvoxelsize, yvoxelsize, zvoxelsize,...
                                    xsize, ysize, zsize);
    clear pointinworld                            
    pointinImage = round(pointinImage);

    [num, ndim] = size(pointinImage);

    for k = 1:num
        
        tx = pointinImage(k, 1);
        ty = pointinImage(k, 2);
        tz = pointinImage(k, 3);
        
        if ( (tx<1) | (tx>xsize) |...
                (ty<1) | (ty>ysize) |...
                (tz<1) | (tz>zsize) )
            continue;
        end
        x = points(k+head-1, 1);
        y = points(k+head-1, 2);
        z = points(k+head-1, 3);
        
        if ( (targetMask(y, x, z)==true) & (mask(ty, tx, tz)==true) )
            newMask(y, x, z) = 1;
        end

    end
    head = ender + 1;
    ender = floor((tt+1)*numOld/4);
end
return