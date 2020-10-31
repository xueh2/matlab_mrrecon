
function PointPropagation_Only(filename, dx_file, dy_file, LandMarksfiles, pts_Properties, NumOfPts, sizeRatio, slicenumflag, aviName)
% visualize the propagated points using the deformation field

if ( isempty(dir(aviName))==0 )
    disp(['Already exsiting: ' aviName]);
%     return;
end

disp(['PointPropagation_Only :' aviName]);
% ---------------------------------------------------------------
% read images and set the parameters properly

[data, header] = LoadAnalyze(filename, 'Grey');

minD = min(data(:));
maxD = max(data(:));

if ( maxD > 255 )
    data = NormalizeImage(double(data)) .* 255;
end

[Dx, header] = LoadAnalyze(dx_file, 'Real');
[Dy, header] = LoadAnalyze(dy_file, 'Real');

[H, W, Num] = size(data);

% ---------------------------------------------------------------
% LandMarks

pts_Images = cell(NumOfPts, 1);
pts_Deformeds = cell(NumOfPts, 1);

pt_lineColor = cell(NumOfPts, 1);
pt_color = cell(NumOfPts, 1);
target_pt_color = cell(NumOfPts, 1);

for pp=1:NumOfPts

    pts_Property = pts_Properties{pp};
    LandMarksfile = LandMarksfiles{pp};
    
    pt_size(pp) = pts_Property.pt_size;
    pt_lineColor{pp} = pts_Property.pt_lineColor;
    pt_lineWidth(pp) = pts_Property.pt_lineWidth;
    lineflag(pp) = pts_Property.lineflag;
    
    pt_color{pp} = pts_Property.color;
    pt_symbol(pp) = pts_Property.symbol;
    pt_flag(pp) = pts_Property.pt_flag;
    
    target_pt_color{pp} = pts_Property.target_color;
    target_pt_symbol(pp) = pts_Property.target_symbol;
    target_flag(pp) = pts_Property.target_flag;
    
    if ( isempty(dir(LandMarksfile))==0 )
        [pts, numCell, cells] = VTKFile2PointCell(LandMarksfile);    
        pts_Image = world2image_DITK2(double(pts(:,2:4)), header);
 
        NumofPts = size(pts, 1);

        pts_Image = pts_Image(:, 1:2);

        pts_Deformed = UseDeformation(pts_Image, Dx, Dy, header);  
        % ---------------------------------------------------------------
        % size Ratio
        C = imresize(data(:, :, 1), sizeRatio, 'bicubic');

        [newH, newW] = size(C);

        pts_Image = pts_Image .* sizeRatio;
        pts_Image(:, 2) = newH - pts_Image(:, 2);

        pts_Deformed = pts_Deformed .* sizeRatio;

        for i=1:Num
            pts_Deformed(:, 2, i) = newH - pts_Deformed(:, 2, i);
        end

        pts_Images{pp} = pts_Image;
        pts_Deformeds{pp} = pts_Deformed;    
    else
        pts_Images{pp} = [];
        pts_Deformeds{pp} = [];    
    end
end
% ---------------------------------------------------------------

f = figure;
global imageRendered;
imageRendered = 0;
for kk=1:header.zsize
    
    kk
    slice = data(:, :, kk);

    slice = imresize(slice, sizeRatio, 'bicubic');
    slice = flipdim(slice,1);
    
    % render the contour and slice
    if ( slicenumflag )
        slicenum = kk;
    else
        slicenum = 0;
    end
    
    for pp=1:NumOfPts

        pts_Image = pts_Images{pp};
        pts_Deformed = pts_Deformeds{pp};  
        
        % ---------------------------------------
        if ( isempty(pts_Deformed) )
            continue;
        end
        
        current_pts_Deformed = pts_Deformed(:, :, kk);
        f = RenderDeformationLines2(slice, pts_Image(:,1:2), current_pts_Deformed, pt_lineWidth(pp), pt_lineColor{pp}, lineflag(pp), ...
            pt_size(pp), pt_color{pp},  pt_symbol(pp), pt_flag(pp), ...
            target_pt_color{pp},  target_pt_symbol(pp), target_flag(pp), ...
            slicenum, f);
    end 

    h = get(f, 'CurrentAxes');
    M(kk) = getframe(h);
    
    cla(h);
    imageRendered = 0;
end    

movie2avi(M, aviName, 'Compression', 'None', 'fps', 8);
close(f);

function  pts_Deformed = UseDeformation(pts_Image, Dx, Dy, header)
% use the deformation

num = size(pts_Image, 1);

numSlice = header.zsize;

pts_Deformed = zeros(num, 2, numSlice); % px, py, numofslice

Dx_2D = Dx(:, :, 1);
Dy_2D = Dy(:, :, 1);

for kk=1:numSlice
    % for every slice
    Dx_2D = Dx(:, :, kk);
    Dy_2D = Dy(:, :, kk);
    
    for tt=1:num
        % for every point
        x = pts_Image(tt, 1);
        y = pts_Image(tt, 2);
        
        local_dx = interp2(Dx_2D, x+1, y+1, 'linear');
        local_dy = interp2(Dy_2D, x+1, y+1, 'linear');
    
        dx = x + local_dx;
        dy = y + local_dy;
        
        if ( dx<1 )
            dx = 1;
        end
        
        if ( dx>header.xsize-1 )
            dx =header.xsize-1;
        end
        
        if ( dy<1 )
            dy = 1;
        end
        
        if ( dy>header.ysize-1 )
            dy =header.ysize-1;
        end
        
        pts_Deformed(tt, 1, kk) = dx;
        pts_Deformed(tt, 2, kk) = dy;
               
    end
end