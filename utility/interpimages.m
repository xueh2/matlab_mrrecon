function out = interpimages(method, imsize);
% function out = interpimages(method, imsize);
% function to get the image data from a figure created by imagescn.m
% interpolate, and re-display the zoomed images
%
% Interpolation method:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - bilinear interpolation (default)
%       'cubic'   - bicubic interpolation
%       'spline'  - spline interpolation
%
% Usage:
% 	out = interpimages(method,imsize);
% 	out = interpimages(method); % defaults to imsize=256
% 	out = interpimages;
% 	interpimages

%      ***************************************
%      *  Peter Kellman  (kellman@nih.gov)   *
%      *  Laboratory for Cardiac Energetics  *
%      *  NIH NHLBI                          *
%      ***************************************

if nargin==0; method='spline'; imsize=256; end % set default method and imsize
if nargin==1; imsize=256; end % set default imsize

% get image data from current figure window
[I, xlim, ylim, clim, axisposition]=getimagedata;
p = colormap;
n_axes=size(xlim,1);


% convert image data if not single or double datatype
datatype=class(I);
switch datatype
    case {'single', 'double'}
        % do nothing
    otherwise
        I=single(I);
end

% crop images to match image zoom using xlim and ylim
if ndims(I)==2
    Icrop{1}=I(ceil(ylim(1)):floor(ylim(2)),ceil(xlim(1)):floor(xlim(2)));
elseif ndims(I)==3
    if n_axes==1
        Icrop{1}=I(ceil(ylim(1)):floor(ylim(2)),ceil(xlim(1)):floor(xlim(2)),:);
    else
        for i=1:n_axes
            Icrop{i}=I(ceil(ylim(i,1)):floor(ylim(i,2)),ceil(xlim(i,1)):floor(xlim(i,2)),i);
		end
    end
else
    for i=1:n_axes
        Icrop{i}=I(ceil(ylim(i,1)):floor(ylim(i,2)),ceil(xlim(i,1)):floor(xlim(i,2)),:,i);
	end
end

% check if cropped image size >= imsize
for i=1:n_axes
    if size(Icrop{i},1)>=imsize; return; end
	if size(Icrop{i},2)>=imsize; return; end
end

% determine aspect ratio (same for all axes)
aspect_ratio=size(Icrop{1},1)/size(Icrop{1},2);

% interpolate images
if ndims(I)==2 | (ndims(I)==3 && n_axes > 1)
	for i=1:n_axes
        x=linspace(0,1,size(Icrop{i},2));
        y=linspace(0,1,size(Icrop{i},1))';
        if aspect_ratio >=1
            xi=linspace(0,1,floor(imsize/aspect_ratio));
            yi=linspace(0,1,imsize)';
        else
            xi=linspace(0,1,imsize);
            yi=linspace(0,1,floor(imsize*aspect_ratio))';
        end
%         Iinterp(:,:,i)=interp2(x,y,Icrop{i},xi,yi,method);
        
        Iinterp(:,:,i) = Matlab_gt_resize_2D_image(double(Icrop{i}), numel(yi), numel(xi), 5);
	end
else
	for i=1:n_axes
        x=linspace(0,1,size(Icrop{i},2));
        y=linspace(0,1,size(Icrop{i},1))';
        if aspect_ratio >=1
            xi=linspace(0,1,floor(imsize/aspect_ratio));
            yi=linspace(0,1,imsize)';
        else
            xi=linspace(0,1,imsize);
            yi=linspace(0,1,floor(imsize*aspect_ratio))';
        end
        tmp=Icrop{i};
        for frame=1:size(tmp,3)
%             Iinterp(:,:,frame,i)=interp2(x,y,tmp(:,:,frame),xi,yi,method);
            Iinterp(:,:,frame, i) = Matlab_gt_resize_2D_image(double(tmp(:,:,frame)), numel(yi), numel(xi), 5);
        end
	end
end

% get figure width
set(gcf,'Units','Inches')
figposition=get(gcf,'Position');
FigureWidth=figposition(3);

% determine row and column layout from axis positions
for i=1:n_axes
    left(i)=axisposition(i,1);
    bottom(i)=axisposition(i,2);
end
if n_axes==1
    cols=1;
else
	for cols=2:n_axes
        if left(cols)==left(1); cols=cols-1; break; end
	end
end
if n_axes==1 | ((n_axes>1) && (n_axes==cols))
    rows=1;
else
	for rows=cols+1:cols:n_axes
        if bottom(rows)==bottom(1); break; end
	end
    rows=ceil(rows/cols);
end

% display interpolated images
figure;
if ndims(Iinterp)==2
	imagescn(Iinterp,[],[rows cols],FigureWidth);
elseif ndims(Iinterp)==3
    if n_axes>1
        imagescn(Iinterp,[],[rows cols],FigureWidth);
    else
        imagescn(Iinterp,[],[rows cols],FigureWidth,3);
	end
elseif ndims(Iinterp)==4
    imagescn(Iinterp,[],[rows cols],FigureWidth,3);
end
colormap(p);

% set window level same as original images
figurehandle=gcf;
h_axes=flipud(findobj(figurehandle,'type','axes'));
for i=1:length(h_axes)
    set(h_axes(i),'clim',clim(i,:));
end

% interpolated output
if nargout>0
    out=Iinterp;
end

return