function [I, xlim, ylim, clim, position]=getimagedata(figurehandle);
% function [I, xlim, ylim, clim, position]=getimagedata(figurehandle);
% function to get the image data from a figure created by imagescn.m
% 
% usage:
% 	[I, xlim, ylim, clim, position]=getimagedata; % gets current figure window
% 	[I, xlim, ylim, clim, position]=getimagedata(figurehandle);

%      ***************************************
%      *  Peter Kellman  (kellman@nih.gov)   *
%      *  Laboratory for Cardiac Energetics  *
%      *  NIH NHLBI                          *
%      ***************************************

if nargin<1
    figurehandle=gcf;
end

h_axes=flipud(findobj(figurehandle,'type','axes'));

if isappdata(h_axes(1),'ImageData') % check if there is a "temporal" dimension
	for i=1:length(h_axes)
		I(:,:,:,i)=getappdata(h_axes(i),'ImageData');
%         I{i}=getappdata(h_axes(i),'ImageData');
		xlim(i,:)=get(h_axes(i),'xlim');
		ylim(i,:)=get(h_axes(i),'ylim');
		clim(i,:)=get(h_axes(i),'clim');
        position(i,:)=get(h_axes(i),'position');
	end   
else
	for i=1:length(h_axes)
        h_image=findobj(h_axes(i),'type','image');
	    I(:,:,i)=get(h_image,'cdata');
		xlim(i,:)=get(h_axes(i),'xlim');
		ylim(i,:)=get(h_axes(i),'ylim');
		clim(i,:)=get(h_axes(i),'clim');
        position(i,:)=get(h_axes(i),'position');
	end
end

return