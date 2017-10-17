function setimagescn(xlim,ylim,clim);
% function setimagescn(xlim,ylim,clim);
% function to set the image data for a figure created by imagescn.m
% 
% setimagescn may be used in conjunction with getimagedata to copy settings
%
% usage:
% 	setimagescn(xlim,ylim,clim); % sets axis properties for current figure window
%   setimagescn(xlim,ylim);

%      ***************************************
%      *  Peter Kellman  (kellman@nih.gov)   *
%      *  Laboratory for Cardiac Energetics  *
%      *  NIH NHLBI                          *
%      ***************************************

if nargin<2; return; end

h_axes=flipud(findobj(gcf,'type','axes'));

% set xlim
if ~isempty(xlim)
    if size(xlim,1)~=length(h_axes)
        for i=1:length(h_axes)
            set(h_axes(i),'xlim',xlim(1,:));
        end   
    else
        for i=1:length(h_axes)
            set(h_axes(i),'xlim',xlim(i,:));
        end   
    end
end

% set ylim
if ~isempty(ylim)
    if size(ylim,1)~=length(h_axes)
        for i=1:length(h_axes)
            set(h_axes(i),'ylim',ylim(1,:));
        end   
    else
        for i=1:length(h_axes)
            set(h_axes(i),'ylim',ylim(i,:));
        end   
    end
end
    
% set clim
if ~isempty(clim)
    if nargin==3
        if size(clim,1)~=length(h_axes)
            for i=1:length(h_axes)
                set(h_axes(i),'clim',clim(1,:));
            end   
        else
            for i=1:length(h_axes)
                set(h_axes(i),'clim',clim(i,:));
            end   
        end
    end
end

return