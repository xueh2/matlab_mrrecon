function imArrayplayer(im,delay_length)
% The function displays any animated GIF's in a figure window
%
% Demo : gifplayer; %plays the animated crystal.gif file
%
% Usage:
% ex: gifplayer('animated.gif',0.1); %name of the gif file and the delay in
% which to update the frames in the GIF file
%
% Vihang Patil, Oct 2006
% Copyright 2006-2007 Vihang Patil
% Email: vihang_patil@yahoo.com
% Created: 17th Oct 2006
%
% Revision:
% Date: 19th Oct 2006..Removed the setappdata and getappdata and used
% functions handling property. Multiple Gif's can be called upon which can
% be opened in new figure window.
% ex: figure;gifplayer;
% ex: figure;gifplayer('abcd.gif',0.1); and so on
% 
% P.N: PLease make sure to close the existing window in which the gif is
% currently being played and open a separate window for another GIF
% image.If another GIF is opened in the same window then the first timer
% continues to run even if you close the figure window.


if nargin<1
    delay_length = 0.2;%frame will be updated after 0.2 sec
elseif nargin <2
    delay_length = 0.2;%frame will be updated after 0.2 sec
end


%[handles.im,map] = imread(gif_image,'frames','all'); %read all frames of an gif image
s = size(im);
im2 = zeros(s(1), s(2),3,s(3));
im2(:,:,1,:) = reshape(im, [s(1), s(2), 1, s(3)]);
im2(:,:,2,:) = im2(:,:,1,:);
im2(:,:,3,:) = im2(:,:,1,:);
handles.im = im2;
handles.len = size(handles.im,4); % number of frames in the gif image
figure;
handles.h1 = imshow(uint8(handles.im(:,:,:,1)), 'Border', 'tight');%loads the first image along with its colormap
handles.count = 1;% intialise counter to update the next frame
set(handles.h1,'UserData',handles.count);%store the value of the counter in the image handles Userdata
%     setappdata(0,'h1',h1); %pass variable to Tmrfcn
%     setappdata(0,'im',im);%pass variable to Tmrfcn
%     setappdata(0,'count',count);%pass variable to Tmrfcn
%     setappdata(0,'len',len);%pass variable to Tmrfcn
handles.tmr = timer('TimerFcn', {@TmrFcn,handles},'BusyMode','Queue',...
    'ExecutionMode','FixedRate','Period',delay_length); %form a Timer Object
start(handles.tmr); %starts Timer

set(gcf,'CloseRequestFcn',{@CloseFigure,handles});


function TmrFcn(src,event,handles)
% function TmrFcn(varargin)
% handles=guihandles;
% h1 = handles.h1;
% im = handles.im;
% count = handles.count;
% len = handles.len;
% h1 = getappdata(0,'h1'); %gets the value of the variable from the parent fcn
% im = getappdata(0,'im'); %gets the value of the variable from the parent fcn
% count = getappdata(0,'count');%gets the value of the variable from the parent fcn
% len = getappdata(0,'len');%gets the value of the variable from the parent fcn

handles.count = get(handles.h1,'UserData');%gets the value of the variable from the parent fcn
set(handles.h1,'CData',uint8(handles.im(:,:,:,handles.count))); %update the frame in the axis
handles.count = handles.count + 1; %increment to next frame
set(handles.h1,'UserData',handles.count); %stores the count in Image userdata

if handles.count > handles.len %if the last frame is achieved intialise to first frame
    handles.count = 1;
    set(handles.h1,'UserData',handles.count);%stores the count in Image userdata
end
% setappdata(0,'count',count);%update the counter



function CloseFigure(src,event,handles)
% function CloseFigure(varargin)
stop(handles.tmr);delete(handles.tmr);%removes the timer from memory
closereq;
% d=timerfind;% stop(d);delete(d);

