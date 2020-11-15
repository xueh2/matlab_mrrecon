function varargout = imuiResize(varargin)
% IMUIRESIZE Application M-file for imuiResize.fig
%    FIG = IMUIRESIZE launch imuiResize GUI.
%    IMUIRESIZE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 24-Jul-2002 12:57:32

if ~ischar(varargin{1})  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	movegui(fig, 'center')
	CX = varargin{1};
	temp = size(CX);
	imgsize = [temp(2), temp(1)];

	set(handles.edtXPixel, 'String', num2str(imgsize(1)))
	set(handles.edtYPixel, 'String', num2str(imgsize(2)))
	set(handles.edtXPercent, 'String', num2str(100))
	set(handles.edtYPercent, 'String', num2str(100))

	set(handles.edtXPixel, 'UserData', imgsize(1))
	set(handles.edtYPixel, 'UserData', imgsize(2))
	set(handles.edtXPercent, 'UserData', 100)
	set(handles.edtYPercent, 'UserData', 100);

	set(handles.chkLock, 'Value', 1)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	set(fig, 'Visible', 'on')
	guidata(fig, handles);

	waitfor(fig, 'Visible')

	if get(handles.btnApply, 'UserData')
		h = waitfig('Resizing...');

		%lastXPixel = get(handles.edtXPixel, 'UserData');
		XPixel = str2num(get(handles.edtXPixel, 'String'));
		%lastYPixel = get(handles.edtYPixel, 'UserData');
		YPixel = str2num(get(handles.edtYPixel, 'String'));

		MethodString = {'nearest', 'bilinear', 'bicubic'};
		Method = MethodString{get(handles.popMethod, 'Value')};

		varargout{1} = imresize(CX,	[YPixel , XPixel], Method);
		close(h)
	else
		varargout{1} = -1;
	end
	delete(fig)


elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


% --------------------------------------------------------------------
function varargout = edtXPixel_Callback(h, eventdata, handles, varargin)

lastXPixel = get(handles.edtXPixel, 'UserData');
XPixel = str2num(get(handles.edtXPixel, 'String'));
set(handles.edtXPercent, 'String', num2str(round(XPixel/lastXPixel*100)) )
if get(handles.chkLock, 'Value')
	lastYPixel = get(handles.edtYPixel, 'UserData');
	YPixel = round(XPixel / lastXPixel * lastYPixel);
	set(handles.edtYPixel, 'String', num2str(YPixel))
	set(handles.edtYPercent, 'String', num2str(round(YPixel/lastYPixel*100)) )
end
% --------------------------------------------------------------------
function varargout = edtYPixel_Callback(h, eventdata, handles, varargin)

lastYPixel = get(handles.edtYPixel, 'UserData');
YPixel = str2num(get(handles.edtYPixel, 'String'));
set(handles.edtYPercent, 'String', num2str(round(YPixel/lastYPixel*100)) )
if get(handles.chkLock, 'Value')
	lastXPixel = get(handles.edtXPixel, 'UserData');
	XPixel = round(YPixel / lastYPixel * lastXPixel);
	set(handles.edtXPixel, 'String', num2str(XPixel))
	set(handles.edtXPercent, 'String', num2str(round(XPixel/lastXPixel*100)) )
end


% --------------------------------------------------------------------
function varargout = edtXPercent_Callback(h, eventdata, handles, varargin)
lastXPixel = get(handles.edtXPixel, 'UserData');
XPercent = str2num(get(handles.edtXPercent, 'String'));
set(handles.edtXPixel, 'String', num2str(round(lastXPixel * XPercent / 100)) )
if get(handles.chkLock, 'Value')
	set(handles.edtYPercent, 'String', num2str(XPercent))
	lastYPixel = get(handles.edtYPixel, 'UserData');
	set(handles.edtYPixel, 'String', num2str(round(lastYPixel * XPercent / 100)))
end


% --------------------------------------------------------------------
function varargout = edtYPercent_Callback(h, eventdata, handles, varargin)

lastYPixel = get(handles.edtYPixel, 'UserData');
YPercent = str2num(get(handles.edtYPercent, 'String'));
set(handles.edtYPixel, 'String', num2str(round(lastYPixel * YPercent / 100)) )
if get(handles.chkLock, 'Value')
	set(handles.edtXPercent, 'String', num2str(YPercent))
	lastXPixel = get(handles.edtXPixel, 'UserData');
	set(handles.edtXPixel, 'String', num2str(round(lastXPixel * YPercent / 100)))
end



% --------------------------------------------------------------------
function varargout = chkLock_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = popMethod_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = btnApply_Callback(h, eventdata, handles, varargin)

set(handles.figResize, 'Visible', 'off')
set(handles.btnApply, 'UserData', 1);

% --------------------------------------------------------------------
function varargout = btnCancel_Callback(h, eventdata, handles, varargin)

set(handles.figResize, 'Visible', 'off')
set(handles.btnApply, 'UserData', 0);

% --------------------------------------------------------------------
function varargout = btnReset_Callback(h, eventdata, handles, varargin)

lastXPixel = get(handles.edtXPixel, 'UserData');
set(handles.edtXPixel, 'String', num2str(lastXPixel) )
lastYPixel = get(handles.edtYPixel, 'UserData');
set(handles.edtYPixel, 'String', num2str(lastYPixel) )
set(handles.edtXPercent, 'UserData', 100)
set(handles.edtYPercent, 'UserData', 100)
set(handles.popMethod, 'Value', 3)
set(handles.chkLock, 'Value', 1)


% --------------------------------------------------------------------
function varargout = figResize_CloseRequestFcn(h, eventdata, handles, varargin)

set(handles.figResize, 'Visible', 'off')
set(handles.btnApply, 'UserData', 0);
