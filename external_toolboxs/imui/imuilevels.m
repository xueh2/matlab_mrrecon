function varargout = imuilevels(varargin)
% IMUILEVELS Application M-file for imuilevels.fig
%    FIG = IMUILEVELS launch imuilevels GUI.
%    IMUILEVELS('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 19-Jul-2002 17:12:00

if ~isstr(varargin{1})

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	movegui(fig, 'center')
	CX = varargin{1};
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[C.L, X] = imhist(CX);
	if isrgb(CX)
		[C.R, X] = imhist( CX(:, :, 1) );
		[C.G, X] = imhist( CX(:, :, 2) );
		[C.B, X] = imhist( CX(:, :, 3) );
	else
		set(handles.popChannel, 'Enable', 'off')
	end
	set(fig, 'CurrentAxes', handles.axeHist)
	handles.histline = line([0 : 255], C.L, ...
		'Color',		[0.2 0.2 0.2]);
	set(handles.axeHist, 'UserData', C, 'XTick', round(linspace(0, 255, 7)))
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	CPREV = thumb(CX);
	set(fig, 'CurrentAxes', handles.axePreview)
	handles.imgPreview = imshow(CPREV, 'notruesize');
	set(handles.axePreview, ...
		'Visible',		'off', ...
		'DrawMode',		'fast')
	set(handles.imgPreview, ...
		'EraseMode',	'xor', ...
		'UserData',		CPREV)

	ud=[0 255 1.0 0 255;
		0 255 1.0 0 255;
		0 255 1.0 0 255;
		0 255 1.0 0 255];
	set(handles.popChannel, 'UserData', ud)
	guidata(fig, handles);
	set(fig, ...
		'HandleVisibility',		'callback', ...
		'Visible',				'on')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	waitfor(fig, 'Visible', 'off');

	if strcmp(get(handles.btnApply, 'UserData'), 'Apply')
		idx = get(handles.popChannel, 'Value');
		ud = get(handles.popChannel, 'UserData');
		h = waitfig('Applying Levels');
		if idx == 1
			CX = imadjust(CX, [ud(1, 1), ud(1, 2)] ./ 255, ...
								[ud(1, 4) ud(1, 5)] ./ 255, ud(1, 3));
		else
			CX = imadjust(CX, [ud(2, 1), ud(3, 1), ud(4, 1); ud(2, 2), ud(3, 2), ud(4, 2)]./255, ...
				[ud(2, 4), ud(3, 4), ud(4, 4); ud(2, 5), ud(3, 5), ud(4, 5)]./255, ...
				[ud(2, 3), ud(3, 3), ud(4, 3)]);
		end
		varargout{1} = CX;
		delete(h)
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
function varargout = popChannel_Callback(h, eventdata, handles, varargin)
idx = get(handles.popChannel, 'Value');
ud = get(handles.popChannel, 'UserData');
C = get(handles.axeHist, 'UserData');
switch idx
case 1
	set(handles.histline, ...
		'YData', 		C.L, ...
		'Color',		[0.2 0.2 0.2])
case 2
	set(handles.histline, ...
		'YData', 		C.R, ...
		'Color',		'r')
case 3
	set(handles.histline, ...
		'YData', 		C.G, ...
		'Color',		'g')
case 4
	set(handles.histline, ...
		'YData', 		C.B, ...
		'Color',		'b')
end
set(handles.sldBlackPoint, 'Value', ud(idx, 1));
set(handles.sldWhitePoint, 'Value', ud(idx, 2));
set(handles.sldGramma, 'Value', ud(idx, 3));
set(handles.sldLowOutput, 'Value', ud(idx, 4));
set(handles.sldHighOutput, 'Value', ud(idx, 5));
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldBlackPoint_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldWhitePoint_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldGramma_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldLowOutput_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldHighOutput_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = chkPreview_Callback(h, eventdata, handles, varargin)

if get(handles.chkPreview, 'Value')
	UpdatePreview(handles)
else
	CX = get(handles.imgPreview, 'UserData');
	set(handles.imgPreview, 'CData', CX)
end
% --------------------------------------------------------------------
function varargout = btnReset_Callback(h, eventdata, handles, varargin)
set(handles.popChannel, 'Value', 1)
set(handles.sldBlackPoint, 'Value', 0)
set(handles.sldWhitePoint, 'Value', 255)
set(handles.sldGramma, 'Value', 1)
set(handles.sldLowOutput, 'Value', 0)
set(handles.sldHighOutput, 'Value', 255)
set(handles.popChannel, 'UserData',  ...
		[0 255 1.0 0 255;
		0 255 1.0 0 255;
		0 255 1.0 0 255;
		0 255 1.0 0 255])
C = get(handles.axeHist, 'UserData');
set(handles.histline, ...
	'YData', 		C.L, ...
	'Color',		[0.2 0.2 0.2])

CX = get(handles.imgPreview, 'UserData');
set(handles.imgPreview, 'CData', CX)

UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = btnApply_Callback(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Apply')
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = btnCancel_Callback(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Cancel')
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = fig_CloseRequestFcn(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Cancel')
set(handles.fig, 'Visible', 'off')
% ====================================================================
function UpdatePreview(handles)
BlackPoint = get(handles.sldBlackPoint, 'Value');
WhitePoint = get(handles.sldWhitePoint, 'Value');
Gramma = get(handles.sldGramma, 'Value');
LowOutput = get(handles.sldLowOutput, 'Value');
HighOutput = get(handles.sldHighOutput, 'Value');

if BlackPoint > WhitePoint - 0.04
	BlackPoint = WhitePoint - 0.04;
	if BlackPoint < 0
		BlackPoint = 0;
	end
	set(handles.sldBlackPoint, 'Value', BlackPoint)
end

set(handles.txtBlackPoint, 'String', sprintf('%d', round(BlackPoint)))
set(handles.txtWhitePoint, 'String', sprintf('%d', round(WhitePoint)))
set(handles.txtGramma, 'String', sprintf('%.1f', Gramma))
set(handles.txtLowOutput, 'String', sprintf('%d', round(LowOutput)))
set(handles.txtHighOutput, 'String', sprintf('%d', round(HighOutput)))

ud = get(handles.popChannel, 'UserData');
idx = get(handles.popChannel, 'Value');
ud(idx, :) = [BlackPoint WhitePoint Gramma LowOutput HighOutput];
set(handles.popChannel, 'UserData', ud)

if get(handles.chkPreview, 'Value')
	CX = get(handles.imgPreview, 'UserData');
	if idx == 1
		CX = imadjust(CX, [BlackPoint/255 WhitePoint/255], ...
							[LowOutput/255 HighOutput/255], Gramma);
	else
		CX = imadjust(CX, [ud(2, 1), ud(3, 1), ud(4, 1); ud(2, 2), ud(3, 2), ud(4, 2)]./255, ...
			[ud(2, 4), ud(3, 4), ud(4, 4); ud(2, 5), ud(3, 5), ud(4, 5)]./255, ...
			[ud(2, 3), ud(3, 3), ud(4, 3)]);
	end
	set(handles.imgPreview, 'CData', CX)
end



