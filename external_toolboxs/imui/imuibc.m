function varargout = imuibc(varargin)
% IMUIBC Application M-file for imuibc.fig
%    FIG = IMUIBC launch imuibc GUI.
%    IMUIBC('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Jul-2002 13:57:35

if ~isstr(varargin{1})

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);

	movegui(fig, 'center')
	CX = varargin{1};
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

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	guidata(fig, handles);
	set(fig, ...
		'HandleVisibility',		'callback', ...
		'Visible',				'on')

	waitfor(fig, 'Visible', 'off');
	if strcmp(get(handles.btnApply, 'UserData'), 'Apply')
		h = waitfig('Applying Brightness/Contrast');
		valBrightness = get(handles.sldBrightness, 'Value');
		valContrast = get(handles.sldContrast, 'Value');
		CX = BCFCN(CX, valBrightness, valContrast);
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
function varargout = sldBrightness_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldContrast_Callback(h, eventdata, handles, varargin)
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
% --------------------------------------------------------------------
function varargout = btnReset_Callback(h, eventdata, handles, varargin)
set(handles.sldBrightness, 'Value', 0)
set(handles.sldContrast, 'Value', 0)
set(handles.txtBrightness, 'String', '0.00')
set(handles.txtContrast, 'String', '0.00')
set(handles.chkPreview, 'Value', 1)
CX = get(handles.imgPreview, 'UserData');
set(handles.imgPreview, 'CData', CX)

% --------------------------------------------------------------------
function varargout = chkPreview_Callback(h, eventdata, handles, varargin)
if get(handles.chkPreview, 'Value')
	UpdatePreview(handles)
else
	CX = get(handles.imgPreview, 'UserData');
	set(handles.imgPreview, 'CData', CX)
end
% ====================================================================
function UpdatePreview(handles)
valBrightness = get(handles.sldBrightness, 'Value');
valContrast = get(handles.sldContrast, 'Value');

set(handles.txtBrightness, 'String', sprintf('%.2f', valBrightness) )
set(handles.txtContrast, 'String', sprintf('%.2f', valContrast) )

if get(handles.chkPreview, 'Value')
	CX = get(handles.imgPreview, 'UserData');
	CX = BCFCN(CX, valBrightness, valContrast);
	set(handles.imgPreview, 'CData', CX)
end
% ====================================================================
function CY = BCFCN(CX, B, C)
switch C
case 0
	C = 1;
case 1
	C = 1000;
otherwise
	if C > 0
		C = 1 / (1 - C);
	else
		C = 1 + C;
	end
end

LTB = C .* ( linspace(0, 1, 256) - 0.5) + B + 0.5;

%Y= C( X - 0.5 )+B+0.5

LTB(LTB > 1) = 1;
LTB(LTB < 0) = 0;
CY = grayxform(CX, LTB);




