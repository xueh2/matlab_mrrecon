function varargout = imuihs(varargin)
% IMUIHS Application M-file for imuihs.fig
%    FIG = IMUIHS launch imuihs GUI.
%    IMUIHS('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Jul-2002 14:27:19

if ~isstr(varargin{1})

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);
	movegui(fig, 'center')
	CX = varargin{1};
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	CPREV = im2double(thumb(CX));
	set(fig, 'CurrentAxes', handles.axePreview)
	handles.imgPreview = imshow(CPREV, 'notruesize');
	set(handles.axePreview, ...
		'Visible',		'off', ...
		'DrawMode',		'fast')
	set(handles.imgPreview, ...
		'EraseMode',	'xor', ...
		'UserData',		rgb2hsv(CPREV))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	guidata(fig, handles);
	set(fig, ...
		'HandleVisibility',		'callback', ...
		'Visible',				'on')

	waitfor(fig, 'Visible', 'off');
	if strcmp(get(handles.btnApply, 'UserData'), 'Apply')
		h = waitfig('Applying Hue/Saturation');
		valHue = get(handles.sldHue, 'Value');
		valSaturation = get(handles.sldSaturation, 'Value');
		valBrightness = get(handles.sldBrightness, 'Value');
		classcx = class(CX);
		CX = HSFCN(rgb2hsv(CX), valHue, valSaturation, valBrightness);
		switch classcx
		case 'uint8'
			varargout{1} = im2uint8(CX);
		case 'uint16'
			varargout{1} = im2uint16(CX);
		otherwise
			varargout{1} = CX;
		end
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
function varargout = sldHue_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldSaturation_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldBrightness_Callback(h, eventdata, handles, varargin)
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
set(handles.sldHue, 'Value', 0)
set(handles.sldSaturation, 'Value', 0)
set(handles.sldBrightness, 'Value', 0)

set(handles.txtHue, 'String', '0.00')
set(handles.txtSaturation, 'String', '0.00')
set(handles.txtBrightness, 'String', '0.00')

set(handles.chkPreview, 'Value', 1)
CX = get(handles.imgPreview, 'UserData');
set(handles.imgPreview, 'CData', hsv2rgb(CX))
% --------------------------------------------------------------------
function varargout = chkPreview_Callback(h, eventdata, handles, varargin)
if get(handles.chkPreview, 'Value')
	UpdatePreview(handles)
else
	CX = get(handles.imgPreview, 'UserData');
	set(handles.imgPreview, 'CData', hsv2rgb(CX))
end
% ====================================================================
function UpdatePreview(handles)
valHue = get(handles.sldHue, 'Value');
valSaturation = get(handles.sldSaturation, 'Value');
valBrightness = get(handles.sldBrightness, 'Value');
set(handles.txtHue, 'String', sprintf('%.2f', valHue) )
set(handles.txtSaturation, 'String', sprintf('%.2f', valSaturation) )
set(handles.txtBrightness, 'String', sprintf('%.2f', valBrightness) )

if get(handles.chkPreview, 'Value')
	CX = get(handles.imgPreview, 'UserData');
	CX = HSFCN(CX, valHue, valSaturation, valBrightness);
	set(handles.imgPreview, 'CData', CX)
end
% ====================================================================
function CY = HSFCN(CX, H, S, B)
LTBH = linspace(0, 1, 256) + H / 2;
for i = 1 : 256
	if LTBH(i) > 1
		LTBH(i) = LTBH(i) - 1;
	end
	if LTBH(i) < 0
		LTBH(i) = LTBH(i) + 1;
	end
end
LTBS = linspace(0, 1, 256) + sign(S) * abs(S .^ 2);
LTBS(LTBS > 1) = 1; LTBS(LTBS < 0) = 0;
LTBB = linspace(0, 1, 256) + sign(B) * abs(B .^ 2);
LTBB(LTBB > 1) = 1; LTBB(LTBB < 0) = 0;
CX = cat(3, ...
	grayxform( CX(:, :, 1), LTBH ), ...
	grayxform( CX(:, :, 2), LTBS ), ...
	grayxform( CX(:, :, 3), LTBB ) );
CY = hsv2rgb(CX);
