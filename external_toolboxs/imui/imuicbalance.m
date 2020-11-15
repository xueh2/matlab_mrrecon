function varargout = imuicbalance(varargin)
% IMUICBALANCE Application M-file for imuicbalance.fig
%    FIG = IMUICBALANCE launch imuicbalance GUI.
%    IMUICBALANCE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Jul-2002 10:47:57

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
		h = waitfig('Adjusting color balance');
		valRed = get(handles.sldRed, 'Value');
		valGreen =get(handles.sldGreen, 'Value');
		valBlue = get(handles.sldBlue, 'Value');

		LTBR = linspace(0, 1, 256) + valRed * 100/255;
		LTBG = linspace(0, 1, 256) + valGreen * 100/255;
		LTBB = linspace(0, 1, 256) + valBlue * 100/255;

		LTBR(LTBR < 0) = 0;LTBR(LTBR > 1) = 1;
		LTBG(LTBG < 0) = 0;LTBG(LTBG > 1) = 1;
		LTBB(LTBB < 0) = 0;LTBB(LTBB > 1) = 1;
		CX = cat(3, ...
			grayxform(CX(:, :, 1), LTBR), ...
			grayxform(CX(:, :, 2), LTBG), ...
			grayxform(CX(:, :, 3), LTBB));

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
function varargout = sldRed_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldGreen_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sldBlue_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = chkPreserve_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = chkPreview_Callback(h, eventdata, handles, varargin)
if get(handles.chkPreview, 'Value')
	UpdatePreview(handles)
else
	CX = get(handles.imgPreview, 'UserData');
	set(handles.imgPreview, 'CData', CX)
end
% --------------------------------------------------------------------
function varargout = btnApply_Callback(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Apply')
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = fig_CloseRequestFcn(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Cancel')
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = btnCancel_Callback(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Cancel')
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = btnReset_Callback(h, eventdata, handles, varargin)
CX = get(handles.imgPreview, 'UserData');
set(handles.imgPreview, 'CData', CX)

set(handles.sldRed, 'Value', 0)
set(handles.sldGreen, 'Value', 0)
set(handles.sldBlue, 'Value', 0)
set(handles.txtRed, 'String', '0')
set(handles.txtGreen, 'String', '0')
set(handles.txtBlue, 'String', '0')
set(handles.chkPreview, 'Value', 1)
set(handles.chkPreserve, 'Value', 0)
% ====================================================================
function UpdatePreview(handles)
valRed = get(handles.sldRed, 'Value');
valGreen =get(handles.sldGreen, 'Value');
valBlue = get(handles.sldBlue, 'Value');
if get(handles.chkPreserve, 'Value')
	normalizefactor = -(valRed + valGreen + valBlue) / 3;
	valRed = valRed + normalizefactor;
	valGreen = valGreen + normalizefactor;
	valBlue	 = valBlue + normalizefactor;
	if valRed > 1; valRed = 1; end
	if valRed < -1; valRed = -1; end
	if valGreen > 1; valGreen = 1; end
	if valGreen < -1; valGreen = -1; end
	if valBlue > 1; valBlue = 1; end
	if valBlue < -1; valBlue = -1; end
end
set(handles.sldRed, 'Value', valRed)
set(handles.sldGreen, 'Value', valGreen)
set(handles.sldBlue, 'Value', valBlue)
set(handles.txtRed, 'String', sprintf('%d', round(valRed*100)) )
set(handles.txtGreen, 'String', sprintf('%d', round(valGreen*100)) )
set(handles.txtBlue, 'String', sprintf('%d', round(valBlue*100)) )

if get(handles.chkPreview, 'Value')
	CX = get(handles.imgPreview, 'UserData');

	LTBR = linspace(0, 1, 256) + valRed * 100/255;
	LTBG = linspace(0, 1, 256) + valGreen * 100/255;
	LTBB = linspace(0, 1, 256) + valBlue * 100/255;

	LTBR(LTBR < 0) = 0;LTBR(LTBR > 1) = 1;
	LTBG(LTBG < 0) = 0;LTBG(LTBG > 1) = 1;
	LTBB(LTBB < 0) = 0;LTBB(LTBB > 1) = 1;
	CX = cat(3, ...
		grayxform(CX(:, :, 1), LTBR), ...
		grayxform(CX(:, :, 2), LTBG), ...
		grayxform(CX(:, :, 3), LTBB));

	set(handles.imgPreview, 'CData', CX)
end




