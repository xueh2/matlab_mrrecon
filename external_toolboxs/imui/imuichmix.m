function varargout = imuichmix(varargin)
% IMUICHMIX Application M-file for imuichmix.fig
%    FIG = IMUICHMIX launch imuichmix GUI.
%    IMUICHMIX('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Jul-2002 08:27:56

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
	ud = [1 1 1; 1 0 0; 0 1 0; 0 0 1];
	set(handles.popChannel, 'UserData', ud)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	guidata(fig, handles);
	set(fig, ...
		'HandleVisibility',		'callback', ...
		'Visible',				'on')

	waitfor(fig, 'Visible', 'off');

	if strcmp(get(handles.btnApply, 'UserData'), 'Apply')
		idx = get(handles.popChannel, 'Value');
		ud = get(handles.popChannel, 'UserData');
		h = waitfig('Mixing channels');
		if idx == 1
			CX = lincomb3(ud(1, 1) / 3, ud(1, 2) / 3, ud(1, 3) / 3, ...
				CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
%			CX = cat(3, CY, CY, CY);

		else
			CR = lincomb3(ud(2, 1), ud(2, 2), ud(2, 3), ...
				CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
			CG = lincomb3(ud(3, 1), ud(3, 2), ud(3, 3), ...
				CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
			CB = lincomb3(ud(4, 1), ud(4, 2), ud(4, 3), ...
				CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
			CX = cat(3, CR, CG, CB);
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

set(handles.sldRed, 'Value', ud(idx, 1))
set(handles.sldGreen, 'Value', ud(idx, 2))
set(handles.sldBlue, 'Value', ud(idx, 3))
UpdatePreview(handles)
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
function varargout = btnApply_Callback(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Apply')
set(handles.fig, 'Visible', 'off')

% --------------------------------------------------------------------
function varargout = btnCancel_Callback(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Cancel')
set(handles.fig, 'Visible', 'off')

% --------------------------------------------------------------------
function varargout = btnReset_Callback(h, eventdata, handles, varargin)
set(handles.sldRed, 'Value', 1)
set(handles.sldGreen, 'Value', 1)
set(handles.sldBlue, 'Value', 1)
set(handles.txtRed, 'String', '100%')
set(handles.txtGreen, 'String', '100%')
set(handles.txtBlue, 'String', '100%')
ud = [1 1 1; 1 0 0; 0 1 0; 0 0 1];
set(handles.popChannel, 'UserData', ud, 'Value', 2)
set(handles.chkPreview, 'Value', 1)
CX = get(handles.imgPreview, 'UserData');
set(handles.imgPreview, 'CData', CX)

% --------------------------------------------------------------------
function varargout = fig_CloseRequestFcn(h, eventdata, handles, varargin)
set(handles.btnApply, 'UserData', 'Cancel')
set(handles.fig, 'Visible', 'off')

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
valRed = get(handles.sldRed, 'Value');
valGreen =get(handles.sldGreen, 'Value');
valBlue = get(handles.sldBlue, 'Value');

set(handles.txtRed, 'String', sprintf('%d%%', round(valRed*100)) )
set(handles.txtGreen, 'String', sprintf('%d%%', round(valGreen*100)) )
set(handles.txtBlue, 'String', sprintf('%d%%', round(valBlue*100)) )

ud = get(handles.popChannel, 'UserData');
idx = get(handles.popChannel, 'Value');
ud(idx, :) = [valRed, valGreen, valBlue];
set(handles.popChannel, 'UserData', ud)

if get(handles.chkPreview, 'Value')
	CX = get(handles.imgPreview, 'UserData');
	if idx == 1
		CY = lincomb3(ud(1, 1) / 3, ud(1, 2) / 3, ud(1, 3) / 3, ...
			CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
		CX = cat(3, CY, CY, CY);
	else
		CR = lincomb3(ud(2, 1), ud(2, 2), ud(2, 3), ...
			CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
		CG = lincomb3(ud(3, 1), ud(3, 2), ud(3, 3), ...
			CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
		CB = lincomb3(ud(4, 1), ud(4, 2), ud(4, 3), ...
			CX(:, :, 1), CX(:, :, 2), CX(:, :, 3));
		CX = cat(3, CR, CG, CB);
	end
	set(handles.imgPreview, 'CData', CX)

end
% ====================================================================
function T = lincomb3(A, B, C, X, Y, Z)
T = imlincomb( C, Z, 1, imlincomb(A, X, B, Y) );



