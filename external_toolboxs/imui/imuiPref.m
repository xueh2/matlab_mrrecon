function varargout = imuiPref(varargin)
% IMUIPREF Application M-file for imuiPref.fig
%    FIG = IMUIPREF launch imuiPref GUI.
%    IMUIPREF('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 02-Jun-2002 10:35:30

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);

	%====================================================
	if strcmpi(iptgetpref('ImshowTruesize'), 'auto');
		set(handles.truesize, 'Value', 1)
	else
		set(handles.truesize, 'Value', 0)
	end
	if strcmpi(iptgetpref('ImshowAxesVisible'), 'off');
		set(handles.axes, 'Value', 0)
	else
		set(handles.axes, 'Value', 1)
	end
	if strcmpi(iptgetpref('TruesizeWarning'), 'on');
		set(handles.truesizewarning, 'Value', 1)
	else
		set(handles.truesizewarning, 'Value', 0)
	end
	if strcmpi(iptgetpref('ImshowBorder'), 'loose');
		set(handles.loose, 'Value', 1)
		set(handles.tight, 'Value', 0)
	else
		set(handles.loose, 'Value', 0)
		set(handles.tight, 'Value', 1)
	end
	%====================================================
	guidata(fig, handles);
	movegui(handles.imuiPref, 'center')
	set(handles.imuiPref, 'Visible', 'on')
	if nargout > 0
		varargout{1} = fig;
	end

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
function varargout = loose_Callback(h, eventdata, handles, varargin)
set(handles.loose, 'Value', 1)
set(handles.tight, 'Value', 0)
iptsetpref('ImshowBorder', 'loose');
% --------------------------------------------------------------------
function varargout = tight_Callback(h, eventdata, handles, varargin)
set(handles.loose, 'Value', 0)
set(handles.tight, 'Value', 1)
iptsetpref('ImshowBorder', 'tight');
% --------------------------------------------------------------------
function varargout = axes_Callback(h, eventdata, handles, varargin)
if get(h, 'Value')
	iptsetpref('ImshowAxesVisible', 'on');
else
	iptsetpref('ImshowAxesVisible', 'off')
end


% --------------------------------------------------------------------
function varargout = truesize_Callback(h, eventdata, handles, varargin)
if get(h, 'Value')
	iptsetpref('ImshowTruesize', 'auto');
else
	iptsetpref('ImshowTruesize', 'manual')
end


% --------------------------------------------------------------------
function varargout = truesizewarning_Callback(h, eventdata, handles, varargin)
if get(h, 'Value')
	iptsetpref('TruesizeWarning', 'on');
else
	iptsetpref('TruesizeWarning', 'off')
end

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
delete(handles.imuiPref)

% --------------------------------------------------------------------
function varargout = reset_Callback(h, eventdata, handles, varargin)
set(handles.loose, 'Value', 1)
set(handles.tight, 'Value', 0)
set(handles.axes, 'Value', 0)
set(handles.truesize, 'Value', 1)
set(handles.truesizewarning, 'Value', 1)

iptsetpref('ImshowTruesize', 'auto');
iptsetpref('ImshowAxesVisible', 'off');
iptsetpref('TruesizeWarning', 'on');
iptsetpref('ImshowBorder', 'loose');