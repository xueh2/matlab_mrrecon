function varargout = imuistrel(varargin)
% IMUISTREL Application M-file for imuistrel.fig
%    FIG = IMUISTREL launch imuistrel GUI.
%    IMUISTREL('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 23-Jul-2002 13:32:55

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
	movegui(fig, 'center')

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
	
	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);
	handles.line0 = line([1 2], [1 2], ...
						'LineStyle', 		'none', ...
						'Marker',			's', ...
						'MarkerEdgeColor',	[0 0.5 1], ...
						'MarkerFaceColor',	'none', ...						
						'MarkerSize',		4);
	handles.line1 = line([1 2], [1 2], ...
						'LineStyle', 		'none', ...
						'Marker',			's', ...
						'MarkerEdgeColor',	[0 0.5 1], ...
						'MarkerFaceColor',	[0 0.5 1], ...							
						'MarkerSize',		4);
	guidata(fig, handles);
	UpdatePreview(handles)
	set(fig, ...
		'HandleVisibility', 	'callback', ...
		'Visible',				'on')

	waitfor(fig, 'Visible', 'off')

	if get(handles.btnOK, 'UserData')
		varargout{1} = get(handles.axePreview, 'UserData');
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

% ====================================================================
function varargout = rad1_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad2_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad3_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad4_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad5_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad6_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad7_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
function varargout = rad8_Callback(h, eventdata, handles, varargin)
RadionButtonCallback(handles)
% --------------------------------------------------------------------
function RadionButtonCallback(handles)
allhandles = allchild(handles.fig);
set(findobj(allhandles, 'Style', 'radiobutton'), 'Value', 0)
set(gcbo, 'Value', 1)
rtag = get(gcbo, 'Tag');
for h = transpose(allhandles)
	htag = get(h, 'Tag');
	if ~strcmp(htag, 'axePreview')
		if ~strcmp(get(h, 'Style'), 'radiobutton') & ...
				~strcmp(get(h, 'Style'), 'pushbutton')  ...
				& length(htag) > 4
			if strcmp(htag(4), rtag(4))
				set(h, 'visible', 'on')
			else
				set(h, 'visible', 'off')
			end
		end
	end
end
UpdatePreview(handles)
% ====================================================================
function varargout = sld1_1_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld2_1_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld3_1_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld3_2_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld5_1_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld5_2_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld6_2_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld6_3_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld7_1_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld7_2_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
function varargout = sld8_1_Callback(h, eventdata, handles, varargin)
SliderSizeCallback(handles)
% --------------------------------------------------------------------
function SliderSizeCallback(handles)
tag = get(gcbo, 'Tag');
hdltxt = findobj( allchild(handles.fig), 'Tag', ['txt', tag(end-2 : end)] );
val = get(gcbo, 'Value');
set(hdltxt, 'String', sprintf('%d', round(val) ))
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = sld4_1_Callback(h, eventdata, handles, varargin)
val = get(gcbo, 'Value');
set(handles.txt4_1, 'String', sprintf('%d', 3 * round(val) ))
UpdatePreview(handles)
function varargout = sld6_1_Callback(h, eventdata, handles, varargin)
val = get(gcbo, 'Value');
set(handles.txt6_1, 'String', sprintf('%d', 2 * round(val) + 1	 ))
UpdatePreview(handles)
% --------------------------------------------------------------------
function varargout = pop2_2_Callback(h, eventdata, handles, varargin)
UpdatePreview(handles)
% ====================================================================
function UpdatePreview(handles)
allhandles = allchild(handles.fig);
for h = transpose(findobj(allhandles, 'Style', 'radiobutton'))
	if get(h, 'Value')
		htag = get(h, 'Tag');
		streltype = str2num(htag(4));
		break
	end
end
switch streltype
case 1
	R = round(get(handles.sld1_1, 'Value'));
	STREL = strel('diamond', R);
case 2
	R = round(get(handles.sld2_1, 'Value'));
	switch get(handles.pop2_2, 'Value')
	case 1
		N = 0;
	case 2
		N = 4;
	case 3
		N = 6;
	case 4
		N = 8;
	end
	STREL = strel('disk', R, N);
case 3
	LEN = round(get(handles.sld3_1, 'Value'));
	DEG = get(handles.sld3_2, 'Value');
	STREL = strel('line', LEN, DEG);
case 4
	R = round(get(handles.sld4_1, 'Value')) * 3;
	STREL = strel('octagon', R);
case 5
	OFFSET = round([ get(handles.sld5_1, 'Value'), ...
					get(handles.sld5_2, 'Value')] );
	STREL = strel('pair', OFFSET);
case 6
	P = round(get(handles.sld6_1, 'Value')) * 2 + 1;
	V = round([ get(handles.sld6_2, 'Value'), ...
					get(handles.sld6_3, 'Value')] );
	STREL = strel('periodicline', P, V);
case 7
	R = round([ get(handles.sld7_1, 'Value'), ...
					get(handles.sld7_2, 'Value')] );
	STREL = strel('rectangle', R);
case 8
	R = round(get(handles.sld8_1, 'Value'));
	STREL = strel('square', R);
end
set(handles.axePreview, 'UserData', STREL)
nhood = getnhood(STREL);
nsize = size(nhood);
s = sum(sum(nhood));
xdt0 = zeros(1, s);
ydt0 = zeros(1, s);
xdt1 = zeros(1, s);
ydt1 = zeros(1, s);
k = 1;
for i = 1 : nsize(1)
	for j = 1 : nsize(2)
		if nhood(i, j)
			xdt1(k) = i;
			ydt1(k) = j;
		else
			xdt0(k) = i;
			ydt0(k) = j;	
		end	
		k = k + 1;
	end
end
set(handles.line0, ...
	'XData', 	xdt0, ...
	'YData', 	ydt0)
set(handles.line1, ...
	'XData', 	xdt1, ...
	'YData', 	ydt1)

% --------------------------------------------------------------------
function varargout = btnOK_Callback(h, eventdata, handles, varargin)
set(handles.btnOK, 'UserData', 1)
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = btnCancel_Callback(h, eventdata, handles, varargin)
set(handles.btnOK, 'UserData', 0)
set(handles.fig, 'Visible', 'off')
% --------------------------------------------------------------------
function varargout = fig_CloseRequestFcn(h, eventdata, handles, varargin)
set(handles.btnOK, 'UserData', 0)
set(handles.fig, 'Visible', 'off')
