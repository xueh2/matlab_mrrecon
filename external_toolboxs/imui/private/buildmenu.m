function hmnu = buildmenu(hfig)
%BUILDMENU - Private function of IMUI
%IMUI calls BUILDMENU while init GUI
%    Member of IMUI
%    Kotaimen.C, 2002/05 - 2002/07, All Rights Reserved


%%%%Root menus
hmnu.File	= lbm(hfig, '&File',	CbStr('File'),	{}, {});
hmnu.Edit	= lbm(hfig, '&Edit',	'',				{}, {});
hmnu.Adjust	= lbm(hfig, '&Adjust',	'',				{}, {});
%hmnu.Analyze= lbm(hfig, '&Analyze',	'',				{}, {});
hmnu.Filter = lbm(hfig, 'F&ilter',	'',				{}, {});
hmnu.View	= lbm(hfig, '&View',	'',				{}, {});
hmnu.Window	= lbm(hfig, '&Window',	CbStr('Window'),{}, {});
hmnu.Help	= lbm(hfig, '&Help', 	'',				{}, {});


%%%%File menu
hmnu.New	= lbm(hmnu.File, '&New...',		CbStr('New'),		{'off', 'N'}, {});
hmnu.Duplicate=lbm(hmnu.File,'&Duplicate',	CbStr('Duplicate'),	{'off', 'D'}, {});
hmnu.Open 	= lbm(hmnu.File, '&Open...',	CbStr('Open'),		{'on', 'O'},  {});
hmnu.Save 	= lbm(hmnu.File, '&Save as...',	CbStr('Save'),		{'off', 'S'}, {});
hmnu.CurDir = lbm(hmnu.File, 'C&urrent DIR',CbStr('CurDir'),	{}, {});
hmnu.Import = lbm(hmnu.File, '&Import from workspace',		CbStr('Import'),{'on'}, {});
hmnu.Export = lbm(hmnu.File, '&Export to workspace...', 	CbStr('Export'),{}, {});
hmnu.Rename = lbm(hmnu.File, '&Rename', 	CbStr('Rename'),	{'on'}, {});
hmnu.imshow = lbm(hmnu.File, 'i&mshow',		CbStr('imshow'),	{'on'}, {});
hmnu.imageview=lbm(hmnu.File,'image&view',	CbStr('imageview'),	{}, {});
hmnu.Prefences= lbm(hmnu.File, '&Toolbox prefences...',CbStr('Prefences'),		{'on', 'K'}, {});
hmnu.Close	= lbm(hmnu.File, '&Close', 		'close(gcf)' ,	{'on', 'W'}, {});

uimenu(hmnu.CurDir, 'Label', '?');
uimenu(hmnu.Import, 'Label', '?');


%%%%Edit menu
hmnu.Fade	= lbm(hmnu.Edit,	'&Fade...', 		CbStr('Fade'),		{'off', 'F'},{[], 'Fade'});
hmnu.History= lbm(hmnu.Edit,	'&History', 		'',					{},{});
hmnu.Resize	= lbm(hmnu.Edit,	'&Resize...',		CbStr('Resize'),	{'on'},{[], 'Resize'});
hmnu.Crop	= lbm(hmnu.Edit,	'&Crop', 			CbStr('Crop'), 		{},	{[], 'Crop'});
hmnu.Flip	= lbm(hmnu.Edit,	'F&lip && rotate', 	'',					{'off'},{});
hmnu.Split	= lbm(hmnu.Edit,	'&Split into','',						{'on'},{[0 0 0 0 1]});
hmnu.Merge	= lbm(hmnu.Edit,	'&Merge channel...',CbStr('Merge'), 	{'off'},{});
hmnu.Convert= lbm(hmnu.Edit,	'C&onvert to',		'', {'on'}, {});
%%%%Analyze menu
hmnu.fft	 = lbm(hmnu.Edit,	'FFT...',CbStr('FFT2'),	{'on'},			{[0 0 0 1 0], 'FFT'});
hmnu.dct	 = lbm(hmnu.Edit,	'DCT...',CbStr('DCT'),	{},				{[0 0 0 1 0], 'DCT'});
hmnu.impixel = lbm(hmnu.Edit,	'Color &picker', 	CbStr('impixel'),	{'on'},		{[0 0 0 1 1]});
hmnu.improfile=lbm(hmnu.Edit,	'Prof&ile',			CbStr('improfile'),	{},			{[0 0 0 1 1]});
hmnu.Contour  =lbm(hmnu.Edit,	'Conto&ur...', 		CbStr('Contour'),	{'on'},		{[0 0 0 1 1]});
hmnu.Histogram=lbm(hmnu.Edit,	'Histo&gram...', 	CbStr('Histogram'),	{'off', 'H'},{[0 0 0 1 1]});

%%Flip submenu
hmnu.FlipV	= lbm(hmnu.Flip,	'Flip &vertical',	CbStr('FlipVertical'),	{'off'},{[], 'Flip vertical'});
hmnu.FlipH	= lbm(hmnu.Flip,	'Flip &horizontal', CbStr('FlipHorizontal'),{'off'},{[], 'Flip horizontal'});
hmnu.RotateR= lbm(hmnu.Flip,	'Rotate &180', 		CbStr('Rotate_180'),	{'on'}, {[], 'Rotate 180'});
hmnu.RotateL= lbm(hmnu.Flip,	'Rotate 90 &right', CbStr('Rotate_90R'),	{'off'},{[], 'Rotate 90 right'});
hmnu.RotateL= lbm(hmnu.Flip,	'Rotate 90 &left',	CbStr('Rotate_90L'),	{'off'},{[], 'Rotate 90 left'});
hmnu.RotateN= lbm(hmnu.Flip,	'&Numberic...', 	CbStr('Rotate_N'),		{'on'},	{[], 'Rotate numberic'});
%%Split submenu
hmnu.SplitRGB = lbm(hmnu.Split, '&R-G-B channels',	CbStr('Split_RGB'),	{},{[0 0 0 0 1], 'Split into R-G-B channels'});
hmnu.SplitHSV = lbm(hmnu.Split, '&H-S-V channels',	CbStr('Split_HSV'),	{},{[0 0 0 0 1], 'Split into H-S-V channels'});
hmnu.SplitYCbCr = lbm(hmnu.Split, '&Y-Cb-Cr channels',	CbStr('Split_YCbCr'),	{},{[0 0 0 0 1], 'Split into Y-Cb-Cr channels'});
%%Convert to submenu
hmnu.ToBinary=lbm(hmnu.Convert,	'&Binary', 			CbStr('ToBinary'),	{}, {[0 0 0 1 1], 'Convert to binary'});
hmnu.ToGray	= lbm(hmnu.Convert,	'&Grayscale',		CbStr('ToGray'),	{}, {[0 1 0 0 1], 'Convert to grayscale'});
hmnu.ToRGB	= lbm(hmnu.Convert,	'&RGB Color', 		CbStr('ToRGB'),		{}, {[0 1 0 1 0], 'Convert to RGB color'});

hmnu.ToUINT8= lbm(hmnu.Convert,	'uint &8', 			CbStr('ToUINT8'),	{'on'}, {[], 'Convert to UINT8'});
hmnu.ToUINT16=lbm(hmnu.Convert,	'uint 1&6', 		CbStr('ToUINT16'),	{}, {[], 'Convert to UINT16'});
hmnu.ToDOUBLE=lbm(hmnu.Convert,	'&double', 			CbStr('ToDOUBLE'),	{}, {[], 'Convert to DOUBLE'});


%%%%Adjust menu
hmnu.Invert	= lbm(hmnu.Adjust, '&Invert', 			CbStr('Invert'), 	{'off', 'I'}, 	{[], 		  'Invert'});
hmnu.Equalize=lbm(hmnu.Adjust, '&Equalize', 		CbStr('Equalize'), 	{}, 			{[0 0 0 1 1], 'Histogram equalize'});
hmnu.AutoLevles=lbm(hmnu.Adjust, '&Auto levels', 	CbStr('AutoLevels'),{}, 			{[0 0 0 1 1], 'Auto levels'});
hmnu.Threshold=lbm(hmnu.Adjust, '&Threshold...', 		CbStr('Threshold'), {}, 			{[0 0 0 1 1], 'Threshold'});
hmnu.BC		= lbm(hmnu.Adjust, 'Brightness/&Contrast...',CbStr('BC'),	{'on'},				{[0 0 0 1 1], 'Brightness/Contrast'});
hmnu.Levels = lbm(hmnu.Adjust, '&Levels...', 		CbStr('Levels'), 	{'off',  'L'}, 	{[0 0 0 1 1], 'Levels'});
hmnu.Curves = lbm(hmnu.Adjust, 'Cur&ves...', 		CbStr('Curves'), 	{'off', 'M'}, 	{[0 0 0 1 1], 'Curves'});
hmnu.Hue 	= lbm(hmnu.Adjust, '&Hue/Saturation...',CbStr('HS'), 	{'on', 	'U'}, 	{[0 0 0 0 1], 'Hue/Saturation'});
hmnu.Desaturate=lbm(hmnu.Adjust, '&Desaturate',		CbStr('Desaturate'),{}, 			{[0 0 0 0 1], 'Desaturate'});
hmnu.ColorBalance=lbm(hmnu.Adjust, 'Color &balance...',CbStr('ColorBalance'),{'off','B'},{[0 0 0 0 1],'Color balance'});
hmnu.ChannelMixer=lbm(hmnu.Adjust, 'Channel mi&xer...',CbStr('ChannelMixer'),{}, 		{[0 0 0 0 1], 'Channel mixer'});
hmnu.Convolution=lbm(hmnu.Adjust, 'Convolutio&n...',	CbStr('Convolution'),{'on'}, 	{[0 0 0 1 1], 'Convolution'});
hmnu.Equalizer 	= lbm(hmnu.Adjust, 'Equalizer...',		CbStr('Equalizer'),	{'off', 'E'},{[0 0 0 1 1], 'Equalizer'});

%%%%Fliter menu
hmnu.class_blursharpen = lbm(hmnu.Filter, 'Blur && sharpen', '', {}, {});
hmnu.class_noise = lbm(hmnu.Filter, 'Noise', '', {}, {});
hmnu.class_edge = lbm(hmnu.Filter, 'Edge detection', '', {}, {});
hmnu.class_morphological = lbm(hmnu.Filter, 'Morphological', '', {}, {});
hmnu.class_others = lbm(hmnu.Filter, 'Others', '', {}, {});
D = dir(fullfile(MyLocation, 'private', 'filter_*.m'));
for i = 1 : length(D)
		P = feval( D(i).name(1 : end - 2) );
%		warning(sprintf(['Plug-in filter error in ', char(D(i).name),'.\n' lasterr]))
		switch lower(P.Class)
		case 'blur&sharpen'
			h_father = hmnu.class_blursharpen;
		case 'noise'
			h_father = hmnu.class_noise;
		case 'edge'
			h_father = hmnu.class_edge;
		case 'morphological'
			h_father = hmnu.class_morphological;
		otherwise
			h_father = hmnu.class_others;
		end
		h = lbm(h_father, 	P.FilterName, CbStr('FilterChild'), {}, ...
				{P.AvailableImageType, ['Filter : ',P.FilterName ]}, str2func(D(i).name(1 : end - 2)));
		eval(['hmnu.Filter_', num2str(i), ' = h;'])

end

%%%%View menu
hmnu.Hide	= lbm(hmnu.View,	'&Hide/Show menubar [SPACE]',CbStr('HideMenu'),	{},{});
hmnu.Zoom	= lbm(hmnu.View,	'&Zoom [Z]',		CbStr('Zoom'),		{'on'},{});
hmnu.Grid	= lbm(hmnu.View,	'&Grid [G]',		CbStr('Grid'),		{'off'},	{});
hmnu.ZoomIn	= lbm(hmnu.View	,	'Zoom &in [+]',		CbStr('ZoomIn'),	{'on'},	{});
hmnu.ZoomOut= lbm(hmnu.View,	'Zoom &out [-]',		CbStr('ZoomOut'),	{},		{});
hmnu.ZoomFit= lbm(hmnu.View,	'Fit &screen [*]',	CbStr('ZoomFit'),	{'on'},	{});
hmnu.ZoomFig= lbm(hmnu.View,	'Fit &figure [.]',	CbStr('ZoomFig'),	{},		{});
hmnu.ZoomTrue=lbm(hmnu.View,	'&Truesize [/]',	CbStr('ZoomTrue'),	{},		{});
hmnu.GridClr= lbm(hmnu.View,	'G&ridline color',	'',	{'on'},{});
hmnu.FigureClr = lbm(hmnu.View, '&Figure color',	'',	{},{});
hmnu.GridStyle = lbm(hmnu.View,	'Grid&line style',	'',	{'on'},{});
%%Figure Color sub-menu
hmnu.FigureClr0 = lbm(hmnu.FigureClr,  '&0 Custom...', 		CbStr('FigureColor'), {}, {});
hmnu.FigureClrD = lbm(hmnu.FigureClr,  '&1 System',			CbStr('FigureColor'), {'on'}, {},...
																 get(0,'defaultUicontrolBackgroundColor'));
hmnu.FigureClr1 = lbm(hmnu.FigureClr,  '&A White', 			CbStr('FigureColor'), {'on'}, {}, [ 1  1  1]);
hmnu.FigureClr2 = lbm(hmnu.FigureClr,  '&B Light gray', 	CbStr('FigureColor'), {}, {}, [.8 .8 .8]);
hmnu.FigureClr3 = lbm(hmnu.FigureClr,  '&C Medium gray', 	CbStr('FigureColor'), {}, {}, [.6 .6 .6]);
hmnu.FigureClr4 = lbm(hmnu.FigureClr,  '&D Dark gray', 		CbStr('FigureColor'), {}, {}, [.3 .3 .3]);
hmnu.FigureClr5 = lbm(hmnu.FigureClr,  '&E Black', 			CbStr('FigureColor'), {}, {}, [ 0  0  0]);
%%Gridline Color sub-menu
hmnu.GridClr0 = lbm(hmnu.GridClr,  '&0 Custom...', 		CbStr('GridColor'), {}, {});
hmnu.GridClr1 = lbm(hmnu.GridClr,  '&1 White', 			CbStr('GridColor'), {'on'}, {}, [ 1  1  1]);
hmnu.GridClr2 = lbm(hmnu.GridClr,  '&2 Gray', 			CbStr('GridColor'), {}, {}, 	[.5 .5 .5]);
hmnu.GridClr3 = lbm(hmnu.GridClr,  '&3 Black', 			CbStr('GridColor'), {}, {}, 	[ 0  0  0]);
hmnu.GridClr4 = lbm(hmnu.GridClr,  '&4 Red', 			CbStr('GridColor'), {}, {}, 	[ 1  0  0]);
hmnu.GridClr5 = lbm(hmnu.GridClr,  '&5 Green', 			CbStr('GridColor'), {}, {}, 	[ 0  1  0]);
hmnu.GridClr6 = lbm(hmnu.GridClr,  '&6 Blue', 			CbStr('GridColor'), {}, {}, 	[ 0  0  1]);
hmnu.GridClr7 = lbm(hmnu.GridClr,  '&7 Light yellow',	CbStr('GridColor'), {}, {}, 	[ 1  1 .8]);
hmnu.GridClr8 = lbm(hmnu.GridClr,  '&8 Light blue',		CbStr('GridColor'), {}, {}, 	[ 0 .7 .9]);
hmnu.GridClr9 = lbm(hmnu.GridClr,  '&9 Dark green', 	CbStr('GridColor'), {}, {}, 	[ 0 .5  0]);
hmnu.GridClr10= lbm(hmnu.GridClr,  '&A Navy blue',		CbStr('GridColor'), {}, {}, 	[ 0  0 .5]);
hmnu.GridClr11= lbm(hmnu.GridClr,  '&B Dark red',		CbStr('GridColor'), {}, {}, 	[.5  0  0]);
%%Gridline Style sub-menu
hmnu.GridStyle1 = lbm(hmnu.GridStyle,  '&1 Solid lines', 	CbStr('GridLineStyle'), {}, {}, '-' );
hmnu.GridStyle2 = lbm(hmnu.GridStyle,  '&2 Dashed lines', 	CbStr('GridLineStyle'), {}, {}, '--');
hmnu.GridStyle3 = lbm(hmnu.GridStyle,  '&3 Dotted lines', 	CbStr('GridLineStyle'), {}, {}, ':' );
hmnu.GridStyle4 = lbm(hmnu.GridStyle,  '&4 Dash-dot lines', CbStr('GridLineStyle'), {}, {}, '-.');

%%%%Window menu
uimenu(hmnu.Window, 'Label', '?');

%%%%Help menu
hmnu.Contents = lbm(hmnu.Help,		'&MATLAB Help', 'doc', {},{});
hmnu.Contents = lbm(hmnu.Help,		'&Image Processing Toolbox Help', 'doc images', {'on'},{});
hmnu.Contents = lbm(hmnu.Help,		'&About IMUI', CbStr('About'), {'on'},{});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Local Build Menu
%%%%varin1 = {MenuSeparator, MenuAccelerator}
function varargout = lbm(...
						ParentHandle, MenuLabel, CallbackStr, ....
						varin1, varin2, varargin)

switch length(varin1)
case 0
	MenuSeparator	= 'off';
	MenuAccelerator	= '';
case 1
	MenuSeparator	= varin1{1};
	MenuAccelerator	= '';
case 2
	MenuSeparator	= varin1{1};
	MenuAccelerator	= varin1{2};
end
switch length(varin2)
case 0
	AvailableImageType	= [0 1 1 1 1];
	HistoryActionName	= '';
case 1
	AvailableImageType	= varin2{1};
	HistoryActionName	= '';
case 2
	if isempty(varin2{1})
		AvailableImageType	= [0 1 1 1 1];
	else
		AvailableImageType	= varin2{1};
	end
	HistoryActionName	= varin2{2};
end
if ~isempty(varargin)
	CustomData = varargin{1};
else
	CustomData = '';
end
hmnu = uimenu( ParentHandle, ...
	'Label',				MenuLabel, ...
	'Accelerator',			MenuAccelerator, ...
	'Separator',			MenuSeparator, ...
	'Callback', 			CallbackStr, ...
	'UserData',	...
			struct( ...
				'AvailableImageType',	AvailableImageType, ...
				'ActionName',			HistoryActionName, ...
				'CustomData',			CustomData) ...
								);
if nargout > 0
	varargout{1} = hmnu;
end

%%%%
function stro = CbStr(stri)
stro = ['imui(  ''::::mnucb_', stri, ''' )'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MyLoaction - returns directory where imui.m locates
function R = MyLocation()
W = which('imui.m');
R = W(1 : strfind(W, 'imui.m') - 1);