function imui(varargin)
% IMUI - Image Processing Toolbox Utilities for MATLAB 6.1
% Under development
% Kotaimen.C, 2002/5 - 200?/?, All rights reserved.
% Connect to kotaimen_c@citiz.net.
%   imui
%   imui('ImageFileName')
%   imui ImageFileName
%   imui(ImageData)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-Initialize
if isempty(getImuiFigureHandles)
	V = str2num((version('-release')));
	if V < 12.1
		error('IMUI requires MATLAB 6.1 (release 12.1) or higher.')
	end
	try
		iptgetpref;
	catch
		error('IMUI requires Image Processing Toolbox.')
	end
	D = get(0, 'ScreenDepth');
	if D < 16
		error('IMUI requires a 16-bit Monitor.')
	end
end
if nargin == 0
	PreAction = 'NewImage';
else
	if isstr(varargin{1})
		if length(varargin{1}) >= 4
			if varargin{1}(1 : 4) == '::::'
				PreAction = 'Callback';
			else
				PreAction = 'ImageFile';
			end
		else
			PreAction = 'ImageFile';
		end
	else
		PreAction = 'ImageData';
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
switch PreAction
case 'NewImage'
	CX = im2uint8(ones(300, 400));
	imgTitle = ['image', num2str(length(getImuiFigureHandles) + 1)];
	Action = '::::BuildGUI';
case 'ImageFile'
	try
		imgFileInfo = imfinfo(varargin{1});
	catch
		errstr = sprintf(['Unable to open "',varargin{1}, ...
			'", file may not exist or unknown file type.\n', ...
			'Type "help imread" for more information.']);
		error(errstr)
		return
	end
 	h = waitfig(sprintf('Opening %s', varargin{1}));
	if strcmpi(imgFileInfo.ColorType, 'indexed')
		[CX, CMAP] = imread(varargin{1});
		CX = im2uint8( ind2rgb(CX, CMAP) );
	else
		CX = imread(varargin{1});
	end
	imgType = imagetype(CX);
	imgTitle = varargin{1};
	Action = '::::BuildGUI';
 	delete(h)
case 'ImageData'
	CX = varargin{1};
	imgType = imagetype(CX);
	if nargin > 1
		imgTitle = varargin{2};
	else
		imgTitle = inputname(1);
		if isempty(imgTitle)
			imgTitle = ['image', num2str(length(getImuiFigureHandles) + 1)];;
		end
	end
	Action = '::::BuildGUI';
case 'Callback'
	Action = varargin{1};
	ud = get(gcbf, 'UserData');
otherwise
	error('Invalid PreAction')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback
switch Action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::BuildGUI'
	hdlFig = figure( ...
		'Name',			'imui', ...
		'Units',		'pixel', ...
		'Menubar',		'none', ...
		'NumberTitle',	'off', ...
		'Resize',		'on', ...
		'Visible',		'off', ...
		'Interruptible','off', ...
		'Busyaction', 	'queue', ...
		'Tag',			'::::imuiFigure');
	hdlTxt = uicontrol(hdlFig, ...
		'Style',		'text', ...
		'Units',		'normalized', ...
		'Position',		[0.1 0.1/3 0.8, 0.1/3], ...
		'FontName',		'Courier New', ...
		'FontSize',		9, ...
		'FontWeight',	'bold', ...
		'String',		'');
	hdlAxe = axes( ...
		'Units',		'normalized', ...
		'Position',		[0 0 1 1]);

	hdlImg = imshow(im2uint8(CX), 'notruesize');

	hdlMenu = buildmenu( hdlFig );

	set(hdlAxe, ...
		'Visible', 		'off', ...
		'DrawMode',		'fast', ...
		'Box',			'on', ...
		'Layer',		'top', ...
		'XAxisLocation','top', ...
		'GridLineStyle',':', ...
		'TickLength',	[0.001 0.025], ...
		'XGrid',		'on', ...
		'YGrid',		'on', ...
		'FontSize',		8)
	set(hdlImg, ...
		'UserData',		CX, ...
		'CDataMapping',	'scaled')

	ifZoom = 0;
	ifGrid = 0;
	clrGrid = [0 0 0.5];
	umtoggle(hdlMenu.GridClr10);
	clrFigure = [.8 .8 .8];
	umtoggle(hdlMenu.FigureClr2);
	styGrid = ':';
	umtoggle(hdlMenu.GridStyle3);
	ud = struct( ...
		'FigureHandle',		hdlFig, ...
		'TextHandle',		hdlTxt, ...
		'AxesHandle',		hdlAxe, ...
		'ImageHandle',		hdlImg, ...
		'MenuHandles',		hdlMenu, ...
		'ImageType',		imagetype(CX), ...
		'ImageClass',		class(CX), ...
		'ImageTitle',		imgTitle, ...
		'ImageSize',		[size(CX, 2), size(CX, 1)], ...
		'ZoomState',		ifZoom, ...
		'GridState',		ifGrid, ...
		'GridColor',		clrGrid, ...
		'FigureColor',		clrFigure, ...
		'GridLineStyle',	styGrid, ...
		'DefaultOversampleMethod',	'bicubic', ...
		'HistoryData',		struct([]), ...
		'HistoryIndex',		0, ...
		'LastHistoryIndex', 0, ...
		'HistoryCount',		0 ...
	);
	truesize( hdlFig )
	movegui(hdlFig, [24, -24] .* (length(getImuiFigureHandles) + 1) )

	set(hdlFig, ...
		'UserData',			ud, ...
		'Visible',			'on', ...
		'HandleVisibility',	'callback', ...
		'Colormap',			gray(256), ...
		'KeyPressFcn',		' imui(  ''::::KeyPress''  ) ', ...
		'CloseRequestFcn',	' imui(  ''::::CloseRequest''  ) ', ...
		'ResizeFcn',		' imui(  ''::::Resize''  ) ', ...
		'WindowButtonMotionFcn', 	'imui( ''::::ButtonMotion'') ', ...
		'WindowButtonDownFcn', 		'imui( ''::::ButtonDown'') ', ...
		'WindowButtonUpFcn', 		'imui( ''::::ButtonUp'') ')

	AppendToHistory(ud, 'New/Create/Open')
	UpdateGUIState(ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% UpdateImage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::UpdateImage'
	NewCX = varargin{2};
	mnuud = get(gcbo, 'UserData');
	NewImgType = imagetype(NewCX);
	NewImgClass = class(NewCX);
	NewImgSize = [size(NewCX, 2), size(NewCX, 1)];
	LastActionName = mnuud.ActionName;

	if all(NewImgSize == ud.ImageSize)% & strcmp(NewImgType, ud.ImageType)
		if strcmp(NewImgClass, 'double')
			NewCX(NewCX > 1) = 1;
			NewCX(NewCX < 0) = 0;
			NewImgType = imagetype(NewCX);
		end

		set(ud.ImageHandle, 'UserData', NewCX)
		DisplayFunction(ud.ImageHandle)

		ud.ImageClass = NewImgClass;
		ud.ImageType = NewImgType;
		ud.ImageSize = NewImgSize;
		AppendToHistory(ud, mnuud.ActionName)

		UpdateGUIState(ud)
	else
		imui(NewCX, [ud.ImageTitle, ' - ' ,mnuud.ActionName])
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RESIZEFCN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::Resize'
	UpdateImageDisplayRatio(ud)
case '::::ButtonDown'
	UpdateImageDisplayRatio(ud)
case '::::ButtonUp'
	UpdateImageDisplayRatio(ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% KEYPRESSFCN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::KeyPress'
	keychar = get(ud.FigureHandle, 'CurrentCharacter');
	switch keychar
	case ' '
		imui('::::mnucb_HideMenu')
	case {'Z', 'z'}
		imui('::::mnucb_Zoom')
	case {'G', 'g', ';'}
		imui('::::mnucb_Grid')
	case {'+', '='}
		imui('::::mnucb_ZoomIn');
	case {'-', '_'}
		imui('::::mnucb_ZoomOut');
	case '/'
		imui('::::mnucb_ZoomTrue');
	case '*'
		imui('::::mnucb_ZoomFit');
	case '.'
		imui('::::mnucb_ZoomFig');
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CLOSEREQUESTFCN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::CloseRequest'
	ud = get(gcf, 'UserData');
	for i = 1 : ud.HistoryCount
		delete([ud.HistoryData(i).FileName, '.mat'])
	end
	delete(gcf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WINDOWBUTTONMOTIONFCN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::ButtonMotion'

	if ud.GridState
		temp = get(ud.AxesHandle, 'CurrentPoint');
		cpX = round( temp(1, 1) );
		cpY = round( temp(1, 2) );

		if all([cpX <= ud.ImageSize(1), cpX >=1 ...
				cpY <= ud.ImageSize(2), cpY >=1])
			CX = get(ud.ImageHandle, 'UserData');
			switch ud.ImageClass
			case 'uint8'
				CX = double(CX(cpY, cpX, :));
			case 'uint16'
				CX = round( double(CX(cpY, cpX, :)) / 255 );
			otherwise
				CX = round(CX(cpY, cpX, :) * 255);
			end

			if strcmp(ud.ImageType, 'RGB')
				textstr = sprintf( 'X : %4d, Y : %4d;  RGB(%3d, %3d, %3d)', ...
					cpX, cpY, CX);
			else
				textstr = sprintf( 'X : %4d, Y : %4d;  Black(%3d) ', ...
					cpX, cpY, CX);
			end
		else
			textstr = 'Out of image';
		end
		set(ud.TextHandle, 'String', textstr)
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% File menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_File'
	if length(pwd) < 60
		set(ud.MenuHandles.CurDir, 'Label', pwd)
	else
		set(ud.MenuHandles.CurDir, 'Label', 'Current DIR')
	end
case '::::mnucb_New'
	prompt = {sprintf('Enter image name :\n'), ...
			sprintf('Enter X size (in pixels) :\n'), ...
			sprintf('Enter Y size (in pixels) :\n')};
	default = {['image', num2str(length(getImuiFigureHandles) + 1)], ...
			num2str(ud.ImageSize(1)), num2str(ud.ImageSize(2))};
	dlgtitle = 'New image';
	answer = inputdlg(prompt, dlgtitle, 1, default);
	if ~isempty(answer)
		imgTitle = answer{1};
		Xsize = round(str2num(answer{2}));
		Ysize = round(str2num(answer{3}));
		if Xsize >= 0 & Ysize >=0
			imui(im2uint8(ones(Ysize, Xsize)), imgTitle)
		end
	end
case '::::mnucb_Duplicate'
	imui(get(ud.ImageHandle, 'UserData'), [ud.ImageTitle, ' - copy'])
case '::::mnucb_Open'
	[filename, pathname] = uigetfile( {...
		'*.jpg;*.tif;*.gif;*.bmp;*.png', 'All image files(*.jpg,*.tif,*.gif,*.bmp,*.png)';
		'*.jpg;*.jpeg', 'JPEG files(*.jpg)';
		'*.gif', 'GIF files(*.gif)';
		'*.tif;*.tiff', 'TIFF files(*.tif)';
		'*.bmp', 'BMP files(*.bmp)';
		'*.png', 'PNG files(*.png)';
		'*.*', 'All Files (*.*)'}, 'Open an image');
	if ischar(filename)
		imui([pathname, filename])
	end
case '::::mnucb_Save'
	[filename, pathname] = uiputfile( ...
		ud.ImageTitle, 'Save as...');

	if ischar(filename)
		try
			imwrite(get(ud.ImageHandle, 'UserData'), ...
				[pathname, filename])
		catch
			errordlg(sprintf('Error while saving file: %s.\n%s', ...
			[pathname, filename], lasterr), 'IMUI', 'modal');
		end
	end
case '::::mnucb_CurDir'
	fileext = {'jpg', 'tif', 'gif', 'bmp', 'png'};
	filelist = '';
	%Build file list
	for i = 1 : length(fileext)
		FS = dir(['*.', char(fileext{i})]);
		if ~isempty(FS)
			for j = 1 : length(FS)
				filelist = strvcat(filelist, FS(j).name);
			end
		end
	end
	filelist = sortrows(filelist);

	hsubmenu = allchild(ud.MenuHandles.CurDir);
	if ~isempty(filelist)
		for i = 1 : size(filelist, 1);
			mnulabel = deblank( filelist(i, :) );
			uimenu( ud.MenuHandles.CurDir, ...
				'Label',		mnulabel, ...
				'Tag',			mnulabel, ...
				'Callback',		'imui(  ''::::mnucb_CurDir_ChildCbFCN''  )');
		end
	else
		uimenu( ud.MenuHandles.CurDir, ...
			'Label',		'There is no image file in current directory', ...
			'Enable',		'off');
	end
	if ~isempty(hsubmenu)
		delete(hsubmenu)
	end
case '::::mnucb_CurDir_ChildCbFCN'
	imui(get(gcbo, 'Tag'))
case '::::mnucb_Import'
	workspacevars = evalin('base', 'whos');
	varlist = '';
	for i = 1 : length(workspacevars)
		varclass = workspacevars(i).class;
		if any( strcmpi(varclass, {'double'; 'uint8'; 'uint16'}) )
			varlist = strvcat(varlist, workspacevars(i).name);
		end
	end
	hsubmenu = get(ud.MenuHandles.Import, 'Children');
	if ~isempty(varlist)
		for i = 1 : size(varlist, 1);

			varinf = evalin('base', ['whos(''',varlist(i, :), ''')']);

			varname = char(varinf.name);
			if ~isempty(varname)
				uimenu( ud.MenuHandles.Import, ...
					'Label',		[varname, ' -- ',...
								mat2str(varinf.size), ', ', char(varinf.class)], ...
					'Tag',			varname, ...
					'Callback',		'imui(  ''::::mnucb_Import_ChildCbFCN''  )');
			end
		end
	else
		uimenu( ud.MenuHandles.Import, ...
			'Label',		'Workspace contains no numberic variable', ...
			'Enable',		'off');
	end
	if ~isempty(hsubmenu)
		delete(hsubmenu)
	end
case '::::mnucb_Import_ChildCbFCN'
	imvar = get(gcbo, 'Tag');
	evalin('base', [  'imui(', imvar, ');'	 ] );
case '::::mnucb_Export'
	prompt = {sprintf('Enter varaible name :\n')};
	default = {strtok(ud.ImageTitle, '. ')};
	dlgtitle = 'Export to workspace';
	answer = inputdlg(prompt, dlgtitle, 1, default);
	if ~isempty(answer)
		varName = answer{1};
		if isvarname(varName)
			assignin('base', varName, get(ud.ImageHandle, 'UserData'));
		else
			prompt = sprintf(['\\bf', varName , '\\rm is not a valid varaible name. \n\n', ...
			'A valid variable name is a character string of letters, digits and', ...
    		'underscores, with length <= 31 and the first character a letter.']);
			h = errordlg(prompt, 'IMUI', ...
				struct('Interpreter', 'tex', 'WindowStyle', 'modal') );
		end
	end
case '::::mnucb_imshow'
	figure
	imshow(get(ud.ImageHandle, 'UserData'));
case '::::mnucb_imageview'
	imageview(get(ud.ImageHandle, 'UserData'));
case '::::mnucb_Prefences'
	imuiPref
case '::::mnucb_Rename'
	prompt = sprintf( 'Enter new name :\n');
	newname = inputdlg(prompt, 'Rename image', 1, {ud.ImageTitle});
	if ~isempty(newname)
		ud.ImageTitle = char(newname);
		UpdateGUIState(ud)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% EDIT menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_Fade'
	CX = imuifade(ud);
	if ~isempty(CX)
		imui('::::UpdateImage', CX)
	end
case '::::mnucb_History_ChildCbFCN'
	set(allchild(ud.MenuHandles.History), 'Checked', 'off')
	set(gcbo, 'Checked', 'on')
	ud.HistoryIndex = get(gcbo, 'UserData');
	mnuHistoryData = ud.HistoryData(ud.HistoryIndex);

	LOADDATA = load(mnuHistoryData.FileName);
	set(ud.ImageHandle, 'UserData', LOADDATA.CX)
	DisplayFunction(ud.ImageHandle)
	ud.ImageType = ImageType(LOADDATA.CX);
	ud.ImageClass = class(LOADDATA.CX);
	UpdateGUIState(ud)
case '::::mnucb_Resize'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuiresize(CX);
	if CX ~= -1
		imui('::::UpdateImage', CX)
	end
case '::::mnucb_Crop'
	DisableAllMenus(ud)
	set(ud.FigureHandle, 'Name', ...
			'IMUI - Crop : click and drag.')
	rect = getrect(ud.FigureHandle);
	CX = get(ud.ImageHandle, 'UserData');
	imui('::::UpdateImage', imcrop(CX, rect));
	EnableAllMenus(ud)

case '::::mnucb_FlipVertical'
	CX = get(ud.ImageHandle, 'UserData');
	if strcmp(ud.ImageType, 'RGB')
		imui('::::UpdateImage', ...
			cat(3, ...
			flipud(CX(:,:,1)), flipud(CX(:,:,2)), flipud(CX(:,:,3)) ) )
	else
		imui('::::UpdateImage', flipud(CX))
	end
case '::::mnucb_FlipHorizontal'
	CX = get(ud.ImageHandle, 'UserData');
	if strcmp(ud.ImageType, 'RGB')
		imui('::::UpdateImage', ...
			cat(3, ...
			fliplr(CX(:,:,1)), fliplr(CX(:,:,2)), fliplr(CX(:,:,3)) ) )
	else
		imui('::::UpdateImage', fliplr(CX))
	end
case '::::mnucb_Rotate_180'
	CX = get(ud.ImageHandle, 'UserData');
	if strcmp(ud.ImageType, 'RGB')
		imui('::::UpdateImage', ...
			cat(3, ...
			rot90(CX(:,:,1), 2), rot90(CX(:,:,2), 2), rot90(CX(:,:,3), 2) ) )
	else
		imui('::::UpdateImage', rot90(CX, 2))
	end
case '::::mnucb_Rotate_90L'
	CX = get(ud.ImageHandle, 'UserData');
	if strcmp(ud.ImageType, 'RGB')
		imui('::::UpdateImage', ...
			cat(3, ...
			rot90(CX(:,:,1)), rot90(CX(:,:,2)), rot90(CX(:,:,3)) ) )
	else
		imui('::::UpdateImage', rot90(CX))
	end
case '::::mnucb_Rotate_90R'
	CX = get(ud.ImageHandle, 'UserData');
	if strcmp(ud.ImageType, 'RGB')
		imui('::::UpdateImage', ...
			cat(3, ...
			rot90(CX(:,:,1), 3), rot90(CX(:,:,2), 3), rot90(CX(:,:,3), 3) ) )
	else
		imui('::::UpdateImage', rot90(CX, 3))
	end
case '::::mnucb_Rotate_N'
	prompt = {sprintf('Enter rotate angle(clockwise, in degree) :\n')};
	default = {'0'};
	dlgtitle = 'Rotate image';
	answer = inputdlg(prompt, dlgtitle, 1, default);
	if ~isempty(answer)
		if all(size(str2num(answer{1}))== [1 1])
			h = waitfig('Rotating');
			CX = get(ud.ImageHandle, 'UserData');
			imui('::::UpdateImage', ...
				imrotate(CX, -str2num(answer{1}), ud.DefaultOversampleMethod) )
			delete(h)
		else
			prompt = sprintf(['\\bf', answer{1} , '\\rm is not a angle. \n']);
			h = errordlg(prompt, 'IMUI', ...
				struct('Interpreter', 'tex', 'WindowStyle', 'modal') );
		end
	end
case '::::mnucb_Split_RGB'
	CX = get(ud.ImageHandle, 'UserData');
	imui(CX(:, :, 1), [ud.ImageTitle, ': R Channel'])
	imui(CX(:, :, 2), [ud.ImageTitle, ': G Channel'])
	imui(CX(:, :, 3), [ud.ImageTitle, ': B Channel'])
case '::::mnucb_Split_HSV'
	h = waitfig('Converting to H-S-V image');
	CX = rgb2hsv(get(ud.ImageHandle, 'UserData'));
	delete(h);
	imui(CX(:, :, 1), [ud.ImageTitle, ': H Channel'])
	imui(CX(:, :, 2), [ud.ImageTitle, ': S Channel'])
	imui(CX(:, :, 3), [ud.ImageTitle, ': V Channel'])
case '::::mnucb_Split_YCbCr'
	h = waitfig('Converting to Y-Cb-Cr image');
	CX = rgb2ycbcr(get(ud.ImageHandle, 'UserData'));
	delete(h);
	imui(CX(:, :, 1), [ud.ImageTitle, ': Y Channel'])
	imui(CX(:, :, 2), [ud.ImageTitle, ': Cb Channel'])
	imui(CX(:, :, 3), [ud.ImageTitle, ': Cr Channel'])
case '::::mnucb_Merge'
	imuimerge( getImuiFigureHandles )
case '::::mnucb_ToBinary'
	h = waitfig('Converting to Binary image');
	CX = get(ud.ImageHandle, 'UserData');
	imui('::::UpdateImage', im2bw(CX, 0.5))
	delete(h)
case '::::mnucb_ToGray'
	h = waitfig('Converting to Grayscale image');
	switch ud.ImageType
	case 'RGB'
		CX = get(ud.ImageHandle, 'UserData');
		imui('::::UpdateImage', rgb2gray(CX))
	case 'Binary'
		CX = get(ud.ImageHandle, 'UserData');
		imui('::::UpdateImage', ind2gray(CX, [1e-6 1e-6 1e-6; 1 1 1]) )
	end
	delete(h)
case '::::mnucb_ToRGB'
	h = waitfig('Converting to RGB Color image');
	switch ud.ImageType
	case 'Gray'
		CX = get(ud.ImageHandle, 'UserData');
		imui('::::UpdateImage', cat(3, CX, CX, CX) )
	case 'Binary'
		CX = get(ud.ImageHandle, 'UserData');
		CX = ind2gray(CX, [1e-6 1e-6 1e-6; 1 1 1]);
		imui('::::UpdateImage', cat(3, CX, CX, CX) )
	end
	delete(h)
case '::::mnucb_ToUINT8'
	h = waitfig('Converting to UINT8');
	imui('::::UpdateImage',  ...
			im2uint8(get(ud.ImageHandle, 'UserData')) )
	delete(h)
case '::::mnucb_ToUINT16'
	h = waitfig('Converting to UINT16');
	imui('::::UpdateImage',  ...
			im2uint16(get(ud.ImageHandle, 'UserData')) )
	delete(h)
case '::::mnucb_ToDOUBLE'
	h = waitfig('Converting to DOUBLE');
	imui('::::UpdateImage',  ...
			im2double(get(ud.ImageHandle, 'UserData')) )
	delete(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ANALYZE menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_FFT2' % Fourier Transform
	h = waitfig('Computing FFT');
	CX = get(ud.ImageHandle, 'UserData');
	CX = log(abs( fftshift(fft2(CX)) ));
	CX = CX./max(max(CX));
	figure
	imshow(CX)
	colormap(jet)
	colorbar
	delete(h)
case '::::mnucb_DCT'
	h = waitfig('Computing DCT');
	CX = get(ud.ImageHandle, 'UserData');
    figure
    imshow(log(abs(dct2(CX))),[])
    colormap(jet)
    colorbar
	delete(h)
case '::::mnucb_impixel' % ColorPicker
	DisableAllMenus(ud)
	set(ud.FigureHandle, 'Name', ...
			'IMUI - Color picker : click on point(s) and press ENTER.')
	[c, r, P] = impixel;
	if ~isempty(P)
		liststr{1} = '    X     Y        Red  Green   Blue';
		liststr{2} = '--------------------------------------';
		switch ud.ImageClass
		case 'uint8'
			P = double(P);
		case 'uint16'
			P = double(P) ./ 255;
		case 'double'
			P = P * 255;
		end
		P = round(P);
		for i = 1 : length(r)
			liststr{i + 2} = [blanks(0), ...
			sprintf('[%4d, ', c(i)), blanks(0), ...
			sprintf('%4d]', r(i)), blanks(6) ...
			sprintf('[%3d,   %3d,   %3d]', P(i, :)) ...
			];
		end
		hfig = dialog( ...
			'Position',			[0 0 300 250], ...
			'Name',				'Color picker', ...
			'HandleVisibility', 'on', ...
			'WindowStyle',		'modal', ...
			'Visible',			'off');
		hcolor = uicontrol( ...
			'Style',			'text', ...
			'Units',			'pixel', ...
			'Position',			[25 10 250 16], ...
			'String',			'Sample color' ...
			);
		hlst = uicontrol( ...
			'Style',			'listbox', ...
			'Units',			'pixel', ...
			'Position',			[0 35 300 215], ...
			'BackgroundColor',	'w', ...
			'FontName',			'Courier New', ...
			'FontSize',			9, ...
			'String',			liststr, ...
			'UserData',			P ./ 255, ...
			'Callback',			'imui(''::::cb_ColorPicker'')' ...
			);
		movegui(hfig, 'center')
		set(hfig, ...
			'Visible', 			'on', ...
			'UserData',			hcolor)
	end
	EnableAllMenus(ud)
case '::::cb_ColorPicker'
	id = get(gcbo, 'Value');
	if id > 2
		clr = get(gcbo, 'UserData');
		set(get(gcf, 'UserData'), ...
			'BackgroundColor', 		clr(id - 2, :), ...
			'ForegroundColor',		1 - clr(id - 2, :))
	else
		set(get(gcf, 'UserData'), ...
			'BackgroundColor', 		get(0, 'DefaultUIControlBackgroundColor'), ...
			'ForegroundColor',		'k')
	end

case '::::mnucb_improfile'
	DisableAllMenus(ud)
	set(ud.FigureHandle, 'Name', ...
			'IMUI - Profile : drag lines and press ENTER')
	try
		improfile(ud.DefaultOversampleMethod)
	end
	EnableAllMenus(ud)
case '::::mnucb_Contour'
	CX = get(ud.ImageHandle, 'UserData');
	h = waitfig('Drawing contours, this make take a little while.');
	if isrgb(CX)
		CX = rgb2gray(CX);
	end
	h2 = figure('Visible', 'off');
	imcontour(CX);
	set(h2, 'Visible', 'on')
	delete(h);
case '::::mnucb_Histogram'
	imuihist(get(ud.ImageHandle, 'UserData'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADJUST menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_Invert'
	h = waitfig('Inverting');
	imui('::::UpdateImage', ...
		imcomplement(get(ud.ImageHandle, 'UserData')));
	delete(h)
case '::::mnucb_Equalize'
	h = waitfig('Equalizing');
	imui('::::UpdateImage', ...
		histeq(get(ud.ImageHandle, 'UserData')));
	delete(h)
case '::::mnucb_AutoLevels'
	h = waitfig('Applying Auto Levels');
	CX = get(ud.ImageHandle, 'UserData');
	imui('::::UpdateImage', ...
		imadjust( CX, stretchlim(CX) ) );
	delete(h)
case '::::mnucb_Threshold'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuifilter(CX, @imuithreshold);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_BC'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuibc(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_Levels'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuilevels(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_Curves'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuicurves(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_HS'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuihs(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_ColorBalance'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuicbalance(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_Desaturate'
	h = waitfig('Desaturating');
	CX = get(ud.ImageHandle, 'UserData');
	CX = rgb2gray(CX);
	CX = cat(3, CX, CX, CX);
	imui('::::UpdateImage',  CX)
	delete(h)
case '::::mnucb_ChannelMixer'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuichmix(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_Convolution'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuiconv(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
case '::::mnucb_Equalizer'
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuieq(CX);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FILTER menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_FilterChild'
	mnuud = get(gcbo, 'UserData');
	CX = get(ud.ImageHandle, 'UserData');
	CX = imuifilter(CX, mnuud.CustomData);
	if CX ~= -1
		imui('::::UpdateImage',  CX)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VIEW menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_HideMenu'
	if isempty(findobj(ud.FigureHandle, ...
					'Type', 'uimenu', 'Visible', 'on', 'Parent', ud.FigureHandle))
		set(findall(ud.FigureHandle, 'Type', 'uimenu', 'Parent', ud.FigureHandle), ...
			'visible', 'on')
	else
		set(findall(ud.FigureHandle, 'Type', 'uimenu', 'Parent', ud.FigureHandle), ...
			'visible', 'off')

	end
case '::::mnucb_Zoom'
	ud.ZoomState = umtoggle(ud.MenuHandles.Zoom);
	UpdateGUIState(ud)
case '::::mnucb_Grid'
	ud.GridState = umtoggle(ud.MenuHandles.Grid);
	UpdateGUIState(ud)
case '::::mnucb_ZoomIn'
	zoom(2)
	UpdateImageDisplayRatio(ud)
case '::::mnucb_ZoomOut'
	zoom(0.5)
	UpdateImageDisplayRatio(ud)
case '::::mnucb_ZoomFit'
	sz = get(0, 'ScreenSize');
	set(ud.FigureHandle, 'Position', [0 0 sz(3) * 0.8 sz(4) * 0.8])
	movegui(ud.FigureHandle, 'center')
	UpdateImageDisplayRatio(ud)
case '::::mnucb_ZoomFig'
	zoom('out')
	UpdateImageDisplayRatio(ud)
case '::::mnucb_ZoomTrue'
	zoom('out')
	truesize
	UpdateImageDisplayRatio(ud)
case '::::mnucb_GridColor'
	set(allchild(ud.MenuHandles.GridClr), 'Checked', 'off');
	umtoggle(gcbo);
	udmnu = get(gcbo, 'UserData');
	if ~isempty(udmnu.CustomData)
		ud.GridColor = udmnu.CustomData;
	else
		ud.GridColor = uisetcolor(ud.GridColor);
	end
	UpdateGUIState(ud)
case '::::mnucb_FigureColor'
	set(allchild(ud.MenuHandles.FigureClr), 'Checked', 'off');
	umtoggle(gcbo);
	udmnu = get(gcbo, 'UserData');
	if ~isempty(udmnu.CustomData)
		ud.FigureColor = udmnu.CustomData;
	else
		ud.FigureColor = uisetcolor(ud.FigureColor);
	end
	UpdateGUIState(ud)
case '::::mnucb_GridLineStyle'
	set(allchild(ud.MenuHandles.GridStyle), 'Checked', 'off');
	umtoggle(gcbo);
	udmnu = get(gcbo, 'UserData');
	ud.GridLineStyle = udmnu.CustomData;
	UpdateGUIState(ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WINDOW menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_Window'
	imuiFigures = sort(getImuiFigureHandles);
	hsubmenu = allchild( ud.MenuHandles.Window );

	uimenu( ud.MenuHandles.Window, ...
		'Label',		'Center &screen', ...
		'Callback',		'movegui(gcf, ''center'')');
	uimenu( ud.MenuHandles.Window, ...
		'Label',		'&Re-arrange', ...
		'UserData',		imuiFigures, ...
		'Callback',		'imui(''::::mnucb_Window_ReArrange'')');
	uimenu( ud.MenuHandles.Window, ...
		'Label',		'&Close all', ...
		'UserData',		imuiFigures, ...
		'Callback',		'close(get(gcbo, ''UserData''))');

	for i = 1 : length(imuiFigures)

		udfig = get(imuiFigures(i), 'UserData');
		h = uimenu( ud.MenuHandles.Window, ...
			'Label',		['&', num2str(i - 1), ' - ', udfig.ImageTitle], ...
			'UserData',		imuiFigures(i), ...
			'Callback',		'figure(get(gcbo,''UserData''))');
		if i == 1
			set(h, 'Separator', 'on')
		end
		if imuiFigures(i) == ud.FigureHandle
			set(h, 'Checked', 'on')
		end
	end
	delete(hsubmenu)
case '::::mnucb_Window_ReArrange'
	imuiFigures = get(gcbo, 'UserData');
	for i = 1 : length(imuiFigures)
		figure(imuiFigures(i))
		movegui(imuiFigures(i), [i ,-i] * 24)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ABOUT menu callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '::::mnucb_About'
	AboutString = sprintf([ ...
			'Thank you very much for trying IMUI !\n', ...
			'\n', ...
			'If you have suggestions, please write to\n', ...
			'kotaimen_c@citiz.net\n', ...
			'\n', ...
			'Concept, GUI Design, M Programming, DIP Tech...\n', ...
			'All by Kotaimen.C\n', ...
			'\n', ...
			'Special Thanks to:\n', ...
			'Xavier\n', ...
			'\n', ...
			'version 0.9.5\n', ...
			'' ...
		]);

	h = dialog( ...
		'Name',			'About IMUI', ...
		'Position',		[0 0 546 275], ...
		'Color',		'k', ...
		'HandleVisibility', 'on', ...
		'Visible',		'off');
	axes( ...
		'Units',		'pixel', ...
		'Position',		[10 10 256 256]);
	logoimg = fullfile(MyLocation, 'private', 'imuilogo.jpg');
	imshow(logoimg, 'notruesize')
	uicontrol( ...
		'Style',		'text', ...
		'BackgroundColor',	'k', ...
		'ForegroundColor',	[0.8 0.8 0.8], ...
		'Units',		'pixel', ...
		'Position',		[276 50 256 256 - 45], ...
		'String',		AboutString);
	uicontrol( ...
		'Style',		'pushbutton', ...
		'BackgroundColor',	'k', ...
		'ForegroundColor',	[0.8 0.8 0.8], ...
		'Units',		'pixel', ...
		'Position',		[367 15 74 24], ...
		'String',		'Close', ...
		'Callback',		'close(gcf)');
	movegui(h, 'center')
	set(h, 'Visible', 'on')
otherwise
	warning(['Inavid action "', Action, '".'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UpdateGUIState - update gui state
function UpdateGUIState(ud)
%%%%Update figure userdata
set(ud.FigureHandle, 'UserData', ud)

%%%%Update Title
hdlMenu = ud.MenuHandles;
switch ud.ImageType
	% cImgType
	% [Reserved, Binary, Indexed, Gray, RGB]
case 'Binary'
	strTitle = ['IMUI - [', ud.ImageTitle, '] - Binary, '];
	cImgType = [-1, 1, -1, -1, -1];
case 'Gray'
	strTitle = ['IMUI - [', ud.ImageTitle, '] - Grayscale, '];
	cImgType = [-1, -1, -1, 1, -1];
case 'RGB'
	strTitle = ['IMUI - [', ud.ImageTitle, '] - RGB Color, '];
	cImgType = [-1, -1, -1, -1, 1];
otherwise
	error('Unknown ImageType')
end

strTitle = [strTitle, sprintf('%d*%d, ', ud.ImageSize), ...
	upper(ud.ImageClass), ' @'];
set(ud.FigureHandle, ...
	'Name',		strTitle)
%%%%Update Menu
temp = struct2cell(ud.MenuHandles);
allMenuHandles = cat(2, temp{:});
for hmnu = allMenuHandles
	udmnu = get(hmnu, 'UserData');
	if any(cImgType == udmnu.AvailableImageType)
		set(hmnu, 'Enable', 'on')
	else
		set(hmnu, 'Enable', 'off')
	end
end
%%%%Update ImageClass
switch ud.ImageClass
case 'uint8'
	set(ud.MenuHandles.ToUINT8, 'Enable', 'off')
case 'uint16'
	set(ud.MenuHandles.ToUINT16, 'Enable', 'off')
case 'double'
	set(ud.MenuHandles.ToDOUBLE, 'Enable', 'off')
end
%%%%Update Fade menu
if ud.HistoryCount == 1
	set(ud.MenuHandles.Fade, 'Enable', 'off')
end
%%%%Update zoom state
if ud.ZoomState
	setptr(ud.FigureHandle, 'glass')
	set(ud.MenuHandles.Crop, 'Enable', 'off')
	set(ud.MenuHandles.impixel, 'Enable', 'off')
	set(ud.MenuHandles.improfile, 'Enable', 'off')
	zoom('on');
else
	setptr(ud.FigureHandle, 'arrow')
	zoom('off');
end
%%%%Update figure state
	set(ud.FigureHandle, ...
		'Color',		ud.FigureColor);

%%%%Update Gridline state
if ud.GridState
	set(ud.AxesHandle, ...
		'Visible',		'on', ...
		'Color',		ud.FigureColor, ...
		'Position',		[0.1 0.1 0.8 0.8], ...
		'XColor',		ud.GridColor, ...
		'YColor',		ud.GridColor, ...
		'GridLineStyle',ud.GridLineStyle)
	set(ud.TextHandle, ...
		'Visible',		'on', ...
		'BackgroundColor', ud.FigureColor, ...
		'ForegroundColor', ud.GridColor, ...
		'String',		'')
%	set(ud.FigureHandle, ...
%		'WindowButtonMotionFcn', 'imui( ''::::ButtonMotion'') ')
else
	set(ud.AxesHandle, ...
		'Visible',		'off', ...
		'Position',		[0 0 1 1])
	set(ud.TextHandle, ...
		'Visible',		'off')
%	set(ud.FigureHandle, ...
%		'WindowButtonMotionFcn', ' ')
end
UpdateImageDisplayRatio(ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imuiFigHdl -returns handles of imui figure
function imuiFigHdl = getImuiFigureHandles

imuiFigHdl =( findobj(allchild(0), 'tag', '::::imuiFigure', 'Visible', 'on') )';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateImageDisplayRatio(ud)
oldname = get(ud.FigureHandle, 'Name');
lastat = findstr(oldname, '@');
if ud.ZoomState
	newname = [oldname(1 : lastat(end)-1), '@Zoom'];
else
	figpos = get(ud.FigureHandle, 'Position');
	axepos = get(ud.AxesHandle, 'Position');
	axelimx = diff(get(ud.AxesHandle, 'Xlim'));
	axelimy = diff(get(ud.AxesHandle, 'Ylim'));
	temp = figpos .* axepos;
	axex = temp(3);
	axey = temp(4);
	if axex > axey
		ZoomRatio = axex / axelimx;
	else
		ZoomRatio = axey / axelimy;
	end
	newname = [oldname(1 : lastat(end)-1), sprintf('@%.1f%%', ZoomRatio * 100)];
end
set(ud.FigureHandle, 'Name', newname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DisableAllMenus(ud)
temp = struct2cell(ud.MenuHandles);
allMenuHandles = cat(2, temp{:});
set(allMenuHandles, 'Enable', 'off')
set(ud.FigureHandle, ...
	'KeyPressFcn',		' ', ...
	'CloseRequestFcn',	' ', ...
	'WindowButtonMotionFcn', ' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EnableAllMenus(ud)
set(ud.FigureHandle, ...
	'KeyPressFcn',		' imui(  ''::::KeyPress''  ) ', ...
	'CloseRequestFcn',	' imui(  ''::::CloseRequest''  ) ', ...
	'WindowButtonMotionFcn', 'imui( ''::::ButtonMotion'') ')
UpdateGUIState(ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MyLoaction - returns directory where imui.m locates
function R = MyLocation()
W = which('imui.m');
R = W(1 : strfind(W, 'imui.m') - 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DisplayFunction(hdlImg)
CIMG = get(hdlImg, 'UserData');
CIMG = im2uint8(CIMG);
if isbw(CIMG)
	CIMG = grayxform(CIMG, [0, ones(1, 255)]);
	set(hdlImg, 'CData', CIMG)
	return
end
if isrgb(CIMG)
	set(hdlImg, 'CData', CIMG)
	return
end
if isgray(CIMG)
	set(hdlImg, 'CData', CIMG)
	return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AppendToHistory(ud, ActionName)
newHistoryCount = ud.HistoryCount + 1;
newHistoryIndex = newHistoryCount;

hmnu = uimenu(ud.MenuHandles.History, ...
	'Label', 		['&', num2str(newHistoryCount), ' - ', ...
						ActionName], ...
	'UserData',		newHistoryCount, ...
	'Callback',		'imui(''::::mnucb_History_ChildCbFCN'')');

set(allchild(ud.MenuHandles.History), 'Checked', 'off')
set(hmnu, 'Checked', 'on')
if newHistoryIndex ~= ud.HistoryIndex + 1;
	set(hmnu, 'Separator', 'on')
end

ud.LastHistoryIndex = ud.HistoryIndex;
ud.HistoryCount = newHistoryCount;
ud.HistoryIndex = newHistoryIndex;
ud.HistoryData(newHistoryCount).ActionName = [num2str(newHistoryCount), ' - ', ActionName];
TempFileName = tempname;
ud.HistoryData(newHistoryCount).FileName = TempFileName;

set(ud.FigureHandle, 'UserData', ud)
assignin('caller', 'ud', ud)
CX = get(ud.ImageHandle, 'UserData');
save(TempFileName, 'CX');