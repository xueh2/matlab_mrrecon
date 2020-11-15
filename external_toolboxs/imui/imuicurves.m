function CY = imuicurves(varargin)
%IMUICURVES Adjust image level using curves
%    Member of IMUI
%    Kotaimen.C, 2002/05 - 2002/07, All Rights Reserved


error(nargchk(1, 1, nargin))

if ~ischar(varargin{1})
	Action = '::::BuildGUI';
	CX = varargin{1};
	CPREV = thumb(CX);
else
	Action = varargin{1};
	ud = get(gcf, 'UserData');
end

switch Action
case '::::BuildGUI'
	ud.fig = dialog( ...
		'Position',			[0 0 605 370], ...
		'BackingStore',		'on', ...
		'Name',				'Curves', ...
		'HandleVisibility', 'on', ...
		'WindowStyle',		'modal', ...
		'Interruptible', 	'off', ...
		'Busyaction', 		'queue', ...
		'Visible',			'off');
	ud.curves = axes( ...
		'Units',			'pixel', ...
		'Position',			[45 40 300 300], ...
		'DrawMode',			'fast', ...
		'Box',				'on', ...
		'XLim',				[0 1], ...
		'YLim',				[0 1], ...
		'XGrid',			'on', ...
		'YGrid',			'on', ...
		'FontSize',			8, ...
		'XTick',			linspace(0, 1, 11), ...
		'YTick',			linspace(0, 1, 11));
	xlabel('\bfInput');ylabel('\bfOutput');
	title('\bf\fontsize{9}Gray Level Transfer Function');
	line([0 1], [0 1], ...
		'LineStyle', 		':', ...
		'Color',		'k');
	ud.linecurve = line( ...
		'XData',			linspace(0, 1, 256), ...
		'YData',			linspace(0, 1, 256), ...
		'LineWidth',		2, ...
		'Color',			[0 0.5 1], ...
		'EraseMode',		'xor');
	ud.linectl = line( ...
		'XData',			linspace(0, 1, 3), ...
		'YData',			linspace(0, 1, 3), ...
		'LineStyle',		'none', ...
		'Marker', 			'o', ...
		'MarkerSize',		8, ...
		'MarkerEdgeColor', 	[0 0 1], ...
		'MarkerFaceColor',	'none', ...
		'EraseMode',		'xor');
	ud.info  = uicontrol( ...
		'Style',			'text', ...
		'Units',			'pixel', ...
		'Position',			[360 250 85 90], ...
		'FontName',			'Courier New', ...
		'FontSize',			9, ...
		'HorizontalAlignment', 'left', ...
		'String',			'');
	ud.move = uicontrol( ...
		'Style',			'toggle', ...
		'Units',			'pixel', ...
		'Position',			[360 210 25 25], ...
		'TooltipString',	'Move control point [Keyboard "1"]', ...
		'String',			'M', ...
		'Value',			1, ...
		'Callback',			'imuicurves(''::::cb_toggles'')');
	ud.add = uicontrol( ...
		'Style',			'toggle', ...
		'Units',			'pixel', ...
		'Position',			[385 210 25 25], ...
		'TooltipString',	'Add control point [Keyboard "2"]', ...
		'String',			'+', ...
		'Callback',			'imuicurves(''::::cb_toggles'')');
	ud.del = uicontrol( ...
		'Style',			'toggle', ...
		'Units',			'pixel', ...
		'Position',			[410 210 25 25], ...
		'TooltipString',	'Delete control point [Keyboard "3"]', ...
		'String',			'-', ...
		'Callback',			'imuicurves(''::::cb_toggles'')');
	ud.reset  = uicontrol( ...
		'Style',			'pushbutton', ...
		'Units',			'pixel', ...
		'Position',			[360 175 75 24 ], ...
		'String',			'Reset', ...
		'Callback',			'imuicurves(''::::cb_reset'')');
	[J, T] = histeq(CX);
	ud.equalize  = uicontrol( ...
		'Style',			'pushbutton', ...
		'Units',			'pixel', ...
		'Position',			[360 145 75 24 ], ...
		'String',			'Equalize', ...
		'Callback',			'imuicurves(''::::cb_equalize'')', ...
		'UserData',			T);
	ud.pop  = uicontrol( ...
		'Style',			'popupmenu', ...
		'Units',			'pixel', ...
		'Position',			[360 105 75 24 ], ...
		'BackgroundColor',	'w', ...
		'String',			{'Nearest', 'Linear', 'Spline', 'PCHIP'}, ...
		'Value',			3, ...
		'Callback',			'imuicurves(''::::cb_pop'')');
	ud.import  = uicontrol( ...
		'Style',			'pushbutton', ...
		'Units',			'pixel', ...
		'Position',			[360 70 75 24 ], ...
		'String',			'Import...', ...
		'Callback',			'imuicurves(''::::cb_import'')');
	ud.export  = uicontrol( ...
		'Style',			'pushbutton', ...
		'Units',			'pixel', ...
		'Position',			[360 40 75 24 ], ...
		'String',			'Export...', ...
		'Callback',			'imuicurves(''::::cb_export'')');

	ud.preview = axes( ...
		'Units',			'pixel', ...
		'Position',			[460 210 128 128]);

	ud.img = imshow(CPREV, 'notruesize');
	set(ud.preview, ...
		'Visible',			'off', ...
		'DrawMode',			'fast')
	set(ud.img, ...
		'UserData', 		CPREV, ...
		'EraseMode',		'xor')
	ud.chk  = uicontrol( ...
		'Style',			'check', ...
		'Units',			'pixel', ...
		'Position',			[496 185 64 24 ], ...
		'Value',			1, ...
		'String',			'Preview', ...
		'Callback',			'imuicurves(''::::cb_chk'')');
	ud.apply  = uicontrol( ...
		'Style',			'pushbutton', ...
		'Units',			'pixel', ...
		'Position',			[488 140 75 24 ], ...
		'FontWeight',		'bold', ...
		'String',			'Apply', ...
		'Callback',			'imuicurves(''::::cb_apply'')');
	ud.cancel  = uicontrol( ...
		'Style',			'pushbutton', ...
		'Units',			'pixel', ...
		'Position',			[488 105 75 24 ], ...
		'FontWeight',		'bold', ...
		'String',			'Cancel', ...
		'Callback',			'imuicurves(''::::cb_cancel'')');

	ud.HotPoint = -1;
	ud.MouseState = 'release';
	movegui(ud.fig, 'center')
	set(ud.fig, ...
		'Visible',				'on', ...
		'UserData',				ud, ...
		'WindowButtonMotion', 	'imuicurves(''::::MouseMotion'')', ...
		'WindowButtonDownFcn',	'imuicurves(''::::MouseDown'')', ...
		'WindowButtonUpFcn',	'imuicurves(''::::MouseUp'')', ...
		'KeyPressFcn',			'imuicurves(''::::KeyPress'')', ...
		'CloseRequestFcn',		'imuicurves(''::::cb_cancel'')' );

	waitfor(ud.fig, 'Visible', 'off')

	if strcmp('apply', get(ud.apply, 'UserData'))
		h = waitfig('Applying curves');
		CY = grayxform(CX, get(ud.linecurve, 'YData'));
		delete(h)
	else
		CY = -1;
	end
	delete(ud.fig)

%=====================================================================
case '::::KeyPress'
	switch get(ud.fig, 'CurrentCharacter')
	case '1'
		set(ud.move,'Value', 1)
		set(ud.add, 'Value', 0)
		set(ud.del, 'Value', 0)
	case '2'
		set(ud.move,'Value', 0)
		set(ud.add, 'Value', 1)
		set(ud.del, 'Value', 0)
	case '3'
		set(ud.move,'Value', 0)
		set(ud.add, 'Value', 0)
		set(ud.del, 'Value', 1)
	end
%=====================================================================
case '::::MouseDown'
	ud.MouseState = 'down';
	set(ud.fig, 'UserData', ud)
	Sensitivity = 0.02;

	temp = get(ud.curves, 'CurrentPoint');

	CPX = temp(1, 1);
	CPY = temp(1, 2);

	CTX = get(ud.linectl, 'XData');
	CTY = get(ud.linectl, 'YData');

	if get(ud.move,'Value') == 1
		mtype = 'move';
	elseif get(ud.add, 'Value') == 1
		mtype = 'add';
	elseif get(ud.del, 'Value') == 1
		mtype = 'del';
	end

	if all([CPX >= 0, CPX <= 1, CPY >=0, CPY <= 1])
		HTX = [];
		HTY = [];
		for i = 1 : length(CTX)
			if all([abs(CTX(i) - CPX) < Sensitivity, ...
						abs(CTY(i) - CPY) < Sensitivity])
				HTX = CTX(i);
				HTY = CTY(i);
				break
			end
		end

		if isempty(HTX)
			ud.HotPoint = -1;
		else
			ud.HotPoint = i;
		end
		switch mtype
		case 'move'
			if ~isempty(HTX)
				setptr(ud.fig, 'closedhand')
			end
		case 'add'
			HTX = [];
			for i = 1 : length(CTX)
				if abs(CTX(i) - CPX) < Sensitivity
					HTX = CTX(i);
					break
				end
			end
			if length(CTX) >= 15
				errordlg('You have already 15 control points.', 'imui', 'modal')
			elseif isempty(HTX)
				if CPX < CTX(1)
					CTX = [CPX, CTX];
					CTY = [CPY, CTY];
				elseif CPX > CTX(end)
					CTX = [CTX, CPX];
					CTY = [CTY, CPY];
				else
					for i = 1 : length(CTX)
						if CPX > CTX(i) & CPX < CTX(i + 1)
							break
						end
					end
					CTX = [CTX(1 : i), CPX, CTX(i + 1 : end)];
					CTY = [CTY(1 : i), CPY, CTY(i + 1 : end)];
				end
				set(ud.linectl, ...
					'XData',	CTX, ...
					'YData', 	CTY)
				UpdateCurve(ud)
			end
		case 'del'
			if ~isempty(HTX)
				if length(CTX) > 2
					CTX = [CTX(1 : i - 1), CTX(i + 1 : end)];
					CTY = [CTY(1 : i - 1), CTY(i + 1 : end)];
					set(ud.linectl, ...
						'XData',	CTX, ...
						'YData', 	CTY)
					UpdateCurve(ud)
				else
					errordlg('You must have at least 2 control points.', 'imui', 'modal')
				end
			end
		end
		set(ud.fig, 'UserData', ud)
	end
%=====================================================================
case '::::MouseUp'
	ud.MouseState = 'release';
	ud.HotPoint = -1;
	set(ud.fig, 'UserData', ud)
	UpdatePreview(ud)
%=====================================================================
case '::::MouseMotion'
	Sensitivity = 0.02;

	temp = get(ud.curves, 'CurrentPoint');

	CPX = temp(1, 1);
	CPY = temp(1, 2);

	CTX = get(ud.linectl, 'XData');
	CTY = get(ud.linectl, 'YData');

	if get(ud.move,'Value') == 1
		mtype = 'move';
	elseif get(ud.add, 'Value') == 1
		mtype = 'add';
	elseif get(ud.del, 'Value') == 1
		mtype = 'del';
	end

	if all([CPX >= 0, CPX <= 1, CPY >=0, CPY <= 1])
		HTX = [];
		HTY = [];
		for i = 1 : length(CTX)
			if all([abs(CTX(i) - CPX) < Sensitivity, ...
						abs(CTY(i) - CPY) < Sensitivity])
				HTX = CTX(i);
				HTY = CTY(i);
				break
			end
		end
		switch mtype
		case 'move'
			LTY =get(ud.linecurve, 'YData');

			if ud.HotPoint == -1
				set(ud.info, 'String', ...
					sprintf('    X:%.2f\n    Y:%.2f\n\n Input:%d\nOutput:%.d', ...
						CPX, CPY, round(CPX*255), round(LTY(round(CPX*255)+ 1)*255) ) )
					setptr(ud.fig, 'hand')
			else
				set(ud.info, 'String', ...
					sprintf('    X:%.2f\n    Y:%.2f\n\n Input:%d\nOutput:%.d\n\nMoving:%d#', ...
						CPX, CPY, round(CPX*255), round(LTY(round(CPX*255)+ 1)*255),ud.HotPoint ) )
				if strcmp(ud.MouseState, 'down')
					setptr(ud.fig, 'closedhand')
					i = ud.HotPoint;
					if ud.HotPoint == 1
					elseif CPX < CTX(i - 1) + Sensitivity
						CPX = CTX(i - 1) + Sensitivity;
					end
					if i == length(CTX)
					elseif CPX >= CTX(i + 1) - Sensitivity
						CPX = CTX(i + 1) - Sensitivity;
					end
					CTX = [CTX(1 : i - 1), CPX, CTX(i + 1 : end)];

					CTY = [CTY(1 : i - 1), CPY, CTY(i + 1 : end)];
					set(ud.linectl, ...
						'XData',	CTX, ...
						'YData', 	CTY)
					UpdateCurve(ud)
				else
					setptr(ud.fig, 'hand')
				end
			end
		case 'add'
			set(ud.info, 'String', ...
				sprintf('    X:%.2f\n    Y:%.2f', CPX, CPY) )
			HTX = [];
			for i = 1 : length(CTX)
				if abs(CTX(i) - CPX) < Sensitivity
					HTX = CTX(i);
					break
				end
			end
			if isempty(HTX)
				setptr(ud.fig, 'add')
			else
				setptr(ud.fig, 'addzero')
			end
		case 'del'
			if isempty(HTX)
				set(ud.info, 'String', ...
					sprintf('    X:%.2f\n    Y:%.2f', CPX, CPY) )
				setptr(ud.fig, 'eraser')
			else
				set(ud.info, 'String', ...
					sprintf('    X:%.2f\n    Y:%.2f\n\nDeleting:%d#', CPX, CPY, i) )
				setptr(ud.fig, 'addpole')
			end
		end
	else
 		set(ud.info, 'String', '')
		setptr(ud.fig, 'arrow')
	end
%=====================================================================
case '::::cb_toggles'
	set(ud.move,'Value', 0)
	set(ud.add, 'Value', 0)
	set(ud.del, 'Value', 0)
	set(gcbo, 	'Value', 1)
%=====================================================================
case '::::cb_pop'
	UpdateCurve(ud)
	UpdatePreview(ud)
%=====================================================================
case '::::cb_reset'
	set(ud.linectl, ...
		'XData',	[0 0.5 1], ...
		'YData', 	[0 0.5 1])
	set(ud.linecurve, ...
		'XData',	linspace(0 ,1, 256), ...
		'YData', 	linspace(0 ,1, 256))
	set(ud.img, 'CData', get(ud.img, 'UserData'))
%=====================================================================
case '::::cb_chk'
	if get(ud.chk, 'Value')
		UpdatePreview(ud)
	else
		set(ud.img, 'CData', get(ud.img, 'UserData'))
	end
%=====================================================================
case '::::cb_equalize'
	histeqdata = get(ud.equalize, 'UserData');
	cpdata = interp1(linspace(0 ,1, 256), histeqdata, linspace(0 ,1, 15), 'spline');
	set(ud.linectl, ...
		'XData', linspace(0, 1, 15), ...
		'YData', cpdata)
	set(ud.pop, 'Value', 3)
	UpdateCurve(ud)
	UpdatePreview(ud)
%=====================================================================
case '::::cb_import'
	prompt = {sprintf(['Enter varaible name or expression:\n', ...
			'Transfer function must be a 1*256 vector with 0~1 values.\n'])};
	default = {'linspace(0, 1, 256)'};
	dlgtitle = 'Import Transfer Function from Workspace';
	answer = inputdlg(prompt, dlgtitle, 1, default);
	if ~isempty(answer)
		try
			val = evalin('base', answer{1});
			if isstr(val)
				val = str2num(val);
			end
			if all(size(val) == [1, 256]) & all(max(val) <= 1) & all(min(val) >= 0)
				histeqdata = val;
				cpdata = interp1(linspace(0 ,1, 256), histeqdata, linspace(0 ,1, 15), 'spline');
				set(ud.linectl, ...
					'XData', linspace(0, 1, 15), ...
					'YData', cpdata)
				set(ud.pop, 'Value', 3)
				UpdateCurve(ud)
				UpdatePreview(ud)
			else
				h = errordlg('Transfer function must be a 1*256 vector with 0~1 values.', ...
					'Curves', 'modal' );
			end
		catch
			prompt = sprintf(['Error while evaluatng:\n\\bf', ...
						answer{1} , '\\rm. \n\n', lasterr]);
			h = errordlg(prompt, 'Curves', ...
				struct('Interpreter', 'tex', 'WindowStyle', 'modal') );
		end

	end
%=====================================================================
case '::::cb_export'
	prompt = {sprintf('Enter varaible name :\n')};
	default = {'gtf'};
	dlgtitle = 'Export Transfer Function to Workspace';
	answer = inputdlg(prompt, dlgtitle, 1, default);
	if ~isempty(answer)
		varName = answer{1};
		if isvarname(varName)
			assignin('base', varName, get(ud.linecurve, 'YData'));
		else
			prompt = sprintf(['\\bf', varName , '\\rm is not a valid varaible name. \n\n', ...
			'A valid variable name is a character string of letters, digits and', ...
    		'underscores, with length <= 31 and the first character a letter.']);
			h = errordlg(prompt, 'Curves', ...
				struct('Interpreter', 'tex', 'WindowStyle', 'modal') );
		end
	end
%=====================================================================
case '::::cb_apply'
	set(ud.apply, 'Userdata', 'apply')
	set(ud.fig, 'Visible', 'off')
%=====================================================================
case '::::cb_cancel'
	set(ud.apply, 'Userdata', 'cancel')
	set(ud.fig, 'Visible', 'off')
end
%=====================================================================
function UpdateCurve(ud)

CTX = get(ud.linectl, 'XData');
CTY = get(ud.linectl, 'YData');
CTX = [-1, CTX, 2];
CTY = [0, CTY, 1];
FPX = linspace(0 ,1, 256);
switch get(ud.pop, 'Value');
case 1
	FPY = interp1(CTX, CTY, linspace(0 ,1, 256), 'nearest');
case 2
	FPY = interp1(CTX, CTY, linspace(0 ,1, 256), 'linear');
case 3
	FPY = spline(CTX, CTY, linspace(0 ,1, 256));
case 4
	FPY = pchip(CTX, CTY, linspace(0 ,1, 256));
end

for i = 1 : 256
	if i/256 < CTX(2)
		FPY(i) = CTY(2);
	end
	if i/256 > CTX(end - 1)
		FPY(i) = CTY(end - 1);
	end
	if FPY(i) < 0
		FPY(i) = 0;
	end
	if FPY(i) > 1
		FPY(i) = 1;
	end
end

set(ud.linecurve, ...
	'XData',	FPX, ...
	'YData', 	FPY)
%=====================================================================
function UpdatePreview(ud)
if get(ud.chk, 'Value')
	CY = grayxform(get(ud.img, 'UserData'), get(ud.linecurve, 'YData'));
	set(ud.img, 'CData', CY)
end	