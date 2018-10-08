function imagescn(I, scale, dims, FigureWidth, timedimension, roi_name, FigureLeft, FigureBottom)
% function imagescn(I,[min max],[rows cols],FigureWidth,timedimension, roi_name)
%
% function to display multiple images
%	I can be 2-d, i.e., single image
%            3-d,       array of images
%            4-d,       2-d array of images
%            5-d,       3-d array of images
%
%   timedimension specifies which dimension is used for time in "movie" mode
%        and activates the movie mode button. if time dimension is not
%        specified then there are 5 dimensions and the last dimension defaults to the time
%        dimension. the time dimension must be between 3 and number of
%        dimensions (max=5).
%
% user specified scale [min max] or is applied to all images (i.e., clims)
% or scale=[min1 max1; min2 max2; ...] will be applied to images 1,2,...,etc.
% defaults to independent scaling of each sub-image to full range (for scale=[])
%
% sub-images are displayed in N1 rows by N2 columns
%   N1 and N2 may be specified in function input as "rows" and "cols"
%		if input values for rows & cols input are such that
%		rows*cols < full number of images, then the 1-st rows*cols images are displayed
%	or N1 and N2 may be calculated by default:
%		for 3-d array (N1 and N2 are calculated to make approx square image array)
%		for 4-d array (N1 and N2 are 3rd and 4th array dimensions, respectively)
%
% FigureWidth sets the width in inches (defaults to 6 inches). It also sets the paper width
% so that exported figures (to files) will have the specified width.
%
% usage:  imagescn(I)
%         imagescn(I,[],[],[],[])
%         imagescn(I,scale)
%         imagescn(I,[],[rows cols])
%         imagescn(I,scale,[rows cols])
%         imagescn(I,[],[],[],timedimension)
%         ...

% written by: 	Peter Kellman  (kellman@nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%           	Dan Herzka  (herzkad@nhlbi.nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%
Nd=ndims(I);
if ~exist('timedimension')& Nd==5; timedimension=5; moviemodeflag=1; end 
if ~exist('timedimension'); timedimension=[]; moviemodeflag=0; end

if(nargin<6) 
    roi_name = []; 
end

if(nargin<7) 
    FigureLeft = -1; 
end
if(nargin<8) 
    FigureBottom = -1; 
end

if nargin>=5
    if isempty(timedimension); moviemodeflag=0; else; moviemodeflag=1;end
    if timedimension>Nd | timedimension<=2; moviemodeflag=0;end
    if moviemodeflag==1
        if Nd==4 & timedimension==3
            I=permute(I,[1 2 4 3]);
        end
        if Nd==5 & timedimension==3
            I=permute(I,[1 2 4 5 3]);
        end
        if Nd==5 & timedimension==4
            I=permute(I,[1 2 3 5 4]);
        end
    end
end
ntimeframes=size(I,Nd);
if Nd==2; % case of single image
    N=1;
	N1=1; N2=1;
elseif Nd==3; % case of array of images
    if moviemodeflag==0
        N=size(I,3);
		N2=ceil(sqrt(N)); N1=ceil(N/N2);
	else
        N=1;N1=1;N2=1;
	end
elseif Nd==4; % case of 2-d array of images
    if moviemodeflag==0
        N1=size(I,3);
        N2=size(I,4);
		N=N1*N2;
	else
        N=size(I,3);
		N2=ceil(sqrt(N)); N1=ceil(N/N2);
	end
elseif moviemodeflag==1 & Nd==5; % case of 2-d array of images
    N1=size(I,3);
    N2=size(I,4);
	N=N1*N2;
end
if exist('dims');
	if length(dims)==2; rows=dims(1);cols=dims(2);
		N1=rows;N2=cols;
	else
		if ~isempty(dims);disp('Error: must enter [rows cols] for dimensions'); return;end
	end
end
if ~exist('scale','var'); scale=[];end
if 1<size(scale,1) & size(scale,1) < min(N,N1*N2)
    disp('scale vector must be either: empty, size 1x2, i.e. [min max],')
    disp('or a matrix [min1 max1;min2 max2;...] with # rows equal to number of plots');
    return
end
                
set(0,'Units','Inches');
scnsize=get(0,'ScreenSize'); % [left,bottom,width,height] % full screen size available
ScreenWidth=scnsize(3); ScreenHeight=scnsize(4); % width & height in inches
Xsize=size(I,2);Ysize=size(I,1); % size of pixels in image (Xsize x Ysize)
border_percent=.01;
deltaX =(border_percent*Xsize); deltaY=(border_percent*Ysize); % calculate pixel borders as % of size
X0size=N2*Xsize+(N2-1)*deltaX; Y0size=N1*Ysize+(N1-1)*deltaY; % full figure size in pixels (before display)
aspect_ratio=Y0size/X0size; % aspect ratio

% center figure on screen with specified figure width in inches
if ~exist('FigureWidth'); FigureWidth=6;end  % figure width in inches (default is 6")
if isempty(FigureWidth); FigureWidth=6;end
FigureHeight=FigureWidth*aspect_ratio; % figure height in inches

if(FigureBottom<0)
    FigureBottom=(ScreenHeight-FigureHeight)/2;
end

if(FigureLeft<0)
    FigureLeft=(ScreenWidth-FigureWidth)/2;
end

fig_handle=gcf;%figure;
set(fig_handle,'Units','Inches')
pos_f = get(fig_handle, 'Position');
set(fig_handle,'Position',[FigureLeft FigureBottom FigureWidth FigureHeight])
% set(fig_handle,'Position',[pos_f(1) pos_f(2) FigureWidth FigureHeight])

% calculate sub-image dimensions in inches
SubImageWidth=FigureWidth*Xsize/X0size;
SubImageHeight=FigureHeight*Ysize/Y0size;
Xborder=FigureWidth*deltaX/X0size;
Yborder=FigureHeight*deltaY/Y0size;

% set background color to be white
set(fig_handle,'Color',[1 1 1]);

% temporary fix for v7 pre-release bug
% v=version;
% if v(1)=='7';
%     set(fig_handle,'Color','none','Renderer','Painters','InvertHardCopy','off');
% end

% calculate sub-image dimensions in normalized units
SubImageWidth=SubImageWidth/FigureWidth;
SubImageHeight=SubImageHeight/FigureHeight;
Xborder=Xborder/FigureWidth;
Yborder=Yborder/FigureHeight;

for k=1:min(N,N1*N2)
	i=ceil(k/N2);
	j=k-(i-1)*N2;
    if Nd>3
    	i0=ceil(k/size(I,4));
		j0=k-(i0-1)*size(I,4);
    end        
	SubImageLeft=(j-1)*(SubImageWidth+Xborder);
	SubImageBottom=(N1-i)*(SubImageHeight+Yborder);
	subplot('Position',[SubImageLeft SubImageBottom SubImageWidth SubImageHeight])
% 	if isempty(scale)
% 		imagesc(I(:,:,k)); axis image; axis off;
%     elseif size(scale,1)==1
% 		imagesc(I(:,:,k),scale); axis image;axis off
%     elseif size(scale,1)==min(N,N1*N2);
%         imagesc(I(:,:,k),scale(k,:)); axis image;axis off
% 	end
	if isempty(scale)
		if Nd>4
            %imagesc(I(:,:,i0,j0)); axis image; axis off;
            imshow(I(:,:,i0,j0), []); axis image; axis off;
		else 		
            imshow(I(:,:,k), []); axis image; axis off;
		end
	elseif size(scale,1)==1
		if Nd>4
            imshow(I(:,:,i0,j0),scale); axis image; axis off;
		else		
            imshow(I(:,:,k),scale); axis image;axis off
		end
	elseif size(scale,1)==min(N,N1*N2);
		if Nd>4
            imshow(I(:,:,i0,j0),scale(k,:)); axis image; axis off;
		else		
            imshow(I(:,:,k),scale(k,:)); axis image;axis off
		end
	end
    if moviemodeflag==1        
        if Nd==3
        	setappdata(gca, 'ImageData', I);
        elseif Nd==4
        	setappdata(gca, 'ImageData', squeeze(I(:,:,k,:)));
        elseif Nd==5
            setappdata(gca, 'ImageData', squeeze(I(:,:,i0,j0,:)));
        end
		setappdata(gca, 'ImageRange', [1 ntimeframes]);
		setappdata(gca, 'ImageRangeAll', [1 ntimeframes]);
		setappdata(gca, 'CurrentImage', 1);
    end
end

% set(fig_handle, 'PaperPosition', [0 0 FigureWidth FigureHeight]); % old
set(fig_handle, 'PaperPosition', [1 1 FigureWidth+1 FigureHeight+1]);

toolbox;
colormap(gray)

if(exist('roi_name') & ~isempty(roi_name))
    Load_ROI(roi_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Change_Current_Axes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the current axes when user clicks on an axes
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_axes= gca;
data_holder = findobj('Tag','ROI_Title_text');
data=  get(data_holder, 'Userdata');
data{4} = h_axes;
set(data_holder,'UserData',data);

% now highlight current axes in the information windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Change_Current_ROI(ROI_info)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the current ROI when user clicks on an active object
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_ROI_circle= gco;
data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');

if nargin==0
    vals = get(h_ROI_circle, 'Userdata');
	data_holder{5} = vals(1:2);
	data_holder{4} = get(h_ROI_circle, 'Parent');
else
   data_holder{5} =  ROI_info;
   UserData = get(findobj(data_holder{1}, 'Tag', 'figROITool'), 'UserData');
   ROI_info_table = UserData{1};
   h_circle = ROI_info_table(ROI_info(1),ROI_info(2)).ROI_Data(1,1);
   data_holder{4} = get(h_circle, 'Parent');
end;
set(findobj('Tag', 'ROI_Title_text'), 'UserData', data_holder);

% now select the current ROI in the information windows
Highlight_Current_ROI(data_holder{5});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Angle_Adjust_Entry;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute once at the beggining of a drag cycle
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

fig = gcf;
set(fig, 'WindowButtonMotionFcn', 'ROI_tool(''ROI_Angle_Adjust'');');
set(fig,'WindowButtonUpFcn', 'ROI_tool(''ROI_Angle_Adjust_Exit'')');
Change_Current_Axes;
ROI_Angle_Adjust;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Angle_Adjust
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_angle = gco;
h_axes = get(h_angle, 'Parent');
point = get(h_axes,'CurrentPoint');

%CALCULATION

val = get(h_angle, 'Userdata');
% 1,2  infor table indexes
% 3-7  handle_values = [h_circle, h_center, h_size, h_angle, h_number]
% 8-12 ROI_values = [center_x, center_y, size_x, size_y, angle]

v1 = [point(1,1) point(1,2)]  -[ val(8), val(9)] ;
v2 = [get(h_angle, 'xdata'), get(h_angle, 'ydata')] - [val(8), val(9)];
% get angle (positive only) and multiply it by direction
d = cross([v1 0],[v2 0]);
% calculate angle between the two...
alpha=  acos(  dot(v1,v2)   /(norm(v1) * norm(v2)) ) * -1*sign(d(3));
%alpha_deg = alpha*180/pi
rotmat = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

% now rotate everything by this amount
c = rotmat*[get(val(3),'xdata') - val(8); get(val(3),'ydata') - val(9)];
set(val(3), 'xdata', c(1,:) + val(8), 'ydata', c(2,:) + val(9));
c = rotmat*[get(val(5),'xdata') - val(8); get(val(5),'ydata') - val(9)];
set(val(5), 'xdata', c(1,:) + val(8), 'ydata', c(2,:) + val(9));
c = rotmat*[get(val(6),'xdata') - val(8); get(val(6),'ydata') -  val(9)];
set(val(6), 'xdata', c(1,:)+ val(8), 'ydata', c(2,:)+ val(9));
p = get(val(7),'Position')';
c = rotmat*( [p(1:2) - [val(8) val(9)]']) ;
set(val(7), 'Position', [ c' + [val(8),val(9)] ,p(3)]);

% update only the current object's userdata... 
%set(val(6),'UserData', [val(1:11), alpha]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Angle_Adjust_Exit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;
h_angle = gco;
val = get(h_angle, 'Userdata');

% 1,2  infor table indexes
% 3-7  handle_values = [h_circle, h_center, h_size, h_angle, h_number]
% 8-12 ROI_values = [center_x, center_y, size_x, size_y, angle]

%CALCULATION
v1 = [1 0];

if isappdata(h_angle, 'alpha0')
	alpha0 = getappdata(h_angle, 'alpha0');
else
	alpha0 = 0;
end;

v2 = [get(h_angle, 'xdata'), get(h_angle, 'ydata')] - [val(8), val(9)] ;
d = cross([v1 0],[v2 0]);
alpha = acos(  dot(v1,v2)   /(norm(v1) * norm(v2)) ) *sign(d(3));
alpha = alpha - alpha0;
%alpha_deg = alpha*180/pi

% update all other objects part of this ROI with correct values
set(val(3:7), 'userdata', [val(1:11), alpha]);

data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
fig = data_holder{1};
fig2 = data_holder{2};
ifig = data_holder{6};
h_all_axes = data_holder{3};
handles = guidata(fig2);
apply_all = get(handles.Link_ROI_togglebutton,'Value');
userdata = get(findobj(fig, 'Tag', 'figROITool'),'Userdata');
ROI_info_table = userdata{1};
update_list = [];

if apply_all
    for i = 1:length(h_all_axes(find(h_all_axes)))        
%        if ~isempty(ROI_info_table(val(1),i).ROI_Exists) & ROI_info_table(val(1),i).ROI_Exists
        if ROI_info_table(val(1),i).ROI_Exists
			
            set(ROI_info_table(val(1),i).ROI_Data(1,1), ...
                'xdata', get(val(3),'xdata') ,...
                'ydata', get(val(3),'ydata'))
            set(ROI_info_table(val(1),i).ROI_Data(1,2), ...
                'xdata', get(val(4),'xdata') ,...
                'ydata', get(val(4),'ydata'))
            set(ROI_info_table(val(1),i).ROI_Data(1,3), ...
                'xdata', get(val(5),'xdata') ,...
                'ydata', get(val(5),'ydata'))
            set(ROI_info_table(val(1),i).ROI_Data(1,4), ...
                'xdata', get(val(6),'xdata') ,...
                'ydata', get(val(6),'ydata'))
            set(ROI_info_table(val(1),i).ROI_Data(1,5), ...
                'position', get(val(7), 'Position'));
            for k = 1:5
                old_val  = get(ROI_info_table(val(1),i).ROI_Data(1,k), 'Userdata');
                set(ROI_info_table(val(1),i).ROI_Data(1,k),...
                    'Userdata', [old_val(1:7), val(8:end)]);  
            end;
            ROI_info_table(val(1),i).ROI_Data(2,1:5) =  val(8:end);
			
			if ~isempty (ROI_info_table(val(1),i).ROI_x_original)
				xdata = get(val(3),'xdata')';
				ydata = get(val(3),'ydata')';
				% downsample the new data to create "original" data
				ROI_info_table(val(1),i).ROI_x_original = xdata(1:Sample_Rate:end);
				ROI_info_table(val(1),i).ROI_y_original = ydata(1:Sample_Rate:end);				
			end;
			
            update_list(size(update_list,1)+1,:) = [val(1), i]; 
        end;
    end;
else
    update_list = [val(1:2)];
	if ~isempty (ROI_info_table(val(1),val(2)).ROI_x_original)
		xdata = get(val(3),'xdata')';
		ydata = get(val(3),'ydata')';
		% downsample the new data to create "original" data
		ROI_info_table(val(1),val(2)).ROI_x_original = xdata(1:Sample_Rate:end);
		ROI_info_table(val(1),val(2)).ROI_y_original = ydata(1:Sample_Rate:end);				
	end;

end;

userdata{1} = ROI_info_table;
set(findobj(fig, 'Tag', 'figROITool'),'Userdata',userdata);

% set the current ROI in storage
Change_Current_ROI(val(1:2));
% update the ROI info in the ROI_info_Table
Update_ROI_Info(update_list);
% Update the info into the 3D string holder
Update_ROI_Info_String(update_list);
% manipulate 3D string into a page for display in listbox
Resort_ROI_Info_Listbox;
% bring the ROI info figures to the front
figure(fig2);
figure(ifig);
figure(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Size_Adjust_Entry;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute once at the beggining of a drag cycle
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

fig = gcf;
set(fig, 'WindowButtonMotionFcn', 'ROI_tool(''ROI_Size_Adjust'');');
set(fig,'WindowButtonUpFcn', 'ROI_tool(''ROI_Size_Adjust_Exit'')');
Change_Current_Axes;
ROI_Size_Adjust;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Size_Adjust
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% global STORAGE_VAR

h_size = gco;
h_axes = get(h_size, 'Parent');
point = get(h_axes,'CurrentPoint');
val = get(h_size, 'Userdata');
% 1,2  info table indexes
% 3-7  handle_values = [h_circle, h_center, h_size, h_angle, h_number]
% 8-12 ROI_values = [center_x, center_y, size_x, size_y, angle]

%CALCULATION
%val(8:12)
%alpha_deg = val(12)*180/pi

% get the new position
center_pt = [val(8) val(9)];
v1 = [point(1,1)           point(1,2)]           - center_pt;   % new position
v2 = [get(val(5),'xdata'), get(val(5),'ydata')]  - center_pt ;  % old position



% Undo rotation v1 -> inverse v1 - val(12)= angle from +x-axis to circle marker 
alpha = val(12);
%alpha = alpha;
rotmat =[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
% undo rotation by alpha on the new point
iv1 = rotmat'*v1';

% d = cross([v2 0],[v1 0])     % angle between old size marker and new size marker                              % after return to horizontal plane
% theta=  acos(  dot(v1,v2)   /(norm(v1) * norm(v2)) ) *sign(d(3));
% theta_deg = theta*180/pi

% determine the skew and use the ratio of the old skew value to the new
% skew value to determine the new coordinates
sx2 = -(iv1(1));
sy2 = -(iv1(2));

%temp_skewmat = [sx2/val(10) 0 ; 0 sy2/val(11)];
skewmat = rotmat * [sx2/val(10) 0 ; 0 sy2/val(11)] * rotmat';

c = skewmat * [get(val(3),'xdata') - center_pt(1); get(val(3), 'ydata') - center_pt(2)];
set(val(3), 'xdata', c(1,:) + center_pt(1), 'ydata', c(2,:) + center_pt(2));

c = (skewmat*[ [get(val(5),'xdata'), get(val(5),'ydata')] - center_pt]')' + center_pt;
set(val(5), 'xdata', c(1), 'ydata', c(2), 'UserData', [val(1:9), sx2, sy2, val(12)]);
c = (skewmat*[  [ get(val(6),'xdata'),get(val(6), 'ydata') ] - center_pt]')'  + center_pt;
set(val(6), 'xdata', c(1), 'ydata', c(2));
p = get(val(7),'Position')';

c = skewmat* ( p(1:2) - center_pt') + center_pt';
set(val(7), 'Position',    [c ;p(3)]);

set(h_size,'Userdata', [val(1:9), sx2, sy2 ,val(12)]);


% DEBUG Trajectory Storage
% STORAGE_VAR(end+1).alpha = alpha_deg;
% STORAGE_VAR(end).v1 = v1;
% STORAGE_VAR(end).v2 = v2;
% STORAGE_VAR(end).rotmat = rotmat;
% STORAGE_VAR(end).iv1 = iv1;
% STORAGE_VAR(end).sx2 = sx2;
% STORAGE_VAR(end).sy2 = sy2;
% STORAGE_VAR(end).temp_skewmat = temp_skewmat;
% STORAGE_VAR(end).sx2 = sx2;
% STORAGE_VAR(end).sy2 = sy2;
% STORAGE_VAR(end).sx2 = sx2;
% STORAGE_VAR(end).sy2 = sy2;
% STORAGE_VAR(end).sx2 = sx2;
% STORAGE_VAR(end).sy2 = sy2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Size_Adjust_Exit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;


h_size = gco;
val = get(h_size, 'Userdata');
set(val(3:7), 'userdata', val);

data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
fig = data_holder{1};
fig2 = data_holder{2};
ifig = data_holder{6};
h_all_axes = data_holder{3};
handles = guidata(fig2);
apply_all = get(handles.Link_ROI_togglebutton,'Value');

userdata = get(findobj(fig, 'Tag', 'figROITool'),'Userdata');
ROI_info_table = userdata{1};
update_list = [];

if apply_all
    for i = 1:length(h_all_axes(find(h_all_axes)))       
        if ~isempty(ROI_info_table(val(1),i).ROI_Exists) & ROI_info_table(val(1),i).ROI_Exists
            set(ROI_info_table(val(1),i).ROI_Data(1,1), ...
                'xdata', get(val(3),'xdata') ,...
                'ydata', get(val(3),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,2), ...
                'xdata', get(val(4),'xdata') ,...
                'ydata', get(val(4),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,3), ...
                'xdata', get(val(5),'xdata') ,...
                'ydata', get(val(5),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,4), ...
                'xdata', get(val(6),'xdata') ,...
                'ydata', get(val(6),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,5), ...
                'position', get(val(7), 'Position'));
            for k = 1:5
                old_val  = get(ROI_info_table(val(1),i).ROI_Data(1,k), 'Userdata');
                set(ROI_info_table(val(1),i).ROI_Data(1,k),...
                    'Userdata', [old_val(1:7), val(8:end)]); 
            end; 
            ROI_info_table(val(1),i).ROI_Data(2,1:5) =  val(8:end);
			
			if ~isempty (ROI_info_table(val(1),i).ROI_x_original)
				xdata = get(val(3),'xdata')';
				ydata = get(val(3),'ydata')';
				% downsample the new data to create "original" data
				ROI_info_table(val(1),i).ROI_x_original = xdata(1:Sample_Rate:end);
				ROI_info_table(val(1),i).ROI_y_original = ydata(1:Sample_Rate:end);				
			end;
			
            update_list(size(update_list,1)+1,:) = [val(1), i];       
        end;
    end;
else
    update_list = [val(1:2)];
	
	if ~isempty(ROI_info_table(val(1),val(2)).ROI_x_original)
		xdata = get(val(3),'xdata')';
		ydata = get(val(3),'ydata')';
		% downsample the new data to create "original" data
		ROI_info_table(val(1),val(2)).ROI_x_original = xdata(1:Sample_Rate:end);
		ROI_info_table(val(1),val(2)).ROI_y_original = ydata(1:Sample_Rate:end);
	end;
	
	
end;

userdata{1} = ROI_info_table;
set(findobj(fig, 'Tag', 'figROITool'),'Userdata',userdata);

Change_Current_ROI(val(1:2));
Update_ROI_Info(update_list);
Update_ROI_Info_String(update_list);
Resort_ROI_Info_Listbox;
figure(fig2);
figure(ifig);
figure(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Pos_Adjust_Entry(origin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute once at the beggining of a drag cycle
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;


fig = gcf;
set(fig, 'WindowButtonMotionFcn', ['ROI_tool(''ROI_Pos_Adjust'',' num2str(origin) ');']);
set(fig,'WindowButtonUpFcn', 'ROI_tool(''ROI_Pos_Adjust_Exit'')');
Change_Current_Axes;
ROI_Pos_Adjust(origin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Pos_Adjust(origin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_pos = gco;
h_axes = get(h_pos, 'Parent');
point = get(h_axes,'CurrentPoint');
val = get(h_pos, 'Userdata');

%CALCULATION


% 1,2  infor table indexes
% 3-7  handle_values = [h_circle, h_center, h_size, h_angle, h_number]
% 8-12 ROI_values = [center_x, center_y, size_x, size_y, angle]

p = get(val(7),'Position')';
% center transform  = new point - old _point
if origin ==1  % call from center plus
    new_center_pt = [point(1,1), point(1,2)] - [val(8) val(9)];
elseif origin ==2  % call from corner number
    new_center_pt= [point(1,1), point(1,2)] - [p(1) p(2)];
end;

d = [get(val(4),'xdata') , get(val(4), 'ydata')] + new_center_pt;
set(val(4), 'xdata',  d(1), 'ydata', d(2));

o = [get(val(3),'xdata')  ;  get(val(3), 'ydata')];
set(val(3), 'xdata', o(1,:) + new_center_pt(1), 'ydata', o(2,:) + new_center_pt(2));
c = [get(val(5),'xdata'), get(val(5),'ydata')] +  new_center_pt;
set(val(5), 'xdata', c(1), 'ydata', c(2));
c = [get(val(6),'xdata'),get(val(6), 'ydata')] + new_center_pt;
set(val(6), 'xdata', c(1), 'ydata', c(2));
c = p(1:2) + new_center_pt';
set(val(7), 'Position',  [c ;p(3)]);

% update info in both number and center
set([val(4), val(7)], 'Userdata', [val(1:7), d(1:2), val(10:12)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Pos_Adjust_Exit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;


h_pos = gco;

% 1,2  infor table indexes
% 3-7  handle_values = [h_circle, h_center, h_size, h_angle, h_number]
% 8-12 ROI_values = [center_x, center_y, size_x, size_y, angle]

% now use the center position as reference for everyone
val = get(h_pos, 'Userdata');
set(val(3:7), 'userdata', val);

data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
fig = data_holder{1};
fig2 = data_holder{2};
ifig = data_holder{6};

h_all_axes = data_holder{3};
handles = guidata(fig2);
apply_all = get(handles.Link_ROI_togglebutton,'Value');
userdata = get(findobj(fig, 'Tag', 'figROITool'),'Userdata');
ROI_info_table = userdata{1};
update_list = [];
if apply_all
    for i = 1:length(h_all_axes(find(h_all_axes)))       
        if ~isempty(ROI_info_table(val(1),i).ROI_Exists) & ROI_info_table(val(1),i).ROI_Exists
            set(ROI_info_table(val(1),i).ROI_Data(1,1), ...
                'xdata', get(val(3),'xdata') ,...
                'ydata', get(val(3),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,2), ...
                'xdata', get(val(4),'xdata') ,...
                'ydata', get(val(4),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,3), ...
                'xdata', get(val(5),'xdata') ,...
                'ydata', get(val(5),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,4), ...
                'xdata', get(val(6),'xdata') ,...
                'ydata', get(val(6),'ydata'));
            set(ROI_info_table(val(1),i).ROI_Data(1,5), ...
                'position', get(val(7), 'Position'));
            for k = 1:5
                old_val  = get(ROI_info_table(val(1),i).ROI_Data(1,k), 'Userdata');
                set(ROI_info_table(val(1),i).ROI_Data(1,k),...
                    'Userdata', [old_val(1:7), val(8:end)]);  
            end;
            ROI_info_table(val(1),i).ROI_Data(2,1:5) =  val(8:end);
            update_list(size(update_list,1)+1,:) = [val(1), i];    
			
			if ~isempty (ROI_info_table(val(1),i).ROI_x_original)
				xdata = get(val(3),'xdata')';
				ydata = get(val(3),'ydata')';
				% downsample the new data to create "original" data
				ROI_info_table(val(1),i).ROI_x_original = xdata(1:Sample_Rate:end);
				ROI_info_table(val(1),i).ROI_y_original = ydata(1:Sample_Rate:end);				
			end;

			
        end;
    end;
else
	%disp(['Not all'])
    update_list = [val(1:2)];
	if ~isempty(ROI_info_table(val(1),val(2)).ROI_x_original)
		xdata = get(val(3),'xdata')';
		ydata = get(val(3),'ydata')';
		% downsample the new data to create "original" data
		ROI_info_table(val(1),val(2)).ROI_x_original = xdata(1:Sample_Rate:end);
		ROI_info_table(val(1),val(2)).ROI_y_original = ydata(1:Sample_Rate:end);
	end;
end;

userdata{1} = ROI_info_table;
set(findobj(fig, 'Tag', 'figROITool'),'Userdata',userdata);

Change_Current_ROI(val(1:2));
Update_ROI_Info(update_list);
Update_ROI_Info_String(update_list);
Resort_ROI_Info_Listbox;
figure(fig2);
figure(ifig);
figure(fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function handles = Prep_Draw_ROI(mode, h_axes);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_all_axes = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
if isempty(h_axes)
	% find axes of interest handle if it wasn't passed in
	h_axes = h_all_axes{4};
end;

% Define Mode posibilties
DRAW = 0; EDIT = 1;

% define cursors
plus_cursor = zeros(16,16)*NaN;
circ_cursor = plus_cursor;
square_cursor = plus_cursor;
ex_button = plus_cursor;

pw = [  8    24    40    56    72   113   114   115   116   117  123   124   125   126   127  168   184   200   216   232];
pb = [  7     9    23    25    39   41    55    57    71    73   97    98    99   100   101 ...
		107   108   109   110   111  129   130   131   132   133  139   140   141   142   143 ...
		167   169   183   185   199  201   215   217   231   233];
cw = [ ...% inner circle  
		23    24    25    37    38 	42    43    52    60    67 	77    83    93    98   110 ...
		114   120   126   130   142   147   157   163   173   180 	188   197   198   202   203 ...
		215   216   217];
cb = [ ... % outer circle
		7     8     9    21    22  	26    27    36    44    51  61    66    78    82    94  ...
		97   111   113   127   129  143   146   158   162   174  179   189   196   204   213  ...
		214   218   219   231   232   233];
ew = [ 1    15    18    30    35    45    52    60    69    75 ...
		86    90   103   105   120   135   137   150   154   165 ...
		171   180   188   195   205   210   222   225   239 ];
eb = [ 2    14    17    19    29    31    34    36    44    46  ...
		51    53    59    61    68    70    74    76    85    87 ...
		89    91   102   104   106   119   121   134  136   138  ...
		149   151   153   155   164   166   170   172   179   181 ...
		187   189   194   196   204   206   209   211   221   223   226   238];

plus_cursor(pw) = 2;
plus_cursor(pb) = 1;
circ_cursor(cw) = 2;
circ_cursor(cb) = 1;
ex_button(ew) = 1; % make button all black...
ex_button(eb) = 1;

% open figure, copy current axes parameters, and draw ROIs if editing
h_temp_figure = openfig('ROI_draw_figure', 'new');
handles = guihandles(h_temp_figure);
handles.Parent_figure = h_temp_figure;
handles.cursors.plus_cursor   = plus_cursor;
handles.cursors.circ_cursor   = circ_cursor;
% handles.cursors.square_cursor = square_cursor;
handles.ButtonDownStrings{1}  = ['ROI_tool(''ROI_Draw_Entry'',' num2str(h_temp_figure), ');' ];
handles.ButtonDownStrings{2}  = ['ROI_tool(''ROI_Push_Point_Entry'',' num2str(h_temp_figure), ');' ];

% initialize other fields used in itearative commands
handles.Spline = [];
handles.Points  = [];
handles.NewPoints = [];
handles.h_spline = [];
handles.h_pixels = [];
handles.h_circle = [];

handles.Show_Pixels = get(handles.Show_Pixels_checkbox, 'Value');
handles.Close_ROI   = get(handles.Close_ROI_checkbox,   'Value');
handles.Point_Drop_Spacing  = str2num(get(handles.Point_Drop_edit, 'String'));

plus_cursor = plus_cursor - 1;  plus_cursor(isnan(plus_cursor)) = 0.8;
circ_cursor = circ_cursor - 1;  circ_cursor(isnan(circ_cursor)) = 0.8;
% square_cursor = square_cursor - 1;  square_cursor(isnan(square_cursor)) = 0.8;
ex_button = ex_button - 1; ex_button(isnan(ex_button))= 0.8;

set(handles.Draw_pushbutton,  'CData', repmat(plus_cursor  , [ 1 1 3]));
set(handles.Push_pushbutton,  'CData', repmat(circ_cursor  , [ 1 1 3]));
% set(handles.Drag_Erase_pushbutton, 'CData', repmat(square_cursor, [ 1 1 3]));
set(handles.Clear_Points_pushbutton,'CData', repmat(ex_button,    [ 1 1 3]));

set(handles.Spline_checkbox, 'Userdata', handles.Polygon_checkbox);
set(handles.Polygon_checkbox, 'Userdata', handles.Spline_checkbox);
set(handles.Temp_Spline_axes, 'Tag', 'Temp_Spline_axes');

% find the axis values and image values of the current axes
h_image = findobj(h_axes, 'Type','Image');
image_values = get(h_image,'CData');
axes_values = get(h_axes, {'xlim', 'ylim', 'clim'});
h_temp_image = imagesc(image_values);

h_temp_axes = get(h_temp_image, 'Parent');
h_temp_axes = handles.Temp_Spline_axes;
%		set(h_temp_axes, 'Tag', 'Temp_Spline_axes');
axis equal; axis off;
hold on;
set(handles.Temp_Spline_axes,'xlim', axes_values{1},'ylim', axes_values{2}, 'clim', axes_values{3});
colormap(get(get(h_axes,'Parent'), 'Colormap'));

% prepare constants for drawing push-cursor
xlim = axes_values{1};
ylim = axes_values{2};
size_x = xlim(2) - xlim(1);
size_y = ylim(2) - ylim(1);
handles.circle_size = min([size_x, size_y]);
handles.circle_ratio = str2num(get(handles.Push_Radius_edit, 'String'))/100;
handles.Color        = 'r';

set(handles.Parent_figure, 'Tag', 'ROI_Draw_figure', ...
	'Pointer', 'custom', 'PointerShapeCData', handles.cursors.plus_cursor, 'PointerShapeHotSpot', [8 8]);
set(h_temp_image, 'ButtonDownFcn', handles.ButtonDownStrings{1});

% store all this info
guidata(handles.Parent_figure, handles);

% Initialize based on incoming parameter
if mode == EDIT
	% Now draw the current ROIsand setup the points
	Prep_Edit_Current_ROI(handles);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Prep_Edit_Current_ROI(handlesDraw);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to prepare editing of an exisitng ROI
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
fig = data_holder{1};
fig2 = data_holder{2};
ifig = data_holder{6};
i_current_ROI = data_holder{5};
h_current_axes = data_holder{4};

colororder = repmat('rbgymcw',1,4);

% disable the other figures, if they exist
if ~isempty(fig2), 
	[handlesDraw.h_objects, handlesDraw.object_enable_states] = Change_Figure_States('Disable', [fig2, ifig]);
end;

UserData = get(findobj(fig, 'Tag', 'figROITool'), 'UserData');
ROI_info_table = UserData{1};
current_ROI_data = ROI_info_table(i_current_ROI(1), i_current_ROI(2));

x_coors = get(current_ROI_data.ROI_Data(1,1), 'Xdata')';
y_coors = get(current_ROI_data.ROI_Data(1,1), 'Ydata')';

% load unsplined points, if they exist
if ~isempty(current_ROI_data.ROI_x_original)
	%disp(['ROI_Tool: Prep_Edit_Current_ROI: Loading un-splined (original) data']);	
	x_coors =current_ROI_data.ROI_x_original;
	y_coors =current_ROI_data.ROI_y_original;
end	


% remove end point if it is the sams as the first point \
% (avoids problems with moving points and fitting splines)
% Should be made more general (i.e. matrix of distances between points, 
% remove all second points with zero distance
thresh = 0.0001;
if ( abs( x_coors(1)- x_coors(end)) < thresh) & ( abs( y_coors(1)- y_coors(end)) < thresh)
	x_coors = x_coors(1:end-1);
	y_coors = y_coors(1:end-1);
end;
	
handlesDraw.Points = [x_coors, y_coors];
handlesDraw.Color = colororder(i_current_ROI(1));
set(handlesDraw.Parent_figure,     'CurrentAxes', handlesDraw.Temp_Spline_axes);
set(handlesDraw.Done_pushbutton,   'Callback', 'ROI_tool(''Edit_ROI_Finish'',''Done'');');
set(handlesDraw.Cancel_pushbutton, 'Callback', 'ROI_tool(''Edit_ROI_Finish'',''Cancel'');');
set(handlesDraw.Parent_figure,     'CloseRequestFcn', 'ROI_tool(''Edit_ROI_Finish'',''Done'');');
set(0, 'CurrentFigure', handlesDraw.Parent_figure);

for i = 1:length(x_coors)	
	pp = plot(x_coors(i),y_coors(i), 'r.', 'Userdata', [x_coors(i), y_coors(i)]);
	set(pp, 'ButtonDownFcn', ['ROI_tool(''ROI_Point_Move_Entry'',' , num2str(handlesDraw.Parent_figure), ',''' , num2str(pp,20), ''');'],...
		'MarkerEdgeColor', handlesDraw.Color);
	%get(pp)
end;

handlesDraw.EDITING = 1;
guidata(handlesDraw.Parent_figure, handlesDraw);

Draw_Spline;
Show_Pixels;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Draw_ROI_Finish(Mode);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to finish creating Freehand ROIs or exit from editing any ROI
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
stored_data = get(handlesDraw.Parent_figure, 'Userdata');
apply_spline = get(handlesDraw.Spline_checkbox,'Value');

[fig, fig2, ifig, h_axes_interest, h_all_axes, h_axes, colororder, Current_ROI_index, Userdata, ...
		ROI_info_table, first_ROI_flag, handlesROI, h_objects, object_enable_states	]  = deal(stored_data{:});
	
% Enable objects after in Draw Mode
if ~isempty(fig2), 
	Change_Figure_States('Enable', [fig2, ifig], h_objects,object_enable_states);
end;

		
% delete (NOT close) the old figure
delete(handlesDraw.Parent_figure);

pts = handlesDraw.Points;
if isempty(pts) | strcmp(Mode, 'Cancel')
	return
end;	

%CALCULATION


% Close the ROI by copying the final point
x = [pts(:,1); pts(1,1)];  y = [pts(:,2) ; pts(1,2)];
xs = x; ys = y;
if apply_spline
	% 	interpolate x & y points into a spline curve
	[xs, ys] = Spline_ROI(x, y);
end;

center_x = mean(xs);
center_y = mean(ys);

maxx = max(xs); maxx = maxx(1);
maxy = max(ys); maxy = maxy(1);
minx = min(xs); minx = minx(1);
miny = min(ys); miny = miny(1);
imaxx = find(xs==maxx); imaxx = imaxx(end);

size_x1 = maxx - center_x;
size_x2 = abs(minx - center_x);
size_y  = abs(miny - center_y);

angle_point_x = xs(imaxx);
angle_point_y = ys(imaxx);

v1 = [1 0];
v2 = [angle_point_x, angle_point_y] - [center_x, center_y] ;
d = cross([v1 0],[v2 0]);
alpha0 = acos(  dot(v1,v2)   /(norm(v1) * norm(v2)) ) *sign(d(3));
% Initialize the starting angle as alpha
alpha = 0;

for i = 1:length(h_axes_interest(:))
	set(fig, 'CurrentAxes', h_axes_interest(i));
	h_axes_index = find(h_all_axes'==h_axes_interest(i));
	set(0, 'CurrentFigure', fig);

	% handle_values = [h_circle, h_center, h_size, h_angle, h_number];
	handle_values = Make_ROI_Elements(...
		xs, ys,...
		colororder(Current_ROI_index),...
		Current_ROI_index,...
		center_x, center_y,...
		center_x - size_x2, center_y - size_y, ...
		angle_point_x, angle_point_y,...
		center_x + size_x1, center_y - size_y, ...
		alpha0);
			
	
	ROI_values = [center_x, center_y, size_x2, size_y, alpha ];
	
	set(handle_values, 'UserData', ...
		[Current_ROI_index, h_axes_index, handle_values, ROI_values ]);
	ROI_info_table(Current_ROI_index,h_axes_index).ROI_Data = ...
		[handle_values ; ...
			ROI_values];
	ROI_info_table(Current_ROI_index, h_axes_index).ROI_Exists = 1;
	
	if apply_spline
		% If spline was applied, save original points
		% Clear original points if simply drawing straight lines
		ROI_info_table(Current_ROI_index,h_axes_index).ROI_x_original = x;
		ROI_info_table(Current_ROI_index,h_axes_index).ROI_y_original = y;
	else
		ROI_info_table(Current_ROI_index,h_axes_index).ROI_x_original = [];
		ROI_info_table(Current_ROI_index,h_axes_index).ROI_y_original = [];	
	end

	update_list(i,:) = [Current_ROI_index, h_axes_index];	
	i_current_ROI = [Current_ROI_index, h_axes_index];
	
end;			

% Now Restore ROI_info_table to its hiding spot
Userdata{1} = ROI_info_table;
set(findobj(fig, 'Tag', 'figROITool'), 'UserData', Userdata);

% call the ROI_info update function: puts data into ROI_info_table
Update_ROI_Info(update_list);

if first_ROI_flag
    % creates figure the first time and creates the string table that is to be
    % used for "publishing" the ROI data
    ifig = Create_ROI_Info_Figure( update_list);
    % published the string into the listbox
    
    % turn on buttons, but turn off print objects
    Change_Object_Enable_State(handlesROI,'Off',1);
    Change_Object_Enable_State(handlesROI,'On',0);
else
    % call function that will take info string table and "publish" it
    
    %Update_ROI_Info(update_list);
    Update_ROI_Info_String(update_list);
end;

Resort_ROI_Info_Listbox;
% update current ROI index
set(findobj('Tag', 'ROI_Title_text'), 'Userdata', { fig, fig2, h_all_axes, h_axes, i_current_ROI, ifig});
Highlight_Current_ROI(i_current_ROI);

figure(fig2);
figure(ifig);
figure(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Edit_ROI_Finish(Mode);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to insert the edited ROI into the ROI_info table
% MODE = 'DONE' - insert ROI
% MODE = 'CANCEL' - do nothing but close the ROI
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
apply_spline = get(handlesDraw.Spline_checkbox,'Value');
pts = handlesDraw.Points;
colororder = repmat('rbgymcw',1,4);
data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');

% if isempty(h_all_axes)
% 	% Check to see if by any chance the parent figure has been closed
% 	% If it has been closed (as demonstrated by the empty handle), close
% 	% the figure and return;
% 	delete(findobj('Tag', 'ROI_Draw_figure'));
% 	return;
% end;
	
fig = data_holder{1};
fig2 = data_holder{2};
ifig = data_holder{6};
h_axes = data_holder{4};
i_current_ROI = data_holder{5};
h_all_axes = data_holder{3};

handlesROI = guidata(fig2);
apply_all = get(handlesROI.Link_ROI_togglebutton,'Value');

% delete (NOT close) the old figure
delete(handlesDraw.Parent_figure);

% Enable objects after in Draw Mode
if ~isempty(fig2), 
	Change_Figure_States('Enable', [fig2, ifig], handlesDraw.h_objects,handlesDraw.object_enable_states);
end;

if  strcmp(Mode, 'Cancel')
	return;
end;	

pts = handlesDraw.Points;
if isempty(pts) 
	% There are no points - delete the current ROI
	Delete_ROI(1);
	return
end;

UserData = get(findobj(fig, 'Tag', 'figROITool'), 'UserData');
ROI_info_table = UserData{1};
h_current_ROI_axes = get(ROI_info_table(i_current_ROI(1), i_current_ROI(2)).ROI_Data(1,1),'Parent');

%CALCULATION

% Close the ROI by copying the final point
x = [pts(:,1); pts(1,1)];  y = [pts(:,2) ; pts(1,2)];
xs = x; ys = y;
if apply_spline
	% 	interpolate x & y points into a spline curve
	[xs, ys] = Spline_ROI(x,y);
end;

center_x = mean(xs);
center_y = mean(ys);

maxx = max(xs); maxx = maxx(1);
maxy = max(ys); maxy = maxy(1);
minx = min(xs); minx = minx(1);
miny = min(ys); miny = miny(1);
imaxx = find(xs==maxx); imaxx = imaxx(1);

size_x1 = maxx - center_x;
size_x2 = abs(minx - center_x);
size_y  = abs(miny - center_y);

angle_point_x = xs(imaxx);
angle_point_y = ys(imaxx);

v1 = [1 0];
v2 = [angle_point_x, angle_point_y] - [center_x, center_y] ;
d = cross([v1 0],[v2 0]);
alpha = acos(  dot(v1,v2)   /(norm(v1) * norm(v2)) ) *sign(d(3));
alpha = 0;

%%WHY????
h_all_axes = h_all_axes';

if apply_all
	h_axes_interest = h_all_axes(find(h_all_axes));
else
	h_axes_interest = h_axes;
end;

update_list = [];
for i = 1:length(h_axes_interest)
    % for each new ROI created, create the ROI with the same parameters
    %    [h_circle, h_center, h_size, h_angle, h_number; ...
    %        center_x, center_y, size_x, size_y, angle];    

	set(fig, 'CurrentAxes', h_axes_interest(i));
	h_axes_index = find(h_all_axes==h_axes_interest(i));
	set(0, 'CurrentFigure', fig);

	exist_current_ROI = ROI_info_table(i_current_ROI(1), h_axes_index).ROI_Exists;

	% if the ROI exists, erase it before writing a new one in its place!
	if exist_current_ROI==1
		old_handles = ROI_info_table(i_current_ROI(1), h_axes_index).ROI_Data;
		delete(old_handles(1,:));
		
		% handle_values = [h_circle, h_center, h_size, h_angle, h_number];
		handle_values = Make_ROI_Elements(...
			xs,ys,...
			colororder(i_current_ROI(1)),...
			i_current_ROI(1), ...
			center_x, center_y, ...
			center_x - size_x2, center_y - size_y, ...
			angle_point_x, angle_point_y, ...
			center_x + size_x1, center_y - size_y);
		
		ROI_values = [center_x, center_y, size_x2, size_y, alpha];
		
		set(handle_values, 'UserData', ...
			[i_current_ROI(1), h_axes_index, handle_values, ROI_values ]);
		ROI_info_table(i_current_ROI(1),h_axes_index).ROI_Data = ...
			[handle_values ; ...
				ROI_values];
		
		if apply_spline
			% If spline was applied, save original points
			ROI_info_table(i_current_ROI(1),h_axes_index).ROI_x_original = x;
			ROI_info_table(i_current_ROI(1),h_axes_index).ROI_y_original = y;	
		else
			% Clear original points if simply drawing straight lines
			ROI_info_table(i_current_ROI(1),h_axes_index).ROI_x_original = [];
			ROI_info_table(i_current_ROI(1),h_axes_index).ROI_y_original = [];	
		end
		
		ROI_info_table(i_current_ROI(1), h_axes_index).ROI_Exists = 1;
		
		% now add the indexes of the created ROIs to the update list
		update_list(size(update_list,1)+1,:) = [i_current_ROI(1), h_axes_index];
		
		%i_current_ROI = [Current_ROI_index, h_axes_index];
	end;
		
end;

% now restore the modified info table
UserData{1} = ROI_info_table;
set(findobj(fig, 'Tag', 'figROITool'), 'UserData', UserData );

Update_ROI_Info(update_list);
Update_ROI_Info_String(update_list);
Resort_ROI_Info_Listbox;
Highlight_Current_ROI(update_list)

figure(fig2);
figure(ifig);
figure(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Draw_Entry(h_figure);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute once at the beggining of a drag cycle
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

global OLD_POINT; 
OLD_POINT = [-10, -10];  % initialize so that first point is automatically drawn.
set(h_figure, 'WindowButtonMotionFcn', 'ROI_tool(''ROI_Draw'');');
set(h_figure, 'WindowButtonUpFcn',     'ROI_tool(''ROI_Draw_Exit'')'   );
ROI_Draw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Draw;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

global OLD_POINT;
h_axes = gca;
handlesDraw = guidata(get(h_axes, 'Parent'));
point_drop_pix_thresh = str2num(get(handlesDraw.Point_Drop_edit, 'String'));
p = get(h_axes,'CurrentPoint');
p = [ p(1,1), p(1,2) ];
d = sqrt ( sum( (  (OLD_POINT - p).^2)  ,2) );
if  d > point_drop_pix_thresh
	% store point in the new queue and plot
	handlesDraw.NewPoints = [handlesDraw.NewPoints; p ];
	pp = plot(p(1),p(2), 'r.', 'Userdata', p);
	set(pp, 'ButtonDownFcn', ['ROI_tool(''ROI_Point_Move_Entry'',' , num2str(get(h_axes, 'Parent')) , ',''' , num2str(pp,20), ''');'], ...
		'MarkerEdgeColor', handlesDraw.Color);
	OLD_POINT = p;
end;
guidata(get(h_axes, 'Parent'), handlesDraw );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Draw_Exit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;
Sort_Points;
Draw_Spline;
Show_Pixels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Point_Move_Entry(h_figure, h_point);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute once at the beggining of a drag cycle
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

sel_type = get(h_figure, 'SelectionType');
handlesDraw = guidata(h_figure);
if strcmp(sel_type, 'normal')
	% Normal selection = move point around
	handlesDraw.old_WBMF = get(h_figure,'WindowButtonMotionFcn');
	handlesDraw.old_WBUF = get(h_figure,'WindowButtonUpFcn');	
	set(h_figure, 'WindowButtonMotionFcn', ['ROI_tool(''ROI_Point_Move'','''     ,h_point,''');']);
	set(h_figure, 'WindowButtonUpFcn',     ['ROI_tool(''ROI_Point_Move_Exit'',''',h_point,''');']);
	ROI_Point_Move(h_point);
	
	guidata(h_figure, handlesDraw);

else
	% remove the point if the point is double-clicked 
	% or right clicked
	h_point = str2num(h_point);
	pos = cell2mat(get(h_point, {'xdata', 'ydata'}));
	idx = find(sum(abs(handlesDraw.Points - repmat(pos, [size(handlesDraw.Points,1),1])),2)==0);
	delete(h_point);
	% remove point from queue
	for i = 1:length(idx)
		handlesDraw.Points = [handlesDraw.Points(1:idx(i)-1,:); handlesDraw.Points(idx(i)+1:end,:)];
	end;
	% store points again
	% update lines and selected pixels
	guidata(h_figure, handlesDraw);

	Draw_Spline;
	Show_Pixels;
end;
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Point_Move(h_point)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_point = str2num(h_point);
h_axes = get(h_point, 'Parent');
p = get(h_axes,'CurrentPoint');
set(h_point, 'Xdata', p(1,1), 'Ydata', p(1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Point_Move_Exit(h_point)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% Change the position of the point in the cue
handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
h_point = str2num(h_point);
h_axes = get(h_point, 'Parent');

set(handlesDraw.Parent_figure, ...
	'WindowButtonMotionFcn', handlesDraw.old_WBMF , ...
	'WindowButtonUpFcn',     handlesDraw.old_WBUF);
handlesDraw.old_WBMF = [];
handlesDraw.old_WBUF = [];

points = handlesDraw.Points;
old_pos = get(h_point, 'Userdata');
new_pos = cell2mat(get(h_point, {'xdata', 'ydata'}));

set(h_point, 'Userdata', new_pos);
idx = find(sum(abs(points - repmat(old_pos, [size(points,1),1])),2)==0);
%idx = idx(1);
handlesDraw.Points(idx,:) = new_pos;

% now store points again
guidata(findobj('Tag', 'ROI_Draw_figure'),handlesDraw);
Draw_Spline;
Show_Pixels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Draw_Spline;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that draws the fitted spline to a series of point on an ROI
% curve. If spline is not desired, a straight line fit (polygon) is
% created.
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
draw_spline = get(handlesDraw.Spline_checkbox, 'Value');

if ~isempty(handlesDraw.h_spline)
	% delete the old line
	delete(handlesDraw.h_spline);
	handlesDraw.h_spline = [];
	handlesDraw.Spline   = [];
end;
% unsampled points
points = handlesDraw.Points;

if ~isempty(points)
	
	close_ROI = get(handlesDraw.Close_ROI_checkbox, 'Value');
	if close_ROI,
		points = [points; points(1,:)];
	end;
	if  draw_spline & ( size(points,1) > 1)	
		% draw or redraw spline
		set(handlesDraw.Parent_figure, 'CurrentAxes', handlesDraw.Temp_Spline_axes);
		[xs, ys] = Spline_ROI(points(:,1),points(:,2));
		handlesDraw.h_spline = plot3(xs, ys, repmat(-1,size(xs)), [handlesDraw.Color, '-']);
		handlesDraw.Spline  = [xs', ys'];
	else
		% draw straight line
		handlesDraw.h_spline = plot3(points(:,1), points(:,2), repmat(-1,size(points(:,1))), [handlesDraw.Color, '-']);
		handlesDraw.Spline  = points;
	end;	
end;
guidata(findobj('Tag', 'ROI_Draw_figure'), handlesDraw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Toggle_Spline_Poly(h_checkbox)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to toggle between drawing splines and polygons
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

val = get(h_checkbox, 'Value');

if val
	set(get(h_checkbox, 'Userdata'), 'Value', 0);
else
	set(get(h_checkbox, 'Userdata'), 'Value', 1);
end;
Draw_Spline;
Show_Pixels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Toggle_Draw_Mode(Mode)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;
% switch between draw mode and edit-push more
handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));;

if strcmp(Mode, 'Push')	
	set(handlesDraw.Parent_figure, 'PointerShapeCData', handlesDraw.cursors.plus_cursor);
	set(findobj(handlesDraw.Parent_figure, 'Type', 'Image'), 'ButtonDownFcn', handlesDraw.ButtonDownStrings{2});
	set(handlesDraw.Parent_figure, 'WindowButtonMotionFcn', 'ROI_tool(''ROI_Push_Move_Cursor'');');
	
	guidata(findobj('Tag', 'ROI_Draw_figure'),handlesDraw);
	Draw_Push_Radius(1);
	ROI_Push_Move_Cursor;

elseif strcmp(Mode, 'Draw')
	set(handlesDraw.Parent_figure, 'PointerShapeCData', handlesDraw.cursors.plus_cursor);
	set(findobj(handlesDraw.Parent_figure, 'Type', 'Image'), 'ButtonDownFcn', handlesDraw.ButtonDownStrings{1});
	set(handlesDraw.Parent_figure, 'WindowButtonMotionFcn', '');
	
	if ~isempty(handlesDraw.h_circle)
		delete(handlesDraw.h_circle);
		handlesDraw.h_circle = [];
	end;
	guidata(findobj('Tag', 'ROI_Draw_figure'),handlesDraw);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Push_Move_Cursor;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_axes = gca;
p = get(h_axes, 'CurrentPoint');
h_circle = findobj(h_axes, 'Tag', 'Cursor_Circle');
xy = get(h_circle, 'UserData');
set(h_circle, 'XData', xy(:,1) + p(1,1), 'YData', xy(:,2) + p(1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Sort_Points;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to sort the newly added points to the end of the cue
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_axes = gca;
handlesDraw = guidata(get(h_axes, 'Parent'));

old_pts = handlesDraw.Points;
new_pts = handlesDraw.NewPoints;

% Only sort if there is more than two points!
if (size(old_pts,1) > 2) %~isempty(old_pts) &  
	% calculate the distance between the first new_pt the last new pt, all the old points
	% to determine the best entry for the new curves
	dist_start = sqrt ( sum (    (old_pts - repmat(new_pts(1  ,:), [ size(old_pts,1),1])).^2   , 2));
	dist_end   = sqrt ( sum (    (old_pts - repmat(new_pts(end,:), [ size(old_pts,1),1])).^2   , 2));
	
	% determine smallest net distance to 2 adjacent points
	% assume user draws first point closest to the point of entry
	min_dist_idx = find(min(dist_start)==dist_start);
	idx1 = mod( min_dist_idx + 1 - 1 , length(dist_start) ) + 1;
	idx2 = mod( min_dist_idx - 1 - 1,  length(dist_start) ) + 1;
	% assume all new points fit between the first new point
	other_min_dist = min([dist_end([idx1,  idx2])]);
	other_min_dist_idx = find(other_min_dist == dist_end);
		
	if other_min_dist_idx  < min_dist_idx
		% if inserting backwards, flip indexes
		t = other_min_dist_idx;
		other_min_dist_idx = min_dist_idx;
		min_dist_idx = t;
	end;	
	
	
	
	if abs(other_min_dist_idx - min_dist_idx) > 1
		% looking at inserting after the last point
		handlesDraw.Points = [old_pts; new_pts];
	elseif other_min_dist_idx == length(old_pts)
		% if the end point is the closest point, just tack on at the end
		handlesDraw.Points = [old_pts; new_pts];
	else
		%insert between the two indexes
		handlesDraw.Points = [ ...
				old_pts(1:min_dist_idx(1),:); ...
				new_pts; ...
				old_pts(min_dist_idx(1)+1:end,:) ];
	end;
	
else
	% not enough points, tack points on at the end without sorting
	handlesDraw.Points = [old_pts; new_pts];
end;
% Store points and clear new pts
handlesDraw.NewPoints = [];
guidata(get(h_axes, 'Parent'),handlesDraw );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Push_Point_Entry(h_figure);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute once at the beggining of a drag cycle
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% Normal selection = move point around
set(h_figure, 'WindowButtonMotionFcn', ['ROI_tool(''ROI_Push_Point'');']);
set(h_figure, 'WindowButtonUpFcn',     ['ROI_tool(''ROI_Push_Point_Exit'');']);
ROI_Push_Point;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Push_Point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
points = handlesDraw.Points;
p = get(handlesDraw.Temp_Spline_axes,'CurrentPoint');
p = [p(1,1), p(1,2)];
if ~isempty(points)
	% if ROI has been cleared, do nothing
	pixel_thresh = handlesDraw.circle_size * handlesDraw.circle_ratio;
	pixel_push   = pixel_thresh * 1.1;	
	% find the points that are within 4 pixels of the center of the cursos
	dist = sqrt(sum( (points - repmat(p,[size(points,1),1]) ).^2, 2) )';
	in_radius = find(dist <= pixel_thresh);
	for i = 1:length(in_radius)
		% move each point to the edge of a circle of radius 4, in the direction
		% of the vector from the center of the circle to the points
		h_current_point = findobj(handlesDraw.Temp_Spline_axes, ...
			'Xdata', points(in_radius(i),1), 'Ydata', points(in_radius(i),2));
		% unit vector pointing in the correct direction
		v = (points(in_radius(i),:) - p) ./ sqrt ( sum  ( (points(in_radius(i),:) - p).^2  , 2));
		% add unit vector * push_radius to center coordinates
		new_p = p + v.* ( pixel_push );
		set(h_current_point, 'Xdata', new_p(1,1), 'Ydata', new_p(1,2), 'UserData', [new_p(1,1), new_p(1,2)]);
		points(in_radius(i),:) = new_p;
	end;
	% now store the poitns again
	handlesDraw.Points = points;
	guidata(findobj('Tag', 'ROI_Draw_figure'), handlesDraw);
	Draw_Spline;
end;
ROI_Push_Move_Cursor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_Push_Point_Exit;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% Update the Pixels 
handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
set(handlesDraw.Parent_figure, 'WindowButtonMotionFcn', 'ROI_tool(''ROI_Push_Move_Cursor'');');
Show_Pixels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Draw_Clear_ROI;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
all_points = findobj(handlesDraw.Temp_Spline_axes,'Type', 'line')';
if ~isempty(handlesDraw.h_circle)
	% don't delete the circular cursor!!!
	i = find(all_points==handlesDraw.h_circle);
	if ~isempty(i)
		all_points = [all_points(1:i-1), all_points(i+1:end)];
	end
end;
delete(all_points);	
handlesDraw.h_spline = [];
handlesDraw.Spline = [];
handlesDraw.Points = [];
handlesDraw.h_pixels = [];
guidata(findobj('Tag', 'ROI_Draw_figure'), handlesDraw );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Show_Pixels;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
handlesDraw.Show_Pixels = get(handlesDraw.Show_Pixels_checkbox, 'Value');

if ~isempty(handlesDraw.h_pixels)
	delete(handlesDraw.h_pixels);
	handlesDraw.h_pixels = [];
end;

if handlesDraw.Show_Pixels & ~isempty(handlesDraw.Points)
	if ~isempty(handlesDraw.Spline)
		xpts = handlesDraw.Spline(:,1);
		ypts = handlesDraw.Spline(:,2);
	else
		xpts = handlesDraw.Points(:,1);
		ypts = handlesDraw.Points(:,2);
	end;
	im = get(findobj(handlesDraw.Temp_Spline_axes,'Type', 'Image'),'CData');		
	% 	xx = repmat([1:size(im,2)], size(im,1),1);
	% 	yy = repmat([1:size(im,1)]',1,          size(im,2));
	[xx,yy] = meshgrid([1:size(im,2)],[1:size(im,1)]);
	rr = inpolygon(xx,yy,xpts,ypts);
	[ii,jj] = find(rr);
	h_pixels = plot3(jj,ii,repmat(-1, size(jj)), 'bs' , 'MarkerSize', 2, 'MarkerFaceColor', 'b');
	handlesDraw.h_pixels = h_pixels;
end;
guidata(findobj('Tag', 'ROI_Draw_figure'), handlesDraw );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Draw_Change_Edit_Value(h_edit,Limits, Default_value)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% Function to insure number entered into editable text
% box is valid
new_val = get(h_edit, 'String');
reject = 1;
try
	x = str2num(new_val);
	if isnumeric(x)
		if (x>Limits(1)) & (x<Limits(2))
			% accept;
			reject = 0;
		end;
	end;
end;
if reject
	set(h_edit, 'String', num2str(Default_value));
end;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Draw_Push_Radius(Mode);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% Function to (re)draw the cursor for the push tool
% if Mode == 1, then first time in, and draw cursor
% if Mode == 0, then already in push mode and redraw
handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
basic_points = 64;	
theta = 0:(360/basic_points):360;
[x,y] = pol2cart(theta*pi/180, repmat(handlesDraw.circle_size * handlesDraw.circle_ratio , size(theta,1), size(theta,2)));

if ~isempty(handlesDraw.h_circle) | (Mode==1)
	delete(handlesDraw.h_circle);
	handlesDraw.h_circle = plot(x-100 ,y-100, 'w-', 'Tag', 'Cursor_Circle', 'Userdata', [x',y']); 
	guidata(findobj('Tag', 'ROI_Draw_figure'), handlesDraw );
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Draw_Change_Radius_Value(h_edit);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

handlesDraw = guidata(findobj('Tag', 'ROI_Draw_figure'));
handlesDraw.circle_ratio = str2num(get(handlesDraw.Push_Radius_edit, 'String'))/100;
guidata(findobj('Tag', 'ROI_Draw_figure'), handlesDraw );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Update_ROI_Info(update_list)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updates the roi info (mean, std, pixs, etc) into the ROI table
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

if ~isempty(update_list)
    h_all_axes = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
    fig = h_all_axes{1};    
    fig2 = h_all_axes{2};   
    h_all_axes = h_all_axes{3};
    handles = guidata(fig2);
    userdata = get(findobj(fig, 'Tag', 'figROITool'),'Userdata');
    ROI_info_table = userdata{1};

	debug_mode = 0;
       
    for i = 1:size(update_list,1)
        % get the handles of the image, and the ROI circle
        h_circle = ROI_info_table(update_list(i,1), update_list(i,2)).ROI_Data(1);
        xpts = get(h_circle, 'xdata');
        ypts = get(h_circle, 'ydata');
        im = get(findobj(get(h_circle, 'Parent'), 'Type', 'Image'), 'CData');

        % check boundary conditions
        %xpts(xpts<=1) = 1;  xpts(xpts>size(im,1)) = size(im,1);
        %ypts(ypts<=1) = 1;  ypts(ypts>size(im,2)) = size(im,2);
        
        % reduce the matrix size        
        min_xpts = min(xpts); max_xpts = max(xpts);
        min_ypts = min(ypts); max_ypts = max(ypts);
         
        % shift indexes
        xpts2 = xpts - floor(min_xpts) ;
        ypts2 = ypts - floor(min_ypts) ;
       
        %reduce image size too cover only points
        im2 = im(floor(min_ypts):ceil(max_ypts), floor(min_xpts):ceil(max_xpts));
        %figure; imagesc(im2)
        %hold on; plot(xpts2, ypts2, 'ro-');
        %axis image
        
        %xx = repmat([1:size(im2,2)], size(im2,1),1);
        %yy = repmat([1:size(im2,1)]',1, size(im2,2));
		  [xx,yy] = meshgrid([1:size(im,2)],[1:size(im2,1)]);

        % Do not use roipoly as it only uses integer vertex coordinates
        %BW = roipoly(im, xpts, ypts);
        % However, because in_polygon uses vector and cross products to determine
        % if point is within polygon, make matrix smaller.
        %tic
        rr = inpolygon(xx,yy,xpts2,ypts2);
        %toc
        [ii,jj] = find(rr);

		if debug_mode
            %plot(jj+floor(min_xpts),ii+floor(min_ypts),'r.');
            f = figure;
            imagesc(im2);
            axis equal; 
            hold on;
			title('Press Any Key To Continue');
			plot(jj+1,ii+1,'r.');
            plot(xpts2+1,ypts2+1,'r-.')
			colormap(get(fig, 'Colormap'));
            pause
			try, close(f); end;
        end;
        ii = ii + 1; jj = jj + 1;        
        ROI_vals = double(im2(sub2ind(size(im2),ii,jj)));
   
        mn  = mean(ROI_vals);
        stdev= std(ROI_vals);
        mins = min(ROI_vals);
        maxs = max(ROI_vals);
        pixels = length(ROI_vals);
            
        ROI_info_table(update_list(i,1), update_list(i,2)).ROI_Info = ...
            [mn, stdev, pixels, mins, maxs];    
    end;
    % restore the ROI_info_table with its new info
    userdata{1} = ROI_info_table;    
    set(findobj(fig, 'Tag','figROITool'), 'UserData', userdata);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Listbox_Change_Current_ROI;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;
h_listbox = gcbo;

% now determine the index of the current ROI
str = get(h_listbox, {'Value', 'String'});
values = str{1}; str = str{2};

if  length(values) >1
    % too many things are highlighted, highlight only last one
    values = values(end);
    set(h_listbox,'Value', values);
end;
% avoid problems is string is empty as all ROIs are deleted 
% and there can't be a current ROI
if ~isempty(str)
    % take first 8 characters, and convert to two numbers 
    ROI_info = fliplr(str2num(str(values,1:8)));
    h_data_holder = findobj('Tag','ROI_Title_text');
    data=  get(h_data_holder, 'Userdata');
    data{5} = ROI_info;
	
	% also update the current axes;
	UserData = get(findobj(data{1}, 'Tag', 'figROITool'), 'UserData');
	ROI_info_table = UserData{1};
	h_circle = ROI_info_table(ROI_info(1), ROI_info(2)).ROI_Data(1,1);
	data{4} = get(h_circle, 'Parent');
	
	set(h_data_holder,'UserData',data);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Resort_ROI_Info_Listbox(h_pushbutton);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to reassemble the String_info_table into a page 
% for display in the listbox; selects the current ROI;
% if called by the toggle button, handle is sent in, if 
% called after the addition of data to the _table,
% then no handle is sent in.
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

if nargin == 0 
    h_pushbutton = findobj('Tag', 'Sort_Order_pushbutton');
end;
a = get(h_pushbutton, {'Value','Userdata', 'String', 'Parent'} );
sort_order = a{1}; new_str= a{2}; cur_str = a{3}; ifig = a{4};
% sort_order = 1 = sort by image, 2 = sort by ROI

% if call from button, toggle button string
if nargin ==1, set(h_pushbutton, 'String', new_str, 'UserData', cur_str); end;

handles = guidata(ifig);
String_info_table = get(handles.ROI_Info_listbox, 'Userdata');
h_data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
fig = h_data_holder{1};    
userdata = get(findobj(fig, 'Tag', 'figROITool'), 'Userdata'); 
ROI_info_table = userdata{1};

if strcmp('Image', new_str)
    String_info_table = permute(String_info_table, [2 1 3]);
end;

st = size(String_info_table);
String_info_table = reshape(String_info_table, st(1)*st(2), st(3));

% now deblank empty rows
% assumes last digit in row is not space in normal strings
if (size(String_info_table, 1)>1)
    g = find(String_info_table(:,size(String_info_table,2))'=='x');
    String_info_table = String_info_table(g,:);
end;
set(handles.ROI_Info_listbox,'String', String_info_table(:,1:end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Save_ROI(pname, fname);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subfunction that saves the ROI_info_table and String_update_table
% into .mat file and/or saves the Page displayed in the listbox to a text file
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_data_holder = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
fig = h_data_holder{1};    
fig2 = h_data_holder{2};   
handles = guidata(fig2);
userdata = get(findobj(fig, 'Tag', 'figROITool'),'Userdata');
ROI_info_table = userdata{1};

ROI_string_table = get(findobj('Tag', 'ROI_Info_listbox'), 'String');

save_mat  = get(handles.Save_MAT_checkbox, 'Value');
% save_text = get(handles.Save_TXT_checkbox, 'Value');
save_matlab_var = get(handles.Save_matlab_var_checkbox, 'Value');
 
if save_mat | save_matlab_var
    % get the filename
    if save_mat
        if nargin < 2
            fname = []; pname = [];
            [fname,pname] = uiputfile('*.mat', 'Save .mat file');
        end;
    end
    if save_matlab_var
       prompt={'Enter the matlab variable name for ROI:'};
       name='export ROI variable';
       numlines=1;
       defaultanswer={'roi'};
       roi_variable_name=char(inputdlg(prompt,name,numlines,defaultanswer));        
    end
    
    % check if there is appdata with image
    h_all_axes = flipud(findobj(fig,'Type','Axes'));
    
    for i =1:size(ROI_info_table,1)
        for j = 1:size(ROI_info_table,2)
            if ROI_info_table(i,j).ROI_Exists
                ROI_info_table(i,j).ROI_x_coordinates = ...
                    get(ROI_info_table(i,j).ROI_Data(1), 'xdata');
                ROI_info_table(i,j).ROI_y_coordinates = ...
                    get(ROI_info_table(i,j).ROI_Data(1), 'ydata');
                P = get(ROI_info_table(i,j).ROI_Data(1,5),'Position');
                ROI_info_table(i,j).Other_coordinates = ...
                    [get(ROI_info_table(i,j).ROI_Data(1,2), 'xdata'),...
                        get(ROI_info_table(i,j).ROI_Data(1,2), 'ydata'),...
                        get(ROI_info_table(i,j).ROI_Data(1,3), 'xdata'),...
                        get(ROI_info_table(i,j).ROI_Data(1,3), 'ydata'),...
                        get(ROI_info_table(i,j).ROI_Data(1,4), 'xdata'),...
                        get(ROI_info_table(i,j).ROI_Data(1,4), 'ydata'),...
                        P(1:2),...
                    ];
                if isappdata(h_all_axes(j), 'ImageData')
                    image_data    = getappdata(h_all_axes(j), 'ImageData');
                    BW=zeros(size(image_data(:,:,1)));
                    BW=roipoly(image_data(:,:,1),ROI_info_table(i,j).ROI_x_coordinates,ROI_info_table(i,j).ROI_y_coordinates);
                    index=find(BW >0);
                    for k=1:size(image_data,3)
                        tmp=image_data(:,:,k);
                        tmp=double(tmp);
                        m(k)=mean(tmp(index));
                        s(k)=std(tmp(index));
                        roi_min(k)=min(tmp(index));
                        roi_max(k)=max(tmp(index));
                    end
                    ROI_info_table(i,j).ROI_mean   = m;
                    ROI_info_table(i,j).ROI_stdev  = s;
                    ROI_info_table(i,j).ROI_pixels = ROI_info_table(i,j).ROI_Info(3);
                    ROI_info_table(i,j).ROI_min    = roi_min;
                    ROI_info_table(i,j).ROI_max    = roi_max;
                else
                    ROI_info_table(i,j).ROI_mean   = ROI_info_table(i,j).ROI_Info(1);
                    ROI_info_table(i,j).ROI_stdev  = ROI_info_table(i,j).ROI_Info(2);
                    ROI_info_table(i,j).ROI_pixels = ROI_info_table(i,j).ROI_Info(3);
                    ROI_info_table(i,j).ROI_min    = ROI_info_table(i,j).ROI_Info(4);
                    ROI_info_table(i,j).ROI_max    = ROI_info_table(i,j).ROI_Info(5);
                end
                
                val = get(ROI_info_table(i,j).ROI_Data(1), 'Userdata');
                ROI_info_table(i,j).ROI_Data(2,1:5) =  val(8:end);
            
                    
            end;
        end;
    end;
    if save_mat
        if ~isempty(fname)
            %eval(['save ''',pname fname , ''' ROI_info_table  ROI_string_table;' ])
            save([pname,fname], 'ROI_info_table', 'ROI_string_table');
        end;
    end
    if save_matlab_var
        assignin('base',roi_variable_name,ROI_info_table);
    end
end;

% if save_text
% 	if nargin < 2
%         fname = []; pname = [];
%         [fname,pname] =  uiputfile('*.txt', 'Save text file'); 
%     end;
%     fid = fopen([pname, [fname, '.txt']],'w');
%     for i = 1:size(ROI_string_table,1)
%         fprintf(fid, '%s\n', ROI_string_table(i,:)) ;
%     end;
%     fclose(fid);
% end;    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Load_ROI(filename);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to load an old ROI.mat file
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;
   
if isequal(filename,0)
    return
else
    ROI = load(filename);
    if ~isfield(ROI, 'ROI_info_table')
        return
    else
        % found file and it contains ROIs
        
        %begin by restoring blank state
        % call delete funtion with erase all ROIs
%         h_all_axes = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
%         fig = h_all_axes{1};    
%         fig2 = h_all_axes{2};
%         h_all_axes = h_all_axes{3};
%         
%         handles = guidata(fig2);
%         g = get(handles.Delete_ROI_popupmenu, 'Value');
%         set(handles.Delete_ROI_popupmenu, 'Value',4);
%         Delete_ROI;
%         set(handles.Delete_ROI_popupmenu, 'Value', g);
%         
%         
%         userdata = get(findobj(fig, 'Tag', 'figROITool'),'Userdata');
%         ROI_info_table = userdata{1};        
%         
%         
%         if isempty([ROI_info_table(:).ROI_Exists])  & size(ROI_info_table,2)==1
%             ROI_info_table = repmat(ROI_info_table, 1, length(h_all_axes(:)));
%         end;
        
        new_ROI_info_table = ROI.ROI_info_table;
        %save RRR1a new_ROI_info_table ROI_info_table;
        
%         if size(new_ROI_info_table,2)>= size(ROI_info_table,2)
%             % there are more images in the new table, get rid of extras
%             new_ROI_info_table = new_ROI_info_table(:,1:size(ROI_info_table,2));
%         else
%             % there are more images in current figure than in original figure,
%             % extend by creating an empty ROI
%             %b= size(new_ROI_info_table,2)+1:size(ROI_info_table,2)
%             new_ROI_info_table(size(new_ROI_info_table,1),size(ROI_info_table,2)).ROI_Exists = [];    
%         end;
        
        Refresh_ROIs(new_ROI_info_table);
        
    end;
    
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Close_Parent_Figure;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to make sure that if parent figure is closed, 
% the ROI info, ROI Tooland ROI Draw figures are closed too.
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

%set(findobj('Tag', 'RT_Info_figure'), 'Closerequestfcn', 'closereq');
try 
    delete(findobj('Tag','RT_Info_figure'));
end;

%set(findobj('Tag', 'ROI_figure'), 'Closerequestfcn', 'closereq');
try 
    delete(findobj('Tag','ROI_figure'));
end;
    
%set(findobj('Tag', 'ROI_Draw_figure'), 'Closerequestfcn', 'closereq');
try 
    delete(findobj('Tag','ROI_Draw_figure'));
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Menu_ROI_Tool;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

hNewMenu = gcbo;
checked=  umtoggle(hNewMenu);
hNewButton = get(hNewMenu, 'userdata');

if ~checked
    % turn off button
    %Deactivate Button;
    set(hNewMenu, 'Checked', 'off');
    set(hNewButton, 'State', 'off' );
else
    %Activate Button;
    set(hNewMenu, 'Checked', 'on');
    set(hNewButton, 'State', 'on' );
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Support Routines %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions called internally function and not as callbacks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function h_axes = Sort_Axes_handles(h_all_axes);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% receives a column vector of handles and 
% returns a matrix depending onthe position of 
% each image on the screen
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% assumes axes are in a grid pattern
% so sort them by position on the figure
for i = 1:length(h_all_axes);
    position(i,:) = get(h_all_axes(i),'Position');
end;

% calculate the different number of row values and the different number of column values to 
% set the matrix size
[hist_pos_y, bins_y] = hist(position(:,1));
[hist_pos_x, bins_x] = hist(position(:,2));
hy = sum(hist_pos_y>0);
hx = sum(hist_pos_x>0) ;
[hist_pos_y, bins_y] = hist(position(:,1), hy);
[hist_pos_x, bins_x] = hist(position(:,2), hx);

%hist_pos_x = fliplr(hist_pos_x);
h_axes = zeros(hx,hy);

sorted_positions = sortrows([position, h_all_axes], [2,1]); % sort x, then y
counter = 0;
for i =1:length(hist_pos_x)
    for j = 1:hist_pos_x(i)
        sorted_positions(j+counter,6) = hx - i + 1;
    end;
    counter = counter + hist_pos_x(i);  
end;

sorted_positions = sortrows(sorted_positions,[1,2]); % sort y, then x
counter = 0;
for i =1:length(hist_pos_y)
    for j = 1:hist_pos_y(i)
        sorted_positions(j+counter,7) = i;
    end;
    counter = counter + hist_pos_y(i);
end;

for i = 1:size(sorted_positions,1)
    h_axes(round(sorted_positions(i,6)),round(sorted_positions(i,7))) = sorted_positions(i,5);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Change_Object_Enable_State(handles, State, Paste_Flag)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

set(handles.Save_ROI_pushbutton, 'Enable', State);
set(handles.Copy_ROI_pushbutton, 'Enable', State);
set(handles.Edit_ROI_pushbutton, 'Enable', State);

set(handles.Delete_ROI_pushbutton, 'Enable', State);

set(handles.Save_MAT_checkbox, 'Enable', State);
set(handles.Save_matlab_var_checkbox, 'Enable', State);
set(handles.Delete_ROI_popupmenu, 'Enable', State);

set(handles.Link_ROI_togglebutton, 'Enable', State);

if Paste_Flag
    set(handles.Paste_ROI_pushbutton, 'Enable', State);
    set(handles.Paste_ROI_checkbox, 'Enable', State);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ifig = Create_ROI_Info_Figure(i_current_ROI)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create and initialize the ROI info figure
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

close_str = [ 'hNewButton = findobj(''Tag'', ''figROITool'');' ...
        ' if strcmp(get(hNewButton, ''Type''), ''uitoggletool''),'....
        ' set(hNewButton, ''State'', ''off'' );' ...
        ' else,  ' ...
        ' ROI_tool(''Deactivate_ROI_Tool'' ,hNewButton);'...
        ' set(hNewButton, ''Value'', 0);',...
        ' end;' ];

ifig = openfig('ROI_info_figure');
set(ifig, 'Name', 'ROI Information' ,...
    'Tag', 'RT_Info_figure' ,...
    'Resize' , 'On',...
    'CloseRequestfcn', close_str...
    );
handles = guihandles(ifig);
guidata(ifig, handles);

% add figure handle into the "storage" place
info_values = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
findobj(info_values{1},   'Tag', 'figROITool');
userdata = get(findobj(info_values{1}, 'Tag', 'figROITool'),'Userdata');
ROI_info_table = userdata{1};
handles_RT_Tool = guidata(info_values{2});

for i = 1:size(ROI_info_table,1)
    for j = 1:size(ROI_info_table,2)
        if ROI_info_table(i,j).ROI_Exists
            current_values = ROI_info_table(i,j).ROI_Info;
            % note that MATLAB automatically pads the strings with empty spaces
            String_info_table(i,j, :)= Convert_ROI_Info([j,i,current_values]);
        end;
    end
end;
        
    % now store the table in the listbox
set(handles.ROI_Info_listbox,'Userdata', String_info_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Highlight_Current_ROI(i_current_ROI);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_listbox = findobj('Tag','ROI_Info_listbox');
String_info_table = get(h_listbox, 'Userdata');
current_page = get(h_listbox, 'String');

if ~isempty(i_current_ROI)
    current_string = squeeze(String_info_table(i_current_ROI(1,1), i_current_ROI(1,2),:))';
    for i = 1:size(current_page, 1)
        if strcmp(current_string(1,1:end-1), current_page(i,:))
            set(h_listbox, 'Value', i);
        end
    end;
else
    % want to set current ROI to blank (due to deletion of current ROI)
    set(h_listbox, 'Value', []);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Info_string = Convert_ROI_Info(Info_numbers)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that converts the information from the ROI info table (numbers)
% into a string that can be ins14erted in a cell array for display in the
% list box.
% temp fixed spacings
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;
a = [     sprintf('%2d', Info_numbers(1))];
a = [a,   sprintf('%4d', Info_numbers(2))];

total_use_digits = 9;
if Info_numbers(6) < 10^5 
    after_decimal_precision_digits = 2;
    total_precision_digits =  5;
else
    after_decimal_precision_digits = 0;
    total_precision_digits = total_use_digits;
end;

a = [a,  FixLengthFormat(Info_numbers(3),total_use_digits, after_decimal_precision_digits)];
a = [a,  FixLengthFormat(Info_numbers(4),total_use_digits-1, after_decimal_precision_digits)];

a = [a,  sprintf('%6s',  num2str(Info_numbers(5), total_precision_digits) )];

a = [a,  FixLengthFormat(Info_numbers(6),total_use_digits-2, after_decimal_precision_digits)];
a = [a,  FixLengthFormat(Info_numbers(7),total_use_digits, after_decimal_precision_digits)];
Info_string = [a, 'x'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Update_ROI_Info_String(update_list)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inserts new data into the string array which is later "published"
% onto the listbox
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_listbox = findobj('Tag','ROI_Info_listbox');
String_info_table = get(h_listbox, 'Userdata');

info_values = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
userdata = get(findobj(info_values{1}, 'Tag', 'figROITool'),'Userdata');
ROI_info_table = userdata{1};

if nargin == 0
    update_list = [];
    % no update list sent, so update all values in the ROI_info_table
    for i = 1:size(ROI_info_table,1)
        for j = 1:size(ROI_info_table,2)
            update_list(size(update_list,1)+1,:) = [i,j];
        end;
    end;    
end

for i = 1:size(update_list,1)
    current_info = [update_list(i,[2,1]), ROI_info_table(update_list(i,1), update_list(i,2)).ROI_Info];
    % if ROI has been deleted, current info second hald will be empty will be empty; Fill with spaces
    if length(current_info)>2
        String_info_table(update_list(i,1), update_list(i,2),:) = Convert_ROI_Info(current_info);
    else
        String_info_table(update_list(i,1), update_list(i,2),:) =' ';
    end;
end;
set(h_listbox, 'Userdata', String_info_table);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function value = FixLengthFormat(num,totalChars, precision)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   By Patrick 'fishing boy' Helm
%	FixLengthFormat(num,totalChars, precision)
%		num is the number to format.	
%		totalChars is number of characters in conversion.
%				The number if right justified if the number
%				of characters is greater than the length of
%				num2str(num)
%		precision is number of digits after decimal
%
tNum = num2str(num, 16);
% find the decimal point index, if it exists
iDecimal = find(tNum == '.');
if (isempty(iDecimal)) 
    % set decimal to end of number if none 
    iDecimal = length(tNum);
    precision = 0;
else
    % add on zeroes until precision requirement is met
    while((length(tNum) - iDecimal) < 16)
        tNum = [tNum,'0'];
    end
end
% insure that even if function fails, blanks are returned
value = blanks(totalChars);
% copy character version onto output, 
% maintaining right justification
if ((iDecimal + precision) <= totalChars)
   startPos = totalChars - (precision + iDecimal) + 1;
   value(startPos:totalChars) = tNum(1:(iDecimal+precision));
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function h_all_axes = Find_All_Axes(varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to find and sort axes in a figure - or - 
% get axes handles if array of image handles is sent in
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

if strcmp(get(varargin{1}(1), 'Type'),'figure')
    % sent in nothing; determine axes handles and sort them into the correct matrix
    h_all_axes = Sort_Axes_handles(findobj(varargin{1}(1),'Type', 'Axes'));  
else
    % sent in the image handles, already sorted into matrix; now find parent axes for each one
    % but don't include them if the image is a Colorbar
    h_images = varargin{1};
    for i =1:size(h_images,1)
        for j = 1:size(h_images,2)
            if h_images(i,j)~= 0
                if ~strcmp( get(h_images(i,j), 'Tag'), 'TMW_COLORBAR')
                    h_all_axes(i,j) = get(h_images(i,j),'Parent');
                end;
            end;
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Refresh_ROIs(ROI_table);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Support function to redraw ROIs from the ROI info table
% assumes a start from the Delete_ROI, case 4 condition;
% Parallels Create_ROI code;
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

% h_all_axes = get(findobj('Tag', 'ROI_Title_text'), 'UserData');
% fig = h_all_axes{1};
% fig2 = h_all_axes{2};
% h_axes = h_all_axes{4};
% i_current_ROI = h_all_axes{5};
% h_all_axes = h_all_axes{3};
% 
% h_all_axes = h_all_axes';
% h_all_axes = h_all_axes(find(h_all_axes));
% 
% handles = guidata(fig2);
% 
% % get old ROI_info table
% Userdata = get(findobj(fig, 'Tag', 'figROITool'), 'UserData');
% 
% ROI_info_table = Init_ROI_info_table;
% 
% temp_h_all_axes = find(h_all_axes);
% % cycle through to create 0's in the Exist field. marks the ROI table as initialized
% [ROI_info_table(1:size(ROI_table,1),1:length(temp_h_all_axes)).ROI_Exists] = deal(0);
first_ROI_flag = 1;    
    

% find empty slots in ROI_info Table
Current_ROI_index = 1;

% colors: red, green, blue, yellow, magenta, cyan, white 
colororder = repmat('rbgymcw',1,4);

ROI_counter = 1;

if ~isfield(ROI_table, 'Other_coordinates')
    % the information necessary for all the other ROIs - IE user specified
	% ROIs that do not have positions available for the action objects
    ROI_table = Insert_Other_Coordinates(ROI_table);
end;
%save RRR2 ROI_table

fig = gcf;
h_axes_index = gca;

ind = 1;
update_list = [];
for i = 1:size(ROI_table,2)
    % for every ROI in the table
    for j = 1:size(ROI_table,1)


        if ROI_table(j,i).ROI_Exists

			
			%CALCULATION

			
            Current_ROI_index = j;
            x = ROI_table(j,i).ROI_x_coordinates;
            y = ROI_table(j,i).ROI_y_coordinates;
            
            hold on
            if(ind==1)
                plot(x, y, 'g');
            else
                plot(x, y, 'k');
            end
            hold off
            
            ind = ind + 1;
        end;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_table = Insert_Other_Coordinates(ROI_table);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Support function to insert "other coordinates" to elliptical ROIs that do not have them:
% need positions for the size, angle and number handles. 
% Places the results of the calculations into the ROI_Data field of each element
% in the ROI info table. Elliptical ROI must have even number of non-equal points. If ellipse is 
% closed, then there must be an odd number of points. Prefereably number of points divisible by 4.
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

debug = 0;
colororder = repmat('rbgymcw',1,4);

if isfield(ROI_table, 'Other_coordinates')
    old_Table = ROI_table;
    ROI_table = rmfield(ROI_table, 'Other_coordinates');
end;

for i = 1:size(ROI_table,2)
    % for every ROI in the table
    if debug==1
        ff = figure;
    end;        
    for j = 1:size(ROI_table,1)
        
        if ROI_table(j,i).ROI_Exists
            
            Current_ROI_index = j;
            x = ROI_table(j,i).ROI_x_coordinates;
            y = ROI_table(j,i).ROI_y_coordinates;
            
            % clear out the field if it does exist
            
            % the information necessary for all the other 
            % objects was not stored. Make an attempt to recover it.
            if debug==1
                plot(1,1);
                hold on;
                set(gca, 'ydir','reverse'); 
                set(gca, 'xlim' ,[ 0 256], 'ylim', [ 0 256]);
                axis('equal');                
            end;
            
            % un-close the ellipse
            if (x(1) == x(end))
                x = x(1:end-1);
                y = y(1:end-1);
            end;
            length_x = length(x);
            
            % fit the ellipse
            %x2 = [x(1)+ 0.1*x(1), x(2:end)];
            aa = fit_ellipse(x,y);
            
            f = aa(6); e = aa(5); d = aa(4);
            c = aa(3); b = aa(2); a = aa(1);
            
            % determine position of center of ROI 
            x_center = mean(x);
            y_center = mean(y); 
                        
            % determine major and minor axes of ellipse
            xdiff = x(1:length_x/2) - x(length_x/2+1:end);
            ydiff = y(1:length_x/2) - y(length_x/2+1:end);
            xydiff = sqrt(xdiff.^2 + ydiff.^2);
            
            % if the coefficient of variation of the distances between opposing points
            % is close to zero (implies circular ROI), then use the standard visual places 
            % for the ROIs: upper left size, upper right number, right angle 
            cv = abs(std(xydiff)/mean(xydiff));
            
            if  cv > 0.01 
                max_xydiff = max(xydiff);
                major_axis = find(xydiff==max_xydiff );
            else
                major_axis = find(x== min(x));
            end;
            
            % the minor axis is 1/4 of the ellipse around
            % assume the minor axis is pi/2 away from major axis. This is not exacly true for 
            % ROIs with num_points/4 
            minor_axis = major_axis + floor(length_x/4);
            
            % take first if there are many equal diameters
            % wrap the minor_axis to the length of points
            major_axis = major_axis(1);
            minor_axis = mod(minor_axis(1),length(x));
            
            %plot([x(major_axis), x(mod(major_axis + length_x/2, length_x))],...
            %    [y(major_axis), y(mod(major_axis + length_x/2, length_x))], 'g-');
            %plot([x(minor_axis), x(mod(minor_axis + length_x/2, length_x))],...
            %    [y(minor_axis), y(mod(minor_axis + length_x/2, length_x))], 'g:');
            
            % determine angle positions the cheap way
            %x_angle1 =  x(mod(major_axis + length_x/2,length_x))
            %y_angle1 =  y(mod(major_axis + length_x/2,length_x))
            %plot(x_angle1, y_angle1, 'kx');
            
            major_slope = ( y(mod(major_axis + length_x/2 - 1, length_x)+1) - y(major_axis) )/ ...
                ( x(mod(major_axis + length_x/2 -1, length_x) +1) - x(major_axis));

            %minor_slope = ( y(mod(minor_axis + length_x/2 - 1, length_x)+1) - y(minor_axis) )/ ...
            %    ( x(mod(minor_axis + length_x/2 -1, length_x) +1) - x(minor_axis))
            R = [ cos(pi/2), sin(pi/2); -sin(pi/2), cos(pi/2)];
            minor_slope = R * [ 1 major_slope]' ;
            minor_slope = minor_slope(2)/ minor_slope(1);
            if isnan(minor_slope), minor_slope = 0; end;
            
            
            xmin = min(x);
            ymin = min(y);
            xmax = max(x);
            ymax = max(y);
            
            % use slope calc intersections of line and ellipse
            s = major_slope;
            if ~isinf(s)
                % intercept of lines crossing ellipse
                % with exactly real roots - i.e. single point intersections 
                t1_major = 1/2/(4*a*c-b^2)*(-2*e*s*b-4*a*e+4*c*s*d+2*d*b+4*(e^2*s*b*a-e*s^2*b*c*d-e*s*b^2*d+a^2*e^2-a*e*d*b+c^2*s^2*d^2+c*s*d^2*b+a*c*e^2*s^2+a*c*d^2-4*a*c*b*s*f-4*a^2*c*f-4*a*c^2*s^2*f+b^3*s*f+b^2*a*f+b^2*c*s^2*f)^(1/2));
                t2_major = 1/2/(4*a*c-b^2)*(-2*e*s*b-4*a*e+4*c*s*d+2*d*b-4*(e^2*s*b*a-e*s^2*b*c*d-e*s*b^2*d+a^2*e^2-a*e*d*b+c^2*s^2*d^2+c*s*d^2*b+a*c*e^2*s^2+a*c*d^2-4*a*c*b*s*f-4*a^2*c*f-4*a*c^2*s^2*f+b^3*s*f+b^2*a*f+b^2*c*s^2*f)^(1/2));
                %    plot([ xmin xmax ], [ xmin*s + t1_major ,xmax*s + t1_major], 'b-')
                %    plot([ xmin xmax ], [ xmin*s + t2_major ,xmax*s + t2_major], 'b-') 
            else
                t1_major = NaN;
                t2_major = NaN;
            end;
            
            
            s = minor_slope;
            % intercept of lines crossing ellipse, for minor slope
            if ~isinf(s)
                t1_minor = 1/2/(4*a*c-b^2)*(-2*e*s*b-4*a*e+4*c*s*d+2*d*b+4*(e^2*s*b*a-e*s^2*b*c*d-e*s*b^2*d+a^2*e^2-a*e*d*b+c^2*s^2*d^2+c*s*d^2*b+a*c*e^2*s^2+a*c*d^2-4*a*c*b*s*f-4*a^2*c*f-4*a*c^2*s^2*f+b^3*s*f+b^2*a*f+b^2*c*s^2*f)^(1/2));
                t2_minor = 1/2/(4*a*c-b^2)*(-2*e*s*b-4*a*e+4*c*s*d+2*d*b-4*(e^2*s*b*a-e*s^2*b*c*d-e*s*b^2*d+a^2*e^2-a*e*d*b+c^2*s^2*d^2+c*s*d^2*b+a*c*e^2*s^2+a*c*d^2-4*a*c*b*s*f-4*a^2*c*f-4*a*c^2*s^2*f+b^3*s*f+b^2*a*f+b^2*c*s^2*f)^(1/2));
                %    plot([ xmin xmax ], [ xmin*s + t1_minor ,xmax*s + t1_minor], 'm-')
                %    plot([ xmin xmax ], [ xmin*s + t2_minor ,xmax*s + t2_minor], 'm-') 
            else
                t1_minor = NaN;
                t2_minor = NaN;
            end;
            
            
            if ~isnan(t1_minor) & ~isnan(t1_major)
                
                x_t1_t1_intersect = (t1_major - t1_minor)  / (minor_slope - major_slope);
                y_t1_t1_intersect = major_slope*x_t1_t1_intersect + t1_major;
                
                x_t1_t2_intersect = (t1_major - t2_minor)  / (minor_slope - major_slope);
                y_t1_t2_intersect = major_slope*x_t1_t2_intersect + t1_major;
                
                x_t2_t1_intersect = (t2_major - t1_minor)  / (minor_slope - major_slope);
                y_t2_t1_intersect = major_slope*x_t2_t1_intersect + t2_major;
                
                x_t2_t2_intersect = (t2_major - t2_minor)  / (minor_slope - major_slope);
                y_t2_t2_intersect = major_slope*x_t2_t2_intersect + t2_major;
                
            else
                % one of the slopes is infinte, therefore the other is zero
                
                x_t1_t1_intersect = xmin;
                y_t1_t1_intersect = ymin;
                
                x_t1_t2_intersect = xmin;
                y_t1_t2_intersect = ymax;
                
                x_t2_t1_intersect = xmax;
                y_t2_t1_intersect = ymin;
                
                x_t2_t2_intersect = xmax;
                y_t2_t2_intersect = ymax;
                
            end;
            x_corners = [...
                    x_t1_t1_intersect,...
                    x_t1_t2_intersect,...
                    x_t2_t1_intersect,...
                    x_t2_t2_intersect ...
            ];
            y_corners = [...
                    y_t1_t1_intersect,...
                    y_t1_t2_intersect,...
                    y_t2_t1_intersect,...
                    y_t2_t2_intersect ...
            ];
            
            %plot(x_corners(1), y_corners(1), 'ms')
            %plot(x_corners(2), y_corners(2), 'mo')
            %plot(x_corners(3), y_corners(3), 'md')
            %plot(x_corners(4), y_corners(4), 'm*')
            
            corners = [x_corners' , y_corners'];
            
            % determine the position of the angle marker by finding the point
            % on the ellipse, where the line with slope major_slope crosses.
            s = major_slope;
            t = y_center - x_center*s;
            if isinf(s) , 
                s = minor_slope; 
                t = y_center - x_center*s;
                xxx_plus = x_center;
                yyy_plus = min(y);
                
                x_angle_plus_180 = xxx_plus;
                y_angle_plus_180 = yyy_plus;
                
            else
                xxx_plus = [ 1/2/(b*s+a+c*s^2)*(-2*c*t*s-d-e*s-b*t+(4*c*t*s*d+d^2+2*d*e*s+2*d*b*t+e^2*s^2-2*e*s*b*t+b^2*t^2-4*b*s*f-4*a*f-4*a*c*t^2-4*a*e*t-4*c*s^2*f)^(1/2))]      ;
                yyy_plus = s*xxx_plus + t;
                
                x_angle_plus_180 = xxx_plus;
                y_angle_plus_180 = yyy_plus ; 
                
                
            end;
            xxx_minus = [ 1/2/(b*s+a+c*s^2)*(-2*c*t*s-d-e*s-b*t-(4*c*t*s*d+d^2+2*d*e*s+2*d*b*t+e^2*s^2-2*e*s*b*t+b^2*t^2-4*b*s*f-4*a*f-4*a*c*t^2-4*a*e*t-4*c*s^2*f)^(1/2))]      ;
            yyy_minus = s*xxx_minus + t;
            
            x_angle = xxx_minus;
            y_angle = yyy_minus;
            
            
            % determine the order of the points
            vec1 = corners - repmat([x_center, y_center],4,1);
            vec1(5,:) =  [x_angle , y_angle] - [x_center, y_center];
            
            [theta, vec_size] = cart2pol(vec1(:,1), vec1(:,2));
            theta = theta*180/pi;
            
            sorted = sortrows([theta, vec1], 1);
            angle_position = find(sorted(:,1)'== theta(5));
            
            x_square = sorted( mod(angle_position - 3,5) +1 ,2) + x_center;
            y_square = sorted( mod(angle_position - 3,5) +1,3) + y_center;
            
            x_number = sorted( mod(angle_position - 2,5) +1,2) + x_center;
            y_number = sorted( mod(angle_position - 2,5) +1,3) + y_center;
            
            size_x = norm( [x_angle - x_center, y_angle - y_center] );
            size_y = norm( [x_angle - x_number, y_angle - y_number] );
            
            
            % calculate angle
            v1 = [1 0];  % basis
            v2 = [x_angle, y_angle] - [x_center, y_center] ;
            d = cross([v1 0],[v2 0]);
            % change this calculation to an atan2 calculation
            alpha = acos(  dot(v1,v2)   /(norm(v1) * norm(v2)) ) *sign(d(3));
            
            
            % close the polygon
			if (x(end)~=x(1)) & (y(end)~=y(1))
				x = [x , x(1)];
				y = [y , y(1)];
			end;
            
            if debug==1
                % now plot circle with basic points, transformed by skew, rotation, and translation
                h_circle = plot(x,y,[colororder(Current_ROI_index),'-']);
                % 'ButtonDownFcn', 'ROI_tool(''Change_Current_ROI'')');
            
                h_center = plot(x_center, y_center , [colororder(Current_ROI_index),'+']);
                % 'ButtonDownFcn', 'ROI_tool(''ROI_Pos_Adjust_Entry'',1)'); 
                
                h_size = plot(x_square, y_square,...
                    [colororder(Current_ROI_index),'s'] );
                %  'ButtonDownFcn', 'ROI_tool(''ROI_Size_Adjust_Entry'')');    
                
                h_angle = plot(x_angle, y_angle,...
                    [colororder(Current_ROI_index),'o']);
                % 'ButtonDownFcn', 'ROI_tool(''ROI_Angle_Adjust_Entry'')');
                
                h_number = text(x_number, y_number, num2str(Current_ROI_index),'color', ...
                    [colororder(Current_ROI_index)], 'HorizontalAlignment', 'center' );
                % 'ButtonDownFcn', 'ROI_tool(''ROI_Pos_Adjust_Entry'',2)'); 
                
            end;
            confirmation = [ROI_table(j,i).ROI_Data(2,:); ...
                    x_center, y_center, size_x, size_y, alpha];
    
            ROI_table(j,i).ROI_Data(2,:) =  [x_center, y_center, size_x, size_y, alpha]
            ROI_table(j,i).Other_coordinates= [x_center, y_center, x_square, y_square, ...
                x_angle, y_angle, x_number, y_number];
            
            
        end;
       
    end;
    
end;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function a = fit_ellipse(x, y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Support function to fit an ellipse from given x,y data points.
% Used to recalculate corner points for elliptical ROIs loaded. 
% a=> [ a1*x^2 + a2*x*y a3*y^2 + a4*x + a5*y + a6 = 0]
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

if size(x,1) == 1
    x = x';
    y = y';
end;
D1 = [x.^2, x.*y, y.^2];                  % quadratic part of the design matrix
D2 = [x, y, ones(size(x))];               % linear part of the design matrix
S1 = D1'*D1;                              % quadratic part of the scatter matrix
S2 = D1'*D2;                              % combined part of the scatter matrix
S3 = D2'*D2;                              % linear part of the scatter matrix
T = -inv(S3)*S2';                         % for getting a2 from a1
M = S1 + S2*T;                            % reduced scatter matrix
M = [M(3,:)./2; - M(2,:); M(1,:)./2;];    % premultiply by inv(C1)
[e_vec, e_val] = eig(M);                    % solve eigensystem
cond = 4*e_vec(1,:).*e_vec(3,:) - e_vec(2,:).^2; % evaluate a'Ca
a1 = e_vec(:,find(cond > 0));              % eigenvector for min. pos. eigenvalue
a = [a1; T*a1];                           % ellipse coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Close_Old_Figure(Figure_Name, Figure_Handle);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

if ~isempty(Figure_Name)
	Figure_Handle = findobj('Tag', Figure_Name);
end;
if ~isempty(Figure_Handle)
	set(Figure_Handle,'CloseRequestFcn', 'closereq');
	try close(Figure_Handle);
	catch delete(Figure_Handle);
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [h_objects, object_enable_states] = Change_Figure_States(State, figure_handles,  h_objects, object_enable_states)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to set the enable states of the uicontrol objects of various
% figures. Figure handles sent in are used to find all the objects, 
% if the object handles are not sent in initially. 
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

if nargin <3
	h_objects = findobj(figure_handles, 'type', 'uicontrol');
	object_enable_states = {};
else
	if length(h_objects) ~= length(object_enable_states)
		disp(['Can''t Enable Objects!']);;
	end;
end;

if strcmp(State, 'Disable')
	object_enable_states = get(h_objects, 'Enable');
	set(h_objects, 'Enable', 'off');
elseif strcmp(State, 'Enable')
	if isempty(object_enable_states)
		set(h_objects, 'Enable', 'On');
	else
		for i=1:length(h_objects)
			set(h_objects(i), 'Enable', object_enable_states{i});
		end;
	end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [xs, ys] = Spline_ROI(x, y);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	interpolate x & y points into a spline curve
%   return xs and ys, upsampled by a factor f
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

f  = 1/Sample_Rate;      % Upsample curve by a factor of 10 
t  = 1:length(x);
ts = 1: f : length(x);
xs = spline(t, x, ts);
ys = spline(t, y, ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function F = Sample_Rate;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ROI_info_table = Init_ROI_info_table;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function to create an empty ROI_info_table with the 
%  correct fields initialized to empty
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

ROI_info_table = struct( ...
	'ROI_Data', [], ...
	'ROI_Exists', [], ...
	'ROI_Info', [], ...
	'ROI_x_original', [],...
	'ROI_y_original', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function  handle_values = Make_ROI_Elements(xs,ys,roi_color,roi_number, center_x, center_y, size_x, size_y, angle_x, angle_y, number_x, number_y, alpha0);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Funciton to create sub-elements of an ROI, including the actual ROI, the
%  resize square, the angle circle, the central plus and the ROI number
%global DB; if DB disp(['ROI_Tool: ', Get_Current_Function]); end;

h_circle = plot(xs,ys,[roi_color,'-'],...
	'ButtonDownFcn', 'ROI_tool(''Change_Current_ROI'')');

h_center = plot(center_x, center_y , [roi_color,'+'], ...
	'ButtonDownFcn', 'ROI_tool(''ROI_Pos_Adjust_Entry'',1)'); 

h_size = plot(size_x , size_y, [roi_color,'s'],...
	'ButtonDownFcn', 'ROI_tool(''ROI_Size_Adjust_Entry'')');

h_angle = plot(angle_x, angle_y, [roi_color,'o'],...
	'ButtonDownFcn', 'ROI_tool(''ROI_Angle_Adjust_Entry'')');
if nargin == 13
	setappdata(h_angle, 'alpha0', alpha0);	
end;

h_number = text(number_x, number_y, num2str(roi_number),...
	'color', roi_color, ...
	'HorizontalAlignment', 'center' , ...
	'ButtonDownFcn', 'ROI_tool(''ROI_Pos_Adjust_Entry'',2)'); 

handle_values = [h_circle, h_center, h_size, h_angle, h_number];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function  func_name = Get_Current_Function;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug function - returns current function name
x = dbstack;
x = x(2).name;
func_name = x(findstr('(', x)+1:findstr(')', x)-1);
