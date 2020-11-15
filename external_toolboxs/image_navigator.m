function varargout = image_navigator(varargin)
%IMAGE_NAVIGATOR M-file for image_navigator.fig
%      IMAGE_NAVIGATOR, by itself, creates a new IMAGE_NAVIGATOR or raises the existing
%      singleton*.
%
%      H = IMAGE_NAVIGATOR returns the handle to a new IMAGE_NAVIGATOR or the handle to
%      the existing singleton*.
%
%      IMAGE_NAVIGATOR('Property','Value',...) creates a new IMAGE_NAVIGATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to image_navigator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IMAGE_NAVIGATOR('CALLBACK') and IMAGE_NAVIGATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IMAGE_NAVIGATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Created by Dr. Angelo Zizzari - KeyRes Technologies
%
% E-mail: angelo.zizzari@keyres-technologies.com
%
% http://www.keyres-technologies.com/
%

% Edit the above text to modify the response to help image_navigator

% Last Modified by GUIDE v2.5 09-Jan-2006 14:55:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image_navigator_OpeningFcn, ...
                   'gui_OutputFcn',  @image_navigator_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before image_navigator is made visible.
function image_navigator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for image_navigator
handles.output = hObject;

a7 = imread('icon7e.bmp');
a8 = imread('icon8e.bmp');
a9 = imread('icon9e.bmp');

set(handles.pushbutton8,'CData',a8);
set(handles.pushbutton9,'CData',a9);
set(handles.pushbutton1,'CData',a7);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes image_navigator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = image_navigator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

w=handles.w_res;

    curr_lim1 = get(handles.axes1,['x','lim']);
    curr_lim2 = get(handles.axes1,['y','lim']);
    
    range_data1 = abs(diff(curr_lim1));
    range_data2 = abs(diff(curr_lim2));
    
    axes(handles.axes1);

    val1_old=handles.val1;

    val1=get(handles.slider1,'Value');

    val1_max=get(handles.slider1,'Max');

    handles.val1=val1;

    delta_pixel1=(val1-val1_old);

    new_lim1 = curr_lim1 + delta_pixel1;

    set(handles.axes1,['x','lim'],new_lim1); 
     
        axes(handles.axes2);
    
        cla

        colormap(gray);     

        imagesc(w); axis image; grid off;
                
        rectangle('Position', [new_lim1(1) curr_lim2(1) new_lim1(2)-new_lim1(1)  range_data2],...
            'EdgeColor','y');
        
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

w=handles.w_res;

    curr_lim1 = get(handles.axes1,['x','lim']);
    curr_lim2 = get(handles.axes1,['y','lim']);
    
    range_data1 = abs(diff(curr_lim1));
    range_data2 = abs(diff(curr_lim2));
    
    axes(handles.axes1);

    val2_old=handles.val2;

    val2=get(handles.slider2,'Value');

    val2_max=get(handles.slider2,'Max');

    handles.val2=val2;

    delta_pixel2=(val2-val2_old);

    new_lim2 = curr_lim2 - delta_pixel2;
        
    set(handles.axes1,['y','lim'],new_lim2);
    
        axes(handles.axes2);
    
        cla

        colormap(gray);     

        imagesc(w); axis image; grid off;
                
        rectangle('Position', [curr_lim1(1) new_lim2(1) range_data1 new_lim2(2)-new_lim2(1)],...
            'EdgeColor','y');
        
    
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% LOAD MAIN IMAGE

err=0;

[FileName,PathName] = uigetfile({'*.tif;*.TIF','TIF Files (*.tif, *.TIF)'; 
    '*.jpg;*.JPG','JPG Files (*.jpg, *JPG)';...
    '*.gif;*.GIF','GIF Files (*.gif, *.GIF)';...
    '*.bmp;*.BMP','BMP Files (*.bmp, *.BMP)';...   
    '*.*','All Files(*.*)'},'Select the input file');

if isequal(FileName,0) || isequal(PathName,0)
       
        warndlg('Error: user pressed cancel. Please, select an image file.','Load Error');
        
else
        [pathstr,name,ext,versn] = fileparts(FileName);
        
        if isequal(ext,'.tif') || isequal(ext,'.jpg') || isequal(ext,'.JPG')...
                || isequal(ext,'.TIF') || isequal(ext,'.BMP') || isequal(ext,'.GIF') ...
                || isequal(ext,'.bmp') || isequal(ext,'.gif'),
            
                set(handles.uipanel1,'Title',strcat('File : ', FileName));

                cd (PathName)

                handles.FileName=FileName;
                           
                w=imread(FileName);

                w_height = size(w,1);   
                w_width = size(w,2);
                
            % View Image      
                                
                figure(handles.figure1); axes(handles.axes1);
                
                cla

                colormap(gray);     

                imagesc(w); axis image; grid off;
            
                handles.w=w;
                handles.w_res=w;
                  
                set(handles.frame1,'Visible','off'); 
                
                figure(handles.figure1); axes(handles.axes2);
                
                cla

                colormap(gray);     

                imagesc(w); axis image; grid off;
                
                set(handles.frame2,'Visible','off'); 
                        
                set(handles.slider1,'Enable','off');
                set(handles.slider2,'Enable','off');

                curr_lim1=[0.5 w_width+0.5];
                curr_lim2=[0.5 w_height+0.5];
                
                set(handles.axes1,['x','lim'],curr_lim1);
                set(handles.axes1,['y','lim'],curr_lim2);

                range_data1 = abs(diff(curr_lim1));
                range_data2 = abs(diff(curr_lim2));  
    
        
            guidata(hObject,handles);
            
        else

                warndlg('Error: invalid file format. Please select an image file.','Load Error');
        
        end
            
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);

zoom(2)

w=handles.w_res;

    m = size(w,1);   
    n = size(w,2);

    curr_lim1 = get(handles.axes1,['x','lim']);
    curr_lim2 = get(handles.axes1,['y','lim']);
    
    range_data1 = abs(diff(curr_lim1));
    range_data2 = abs(diff(curr_lim2));  

    set(handles.slider1,'Max',n-range_data1+0.5);
    set(handles.slider2,'Max',m-range_data2+0.5);   

    set(handles.slider1,'Value',curr_lim1(1));
    handles.val1=curr_lim1(1);
    set(handles.slider2,'Value',(m-range_data2)+1-curr_lim2(1));
    handles.val2=(m-range_data2)+1-curr_lim2(1);
   
    set(handles.slider1,'Min',0.5);
    set(handles.slider2,'Min',0.5);
    
    set(handles.slider1,'Enable','on');
    set(handles.slider2,'Enable','on');
    
        axes(handles.axes2);
    
        cla

        colormap(gray);     

        imagesc(w); axis image; grid off;
                
        rectangle('Position', [curr_lim1(1) curr_lim2(1) range_data1  range_data2],...
            'EdgeColor','y');
    
    
guidata(hObject,handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


axes(handles.axes1);

    curr_lim1 = get(handles.axes1,['x','lim']);
    curr_lim2 = get(handles.axes1,['y','lim']);

    range_data1_old = abs(diff(curr_lim1));
    range_data2_old = abs(diff(curr_lim2));

    w=handles.w_res;

    m = size(w,1);   
    n = size(w,2);

    px=curr_lim1(1)-range_data1_old/2;
    py=curr_lim2(1)-range_data2_old/2;

    range_data1 = 2*range_data1_old;
    range_data2 = 2*range_data2_old;

    curr_lim1 = [px px+range_data1];
    curr_lim2 = [py py+range_data2];  

    set(handles.axes1,['x','lim'],curr_lim1);
    set(handles.axes1,['y','lim'],curr_lim2);


        if ((n-range_data1)>0) && ((m-range_data2)>0),
        
            set(handles.slider1,'Min',0.5);
            set(handles.slider2,'Min',0.5);

            setval1=curr_lim1(1);
            setval2=curr_lim2(1);           
            
            if curr_lim1(1)<1, 
                setval1=0.5; 
                range_data1=2*range_data1_old;
                curr_lim1=[0.5 range_data1+0.5];
                set(handles.axes1,['x','lim'],curr_lim1); 
            end
            
            if curr_lim1(1)>(n-range_data1), 
                range_data1=2*range_data1_old;
                setval1=(n-range_data1+0.5); 
                curr_lim1=[(n-range_data1)+0.5 n+0.5];
                set(handles.axes1,['x','lim'],curr_lim1);
            end
            
            if curr_lim2(1)<1, 
                setval2=0.5; 
                range_data2=2*range_data2_old;
                curr_lim2=[0.5 range_data2+0.5];
                set(handles.axes1,['y','lim'],curr_lim2);
            end
            
            if curr_lim2(1)>(m-range_data2), 
                range_data2=2*range_data2_old;
                setval2=(m-range_data2+0.5); 
                curr_lim2=[(m-range_data2)+0.5 m+0.5];
                set(handles.axes1,['y','lim'],curr_lim2);
            end        
            
            if ((n-range_data1)>0) && ((m-range_data2)>0),

                set(handles.slider1,'Value',setval1);
                handles.val1=setval1;
                set(handles.slider2,'Value',m-range_data2+1-setval2);
                handles.val2=m-range_data2+1-setval2;
                
                set(handles.slider1,'Max',n-range_data1+0.5);
                set(handles.slider2,'Max',m-range_data2+0.5);
                
                set(handles.slider1,'Enable','on');
                set(handles.slider2,'Enable','on');
                
            else
                
                set(handles.slider1,'Enable','off');
                set(handles.slider2,'Enable','off');
                
            end
            
    
        else

            curr_lim1=[0.5 n+0.5]; curr_lim2=[0.5 m+0.5];
            range_data1=n; range_data2=m;
            
            set(handles.axes1,['x','lim'],curr_lim1); 
            set(handles.axes1,['y','lim'],curr_lim2);
            
            set(handles.slider1,'Enable','off');
            set(handles.slider2,'Enable','off');
        
        end
            
    axes(handles.axes2);

    cla

    colormap(gray);     

    imagesc(w); axis image; grid off;

    rectangle('Position', [curr_lim1(1) curr_lim2(1) range_data1  range_data2],...
        'EdgeColor','y');
        
guidata(hObject,handles);
            