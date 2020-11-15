function varargout = CQmark(varargin)
% CQMARK M-file for CQmark.fig
%      CQMARK, by itself, creates a new CQMARK or raises the existing
%      singleton*.
%
%      H = CQMARK returns the handle to a new CQMARK or the handle to
%      the existing singleton*.
%
%      CQMARK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CQMARK.M with the given input arguments.
%
%      CQMARK('Property','Value',...) creates a new CQMARK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CQmark_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CQmark_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help CQmark

% Last Modified by GUIDE v2.5 21-Dec-2006 00:44:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CQmark_OpeningFcn, ...
                   'gui_OutputFcn',  @CQmark_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before CQmark is made visible.
function CQmark_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CQmark (see VARARGIN)

% Choose default command line output for CQmark
handles.output = hObject;

handles.path = 1;  % Number of traced path to be shown

handles.st = 1;
handles.en = 1;
%set(handles.axes1,'nextplot','add','ydir','reverse','layer','top');%, 'XLim', '[1 handles.sizec]', 'YLim','[1 handles.sizer]');
handles.marker=[];
%variable for maker handles and locations
handles.cp=[];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CQmark wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CQmark_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1}=handles.marker;


% Get default command line output from handles structure
%varargout{1} = handles.output;

% --- Executes on slider movement.
function slide_frame_Callback(hObject, eventdata, handles)
% hObject    handle to slide_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.oldn=handles.n;
sv = get(handles.slide_frame,'Value');
handles.n=round(sv*max(size(handles.mov)));
if handles.n ==0, handles.n=1;end
set(handles.frame_no,'String',handles.n);
guidata(hObject, handles);

imageupdate(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slide_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slide_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function frame_no_Callback(hObject, eventdata, handles)
% hObject    handle to frame_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_no as text
%        str2double(get(hObject,'String')) returns contents of frame_no as a double

handles.oldn=handles.n;
handles.n = round(str2double(get(handles.frame_no,'String')));

if handles.n>max(size(handles.mov))
    handles.n = max(size(handles.mov));
    set(handles.frame_no,'String',handles.n);
end
set(handles.slide_frame,'Value',handles.n/max(size(handles.mov)));
guidata(hObject, handles);
imageupdate(hObject,handles);

% --- Executes during object creation, after setting all properties.
function frame_no_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Function imageupdate
function imageupdate(hObject,handles)


[handles.I,Map] = frame2im(handles.mov(1,handles.n));
  if get(handles.pushbutton5,'value')
      handles.I=imcomplement(handles.I);  %% invert
     
  end
    delete(handles.marker(handles.oldn).markerHandles);
    handles.marker(handles.oldn).markerHandles=[];
    set(handles.imageHandle,'cData',handles.I);
   
    if length(handles.marker(handles.n).xlocation)
        for indx=1:length(handles.marker(handles.n).xlocation)
            handles.marker(handles.n).markerHandles(indx,:)=plot(get(handles.imageHandle,'Parent'),...
                handles.marker(handles.n).xlocation(indx),handles.marker(handles.n).ylocation(indx),'ro',...
                handles.marker(handles.n).xlocation(indx),handles.marker(handles.n).ylocation(indx),'r+','MarkerSize',15);
        end
        %handles.marker(handles.n).markerHandles(indx)
    end

guidata(hObject, handles);    


%require callback property
%CQmark('marker_Callback',gcbo,[],guidata(gcbo))
% --- Executes on button press in marker.
function marker_Callback(hObject, eventdata, handles)
% hObject    handle to marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in del_marker.
function del_marker_Callback(hObject, eventdata, handles)
% hObject    handle to del_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.marker(handles.n).xlocation=handles.marker(handles.n).xlocation(1:end-1);
handles.marker(handles.n).ylocation=handles.marker(handles.n).ylocation(1:end-1);
delete(handles.marker(handles.n).markerHandles(end,:));
handles.marker(handles.n).markerHandles=handles.marker(handles.n).markerHandles(1:end-1,:);
length(handles.marker(handles.n).markerHandles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Open_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pack
handles.I = [];
handles.mov = [];

[filename, pathname] = uigetfile('*.mat', 'Pick a mat-file');

if isequal(filename,0)
    disp('User selected Cancel')
else
    file = fullfile(pathname, filename);
    tmpimvar=load(filename);
    handles.mov=tmpimvar.imgt;
    
    disp(['User selected ', file])
    name = ['CQmark - ' file];
    set(CQmark,'Name',name)
    
    handles.nof = max(size(handles.mov));
    handles.n = handles.nof;
    handles.marker=struct('markerHandles',[],'xlocation',[],'ylocation',[]);
    for indx=1:handles.nof
        handles.marker(indx).xlocation=[];
        handles.marker(indx).ylocation=[];
        handles.marker(indx).markerHandles=[];
    end
    set(handles.slide_frame,'Value',1);
    set(handles.frame_no,'String',handles.n);
    [handles.I,Map] = frame2im(handles.mov(1,handles.n));
    set(gcf,'units','pixels')
    set(gca,'units','pixels')
    [handles.sizer,handles.sizec,handles.sizez]=size(handles.I);
    set(gca,'pos',[50,150,handles.sizec,handles.sizer]);
    handles.imageHandle=imagesc(handles.I,'Parent',handles.axes1);
    
    if get(handles.pushbutton5,'value')
        handles.I=imcomplement(handles.I);  %invert
    end
    
    set(handles.figure1,'colormap',Map);
    set(handles.imageHandle,'ButtonDownFcn','CQmark(''image_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
    set(handles.slide_frame,'Enable','on');
    set(handles.frame_no,'Enable','on');
    set(handles.slide_frame,'SliderStep',[1/handles.nof;.1]);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function load_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to load_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pack
disp(' Loading data from Workspace ... ');
handles.mov = [];

handles.mov = evalin('base','imgt');
handles.nof = max(size(handles.mov));
handles.n = handles.nof;
handles.marker=struct('markerHandles',[],'xlocation',[],'ylocation',[]);
    for indx=1:handles.nof
        handles.marker(indx).xlocation=[];
        handles.marker(indx).ylocation=[];
        handles.marker(indx).markerHandles=[];
    end
set(handles.slide_frame,'Value',1);
set(handles.frame_no,'String',handles.n);
[handles.I,Map] = frame2im(handles.mov(1,handles.n));
set(gcf,'units','pixels')
set(gca,'units','pixels')
[handles.sizer,handles.sizec,handles.sizez]=size(handles.I);
set(gca,'pos',[50,150,handles.sizec,handles.sizer]);
handles.imageHandle=imagesc(handles.I,'Parent',handles.axes1); 

if get(handles.pushbutton5,'value')
    handles.I=imcomplement(handles.I);  %invert
end

set(handles.figure1,'colormap',Map);
set(handles.imageHandle,'ButtonDownFcn','CQmark(''image_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
set(handles.slide_frame,'Enable','on');
set(handles.frame_no,'Enable','on');
set(handles.slide_frame,'SliderStep',[1/handles.nof;.1]);
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get(get(handles.imageHandle,'Parent'))
if get(handles.markerButton,'value')
    cp=get(get(handles.imageHandle,'Parent'),'CurrentPoint');
    handles.marker(handles.n).xlocation(end+1)=cp(1);
    handles.marker(handles.n).ylocation(end+1)=cp(3);
    set(get(handles.imageHandle,'Parent'),'NextPlot','add');
    handles.marker(handles.n).markerHandles(end+1,:)=plot(get(handles.imageHandle,'Parent'),handles.marker(handles.n).xlocation(end),handles.marker(handles.n).ylocation(end),'ro',...
                                                     handles.marker(handles.n).xlocation(end),handles.marker(handles.n).ylocation(end),'r+','MarkerSize',15);    
end
guidata(hObject, handles);





% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[filename, pathname] = uigetfile('*.mat', 'Pick a mat-file');
if isequal(filename,0)
    disp('User selected Cancel')
else
    file = fullfile(pathname, filename);
    tmpvar=load(filename);
    handles.marker=tmpvar.output;
    set(handles.axes1,'nextplot','add','ydir','reverse','layer','top');
   
    if length(handles.marker(handles.n).xlocation)
        for indx=1:length(handles.marker(handles.n).xlocation)
            handles.marker(handles.n).markerHandles(indx,:)=plot(handles.axes1,handles.marker(handles.n).xlocation(indx),handles.marker(handles.n).ylocation(indx),'ro',...
                handles.marker(handles.n).xlocation(indx),handles.marker(handles.n).ylocation(indx),'r+','MarkerSize',15);
        end
    %handles.marker(handles.n).markerHandles(indx)
    end
   
 guidata(hObject, handles); 
end
   


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.oldn=handles.n;

if get(handles.pushbutton5,'value')
    handles.I=imcomplement(handles.imageHandle);

end
 guidata(hObject, handles); 
imageupdate(hObject,handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton6
  if get(handles.pushbutton6,'value')
      
     set(get(handles.imageHandle,'Parent'),'xgrid','on','ygrid','on');
  else
      set(get(handles.imageHandle,'Parent'),'xgrid','off','ygrid','off');
  end



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CQgen



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uiputfile('*.mat');
for indx=1:handles.nof
    output(indx).xlocation=handles.marker(indx).xlocation;
    output(indx).ylocation=handles.marker(indx).ylocation;
end

save([path,file],'output');



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nom=str2double(get(handles.edit2,'String'));
stf=str2double(get(handles.edit3,'String'));
nof=str2double(get(handles.edit4,'String'));

if (isnan(nom) || isnan(stf) || isnan(nof))
    
    error('Oops! Enter the NaN fields')
end


for indx=stf:nof
    X(indx-stf+1,:)=handles.marker(indx).xlocation;
    Y(indx-stf+1,:)=handles.marker(indx).ylocation;
end

nop=max(size(X(nof-stf+1,:)));

verteb=(nop)/nom;

if (round(verteb) / verteb)==1

    mxp.xlocation=reshape(X',nom,verteb,nof-stf+1);
    mxp.ylocation=reshape(Y',nom,verteb,nof-stf+1);
else
    error('Oops! The number of markers in each frame is different!!')
end

[file,path]=uiputfile('*.mat');

save([path,file],'mxp');



