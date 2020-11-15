function light_direction_Callback(hObject,eventdata,handles)
% hObject    handle to light_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns light_direction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from light_direction


global V3D_HANDLES

options={'right','left','headlight'};

index=get(handles.light_direction,'Value');

V3D_HANDLES.light_direction=options{index};

camlight(V3D_HANDLES.camlight_handle,V3D_HANDLES.light_direction)
