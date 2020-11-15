close all;
clear all;


% Important! Use the command below to change MatLab working directory to
% your current directory:

cd C:\WSQ_Matlab;

% Please note that MatLab sometimes "forgets" and resets working directory
% to Windows user "temp" directory.
% As for example on Windows XP, MatLab resets the working directory to:
% "C:\Documents and Settings\UserName\Local Settings\Temp\"
% Place a copy of files "WSQ_library.dll", "WSQ_library_read_image.dll" and
% "WSQ_library_write_image.dll" to that directory in order to ensure
% that MatLab always finds them.




% As of MATLAB 6.5.1 (R13SP1) and later versions, it is possible to access functions
% defined in Windows standard dynamic link libraries (.dll) through the MATLAB command
% line using the LOADLIBRARY function. This feature is available for Windows platforms
% only for releases R13SP1 through R14SP1. As of release R14SP2, LOADLIBRARY is
% supported on both Windows and Linux.
% The difficulty arises if you want to use HBITMAP data type in function call to the DLL,
% because MatLab does not natively support HBITMAP data type. 
% The solution is to use wrapper DLL which wraps HBITMAP data type into native MatLab
% uint8 array for image storage using MatLab MEX interface. 
% Such a wrapper interface is provided in DLL files "WSQ_library_read_image.dll"
% and "WSQ_library_write_image.dll". 


 
% Unload WSQ Image Library if it is already loaded
 if libisloaded('WSQ_lib')
 unloadlibrary WSQ_lib;
 end

 % Load WSQ Image Library
 loadlibrary 'WSQ_library' 'WSQ_library.h' 'alias' 'WSQ_lib'; 


% To show functions available in the WSQ Image Library use the command below:
% libfunctionsview('WSQ_lib'); 
 
 

filename = 'sample_image.wsq';

picture_data = WSQ_library_read_image(filename);

figure;
imshow(picture_data);



% Make sure that image data array is in uint8 format.
% This step is not needed if data is already in uint8 format
picture_data = uint8(picture_data); 




% Function SetShowFilePropertiesDialog enables/disables invocation of
% graphic file properties dialog window.
% "file_properties_dialog" denotes integer with possible values:
%        1 - Show file properties dialog
%        0 - Do not show file properties dialog
% Function "void SetShowFilePropertiesDialog(int file_properties_dialog)" 
% has effect only when saving WSQ and TIFF file types.
%
% Make sure that function SetShowFilePropertiesDialog argument is in int32 format.
% MatLab by default stores number in real(double) format,
% thus the explicit conversion to int32 is necessary.
calllib('WSQ_lib', 'SetShowFilePropertiesDialog', int32(1));



filename_write = 'outputimage.bmp';


% Variable "filetype" denotes integer with possible values:
% 1 - WSQ  FBI's Wavelet Scalar Quantization  
% 2 - BMP  Windows Bitmap Graphics  
% 3 - TIFF  Tagged Information File Format (no LZW compression)  
% 4 - PNG  Portable Network Graphics  
% 5 - JPEG  Joint Photographic Experts Group  
% 6 - RGB  Silicon Graphics International (uncompressed)  
% 7 - TGA  Truevision Targa Graphic  
file_type = 2;


% Make sure that variable "file_type" is in int32 format.
% MatLab by default stores number in real(double) format,
% thus the explicit conversion to int32 is necessary
WSQ_library_write_image(picture_data, filename_write, int32(file_type));



% To register "WSQ Image Library" on your computer use the function below:
%
% calllib('WSQ_lib', 'RegisterWSQ' );

% Please note that MatLab sometimes "forgets" and resets working directory
% to Windows user "temp" directory.
% As for example on Windows XP, MatLab resets the working directory to:
% "C:\Documents and Settings\UserName\Local Settings\Temp\"
% Place a copy of the registration file "wsq_license.key" and files
% "WSQ_library.dll", "WSQ_library_read_image.dll" and
% "WSQ_library_write_image.dll" to that directory in order to ensure
% that MatLab always finds them


 

% Unload WSQ Image Library
 unloadlibrary WSQ_lib;
