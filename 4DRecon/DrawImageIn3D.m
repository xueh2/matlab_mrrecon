function [ fighandle, surfh ] = DrawImageIn3D(data, header, fighandle, transparency )
%============================
% [ fighandle, surfh ] = DrawImageIn3D(data, header, fighandle, transparency)
% This is a function for drawing slices in 3D space based on the dicom coordinates
% information.
% Input:
%   data : [column row SLC], [2nd dimension, 1st dimension, 3rd dimension]
%   header : the ftk header format, including positionPatient (tx ty tz) and orientationPatient (row vector; col vector; slc vector)
%   figh = figure handle
%   transparency = alpha for the image slice. 
%                   Use values < 1 to get partial transparency when
%                   plotting multiple slices in one figure
%======================================================================

% set the current point as the center of 3D volume
currentPoint3D = [header.sizeX/2 header.sizeY/2 header.sizeZ/2];                
currentPoint3DInWorld = image2WorldMrFtk(header, currentPoint3D);

xlim3D(1) = header.spacingX*-0.5;
xlim3D(2) = header.spacingX*(header.sizeX-0.5);

ylim3D(1) = header.spacingY*-0.5;
ylim3D(2) = header.spacingY*(header.sizeY-0.5);

zlim3D(1) = header.spacingZ*-0.5;
zlim3D(2) = header.spacingZ*(header.sizeZ-0.5);

minLim = min([xlim3D(1) ylim3D(1) zlim3D(1)]);
maxLim = max([xlim3D(2) ylim3D(2) zlim3D(2)]);

maxXYLim = max([ylim3D(2)-ylim3D(1) xlim3D(2)-xlim3D(1)]);
maxXZLim = max([zlim3D(2)-zlim3D(1) xlim3D(2)-xlim3D(1)]);
maxYZLim = max([ylim3D(2)-ylim3D(1) zlim3D(2)-zlim3D(1)]);





%% Define constants
topl = 1;
topr = 2;
botl = 3;
botr = 4;

X = 1;
Y = 2;
Z = 3;

%% Set defaults for input arguments
if nargin <1
    dicom_filename = '*';
end

if nargin <2
    fighandle = [];
end

if ~ishandle( fighandle )
    msgstr = sprintf('Warning in %s: Second argument is not a valid figure handle. Creating a new one', mfilename);
    disp( msgstr)
     fighandle = [];
end

if nargin <3
   transparency = 1;
end

if isempty( dicom_filename )
    dicom_filename = '*';
end

%save the current path
orig_path = pwd;
 
% Call the filebrowswer UI if the dicom filename contains a wildcard
if ~isempty( findstr( dicom_filename, '*') )
    [ filename, fpath ] = uigetfile( dicom_filename , 'Select DICOM file');
    cd( fpath );
    dicom_filename =filename;
end


%% Read DICOM header and image data
dinfo= dicominfo( dicom_filename );
imdata = dicomread(  dicom_filename );
cd( orig_path )

%% Calculate slice corner positions from the DICOM header info
% Get the top left corner position in XYZ coordinates
pos    = dinfo.ImagePositionPatient;
 
nc = double(dinfo.Columns);
nr = double(dinfo.Rows);

% Get the row and column direction vercors in XYZ coordinates
row_dircos(X) = dinfo.ImageOrientationPatient(1);
row_dircos(Y) = dinfo.ImageOrientationPatient(2);
row_dircos(Z) = dinfo.ImageOrientationPatient(3);
col_dircos(X) = dinfo.ImageOrientationPatient(4);
col_dircos(Y) = dinfo.ImageOrientationPatient(5);
col_dircos(Z) = dinfo.ImageOrientationPatient(6);

% % Check normality and orthogonality of the row and col vectors
% Crownorm = dot(row_dircos, row_dircos);
% Ccolnorm = dot(col_dircos, col_dircos);
% Cdotprod = dot(row_dircos, col_dircos);
% 
% if abs(Cdotprod) > 1e-5
%     warnstr = sprintf('Possible dicominfo error: the dotproduct of the row and col vectors is %f should be 0',Cdotprod );
%     disp(warnstr)
% end

% Calculate image dimensions
row_length = dinfo.PixelSpacing(1) * nr;
col_length = dinfo.PixelSpacing(2) * nc;


%% Set up the corner positions matrix in XYZ coordinates
% Top Left Hand Corner
corners( topl, X) = pos(X);
corners( topl, Y) = pos(Y);
corners( topl, Z) = pos(Z);

% Top Right Hand Corner
corners( topr, X) = pos(X) + row_dircos(X) * row_length;
corners( topr, Y) = pos(Y) + row_dircos(Y) * row_length;
corners( topr, Z) = pos(Z) + row_dircos(Z) * row_length;

% Bottom Left Hand Corner
corners( botl, X) = pos(X) + col_dircos(X) * col_length;
corners( botl, Y) = pos(Y) + col_dircos(Y) * col_length;
corners( botl, Z) = pos(Z) + col_dircos(Z) * col_length;

% Bottom Right Hand Corner
corners( botr, X) = pos(X) + row_dircos(X) * row_length + col_dircos(X) * col_length;
corners( botr, Y) = pos(Y) + row_dircos(Y) * row_length + col_dircos(Y) * col_length;
corners( botr, Z) = pos(Z) + row_dircos(Z) * row_length + col_dircos(Z) * col_length;

%% Prepare the figure
% Select active figure, set hold on to alow multiple slices in one figure
if isempty(  fighandle )
     fighandle = figure;
     colormap( gray );
     newfig = 1;
else
    newfig = 0;
end

figure( fighandle );
hold on;

%Tidy up the figure
% set aspect ratio
daspect( [1,1,1]);
set( gca, 'color', 'none')

%% Display slice
%  normalize image data
imdata = double( imdata );
imdata = imdata / max( imdata(:));
% scale the image
I = imdata*255;
% create an alternative matrix for corner points
A( 1,1 , 1:3 ) = corners( topl, : );
A( 1,2 , 1:3 ) = corners( topr, : );
A( 2,1 , 1:3 ) = corners( botl, : );
A( 2,2 , 1:3 ) = corners( botr, : );
% extract the coordinates for the surfaces
x = A( :,:,X );
y = A( :,:,Y );
z = A( :,:,Z );

% plot surface
surfh = surface('XData',x,'YData',y,'ZData',z,...
'CData', I,...
'FaceColor','texturemap',...
'EdgeColor','none',...
'LineStyle','none',...
'Marker','none',...
'MarkerFaceColor','none',...
'MarkerEdgeColor','none',...
'CDataMapping','direct');    
%set transparency level
set( surfh, 'FaceAlpha', transparency );
% label axes and optimize figure
xlabel('RL');
ylabel('AP');
zlabel('FH');

% if only one slice is in the figure this may flatten the 3rd axis
if ~newfig
    axis tight
end

return

