function h=imsurf(image,upperLeftPoint3,normal,imXDirVec,scale,varargin)
% function h=imsurf(image,upperLeftPoint3,normal,imXDirVec,scale,varargin)
%Plot an image in a 3d plane. The position and scale of the image are controlled by the input
%parameters.
%
%The plane is defined by its normal direction vector (normal), the coordinate of the upper left hand
%corner of the image in the coordinates of the final plot (upperLeftPoint3) and the direction of the
%positive x direction (imXDirVec).
%
%The normal and the imXDirVec have to be orthogonal, if they aren't orthogonal, the function will
%change the imXDirVec vector to make it orthogonal to the normal, in the plane of the two original
%vectors. If they're parallel or zero, the function will fail to draw the image plane.
%
%The scale is defined by scale.
%
%The normal is defined as coming out of the screen in usual viewing.
%
%varargin are paramater value pairs suitable for surface
%
%The colours of the image input are expected to be displayed in the same way that imshow would
%display them.
%
%Author: M Arthington
%Date: 29/08/2010
%Example 1:
% myIm=load('mandrill');
% figure;
% imsurf(myIm.X);
% axis equal
% colormap(gray)
%Example 2:
% myIm=load('mandrill');
% figure;
% hold on;
% imsurf(myIm.X,[],[-1  0 0],[0 -1 0],0.1);
% imsurf(myIm.X,[],[-1 -1 0],[1 -1 0],0.2);
% imsurf(myIm.X,[],[0  -1 0],[1  0 0],0.3);
% axis equal
% view([-35 35]);
% colormap(gray)

if ~(exist('upperLeftPoint3','var') && ~isempty(upperLeftPoint3))
	upperLeftPoint3 = [0 0 0];
end
if ~(exist('normal','var') && ~isempty(normal))
	normal = [0 0 1];
end
if ~(exist('imXDirVec','var') && ~isempty(imXDirVec))
	imXDirVec = [1 0 0];
end
if ~(exist('scale','var') && ~isempty(scale))
	scale = 1;
end

if dot(normal,imXDirVec) ~=0
	information(3,'Making imXDirVec normal to normal');
	imXDirVec = cross(normal,cross(imXDirVec,normal));
end

R=rotationMatrixFromTwoVectors([0 0 1]',-normal);
R2=rotationMatrixFromTwoVectors(R*[1 0 0]',imXDirVec(:));

% x = [0 1; 0 1]*size(image,2);
x = [0 1; 0 1]*size(image,2);
y = [0 0; 1 1]*size(image,1);
z = [0 0; 0 0];

xyz = R2*R*[x(:)';y(:)';z(:)'];
x = xyz(1,:);
x = reshape(x,2,2);
y = xyz(2,:);
y = reshape(y,2,2);
z = xyz(3,:);
z = reshape(z,2,2);

clear xyz
x = x*scale + upperLeftPoint3(1);
y = y*scale + upperLeftPoint3(2);
z = z*scale + upperLeftPoint3(3);

h=surface('XData',x,'YData',y,'ZData',z,'CData',image,'FaceColor','texturemap','EdgeColor','none',varargin{:});
end
function [R,theta]=rotationMatrixFromTwoVectors(vStart,vEnd)
%function [R,theta]=rotationMatrixFromTwoVectors(vStart,vEnd)
%Find the matrix R such that vEnd = R*vStart;
%theta is returned in radians
%Author: M Arthington
%Date: 18/07/2010
%Adapted from: ANGLEAXIS2MATRIX
% Copyright (c) 2008 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/

if numel(vStart) ==2
	vStart(3) = 0;
end
if numel(vEnd) ==2
	vEnd(3) = 0;
end

vStartN = vStart/norm(vStart);
vEndN = vEnd/norm(vEnd);

t = cross(vStartN,vEndN);

theta = acos(dot(vStartN,vEndN));

if all(roundDecPlace(t,5)==0) && theta>eps
	
	%Find an orthogonal vector to the start vector
	t = cross(vStartN,vStartN([2 3 1]));
	while all(~logical(t)) || any(isnan(t))
		t = cross(vStartN,rand(3,1));
	end
	t = theta*t/norm(t);
end
if theta < eps   % If the rotation is very small...
	t = [0 0 0];
	R = [ 1   -t(3) t(2)
		t(3) 1   -t(1)
		-t(2) t(1) 1
		];
	return
end

t = theta*t/norm(t);

% Otherwise set up standard matrix, first setting up some convenience
% variables
t = t/theta;
x = t(1);
y = t(2);
z = t(3);

c = cos(theta); s = sin(theta); C = 1-c;
xs = x*s;   ys = y*s;   zs = z*s;
xC = x*C;   yC = y*C;   zC = z*C;
xyC = x*yC; yzC = y*zC; zxC = z*xC;

R = [ x*xC+c   xyC-zs   zxC+ys
	xyC+zs   y*yC+c   yzC-xs
	zxC-ys   yzC+xs   z*zC+c
	];
end
function [n] = roundDecPlace(input,decPlace)
%function [n] = roundDecPlace(input,decPlace)
%Round a number, input, to a specified number of decimal places, decPlace.
%Works up to 15 significant figures
%This function is matricised (it works on elements of matrices)

tenToDecPlace = 10^decPlace;
n = input.*tenToDecPlace;
n = round(n);

n = n./tenToDecPlace;
end