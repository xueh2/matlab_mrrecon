function CY = thumb(varargin)
%THUMB generates thumbnail of an image
%   ThumbnailImage = THUMB(InputImage)
%   ThumbnailImage = THUMB(InputImage, ThumbnailSize)
%   ThumbnailSize must be an integer between 32 and 384
%    Member of IMUI
%    Kotaimen.C, 2002/05 - 2002/07, All Rights Reserved

error(nargchk(1, 2, nargin))

switch nargin
case 1
	CX = varargin{1};
	% set default thumbnail size
	Tsize = 128;
case 2
	CX = varargin{1};
	Tsize = round(varargin{2});
	if Tsize < 32 | Tsize > 384
		error('ThumbnailSize must between 32 and 384.')
	end
end

Xsize = size(CX);
% Resize to thumbnail size use IMRESIZE
CY = imresize(CX, Tsize / max(Xsize), 'nearest', 11);
