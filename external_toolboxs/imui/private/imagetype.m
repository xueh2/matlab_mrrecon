function imgType = imagetype(CX)
% imagetype -returns the type of image CX
%    Member of IMUI
%    Kotaimen.C, 2002/05 - 2002/07, All Rights Reserved

if isbw(CX)
	imgType = 'Binary';
	return
end
if isgray(CX)
	imgType = 'Gray';
	return
end
if isrgb(CX)
	imgType = 'RGB';
	return
end
imgType = ''; %otherwise return empty string