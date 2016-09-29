
[filename1, pathname1] = uigetfile('../Projects/B1Mapping/Data/Data20130321/SNV/*.IMA');
[filename2, pathname2] = uigetfile('../Projects/B1Mapping/Data/Data20130321/SNV/*.IMA');

name1=fullfile(pathname1,filename1)
name2=fullfile(pathname2,filename2)

M1 = dicomread(name1);
% info1 = dicominfo(name1);

% info1.Width;
% info1.Height;
% figure(1), imshow(M1,[0,4095])

M2 = dicomread(name2);
% info2 = dicominfo(name2);
% imshow(M2,[0,4095])

dM1 = double(M1);
dM2 = double(M2);
Angle = (180 / pi) * abs (acos(dM2 ./ (2* dM1)));
figure(3), imshow(Angle,[0,100])





