
[names, num] = findFILE(home, suffix1);
if ( num < 1 )
    ImageName = [];
    return;
end
ImageName = names{1}

if ( isempty(dir(fullfile(home, SingleSlice))) )
    CreateDirs_inCurrent(home, SingleSlice);
end

if ( isempty(dir(fullfile(home, workingDirs))) )
    CreateDirs_inCurrent(home, workingDirs);
end

if ( isempty(dir(fullfile(home, RegResult))) )
    CreateDirs_inCurrent(home, RegResult);
end

if ( isempty(dir(fullfile(home, ROI))) )
    CreateDirs_inCurrent(home, ROI);
end


[pathstr, name, ext, versn] = fileparts(ImageName);

[patientName, serialNo, slicePosition, sliceNum] = parsePerfusionSliceWithSliceNum([name ext]);
hdrnames = dir(fullfile(home, SingleSlice, '*.hdr'));
numHdr = numel(hdrnames);
if ( (numHdr==0) | (numHdr ~= sliceNum) )
    if ( exist('minW') )
        Save3DTo2D_Resize(home, SingleSlice, ImageName, prefix, suffix2, minW, minH);
    else
        Save3DTo2D(home, SingleSlice, ImageName, prefix, suffix2);
    end
end

[names, num] = findFILE(home, suffix1);
if ( num < 1 )
    ImageName = [];
    return;
end
ImageName = names{1}

filename = fullfile(home, 'RegistrationRecord.txt');
info = ReadRegistrationRecord(filename);
templateNo = info.templateNo;
leftup = info.leftup;
rightdown = info.rightdown;


parametersFileDir_AR = 'C:\hxue\include\parameters2';

parameterFile_Rigid2D = fullfile(parametersFileDir_AR, 'parameters2D.rigid');
parameterFile_Affine2D = fullfile(parametersFileDir_AR, 'parameters2D.affine');
parameterFile_Rigid3D = fullfile(parametersFileDir_AR, 'parameters3D.rigid');
parameterFile_Affine3D = fullfile(parametersFileDir_AR, 'parameters3D.affine');
