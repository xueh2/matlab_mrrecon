
CreateDirs_inCurrent(home, workingDirs);

filename = fullfile(home, RegResult, ['transformed_CT_AR' suffix2]);
if ( isempty(dir(filename)) )
    return;
end

ImageName = filename;

CreateDirs_inCurrent(fullfile(home, workingDirs), SingleSlice);
Save3DTo2D(fullfile(home, workingDirs), SingleSlice, ImageName, prefix, suffix2)

parametersFileDir = 'C:\huixue\work\include\parameters2';

numOfLevels = 2;
dx = 14;
dy = 14;

parameterFile_FFD = cell(numOfLevels, 1);

parameterFile_FFD{1} = fullfile(parametersFileDir, 'parameters2D-10mm.ffd');
parameterFile_FFD{2} = fullfile(parametersFileDir, 'parameters2D-5mm.ffd');
