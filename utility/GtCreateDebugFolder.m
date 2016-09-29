function GtCreateDebugFolder(debugFolder, ut)

if nargin<2
    ut = 'D:\gtuser\mrprogs\install'
end

if nargin<1
    debugFolder = 'DebugOutput'
end

command = ['rmdir /S /Q ' fullfile(ut, debugFolder)];
dos(command);

mkdir(fullfile(ut, debugFolder));
mkdir(fullfile(ut, debugFolder, '0'));
mkdir(fullfile(ut, debugFolder, '1'));
mkdir(fullfile(ut, debugFolder, '2'));
mkdir(fullfile(ut, debugFolder, '3'));