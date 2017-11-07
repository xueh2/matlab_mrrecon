
function Resample_RealRun(home, subdirectory, outputDir, filename, resolution)
% resample the image

% path_abo = fullfile(home, subdirectory, '*.hdr' );
% indir = dir(path_abo) ;
% 
% num = length(indir);
% if ( num == 0 )
%     disp('empty directory');
%     return;
% end
% 
% for i = 1:num
fullname = fullfile(home, subdirectory, filename);
[pathstr,name,ext,versn] = fileparts(fullname);
% [data, header] = LoadAnalyze(fullname, 'Grey');
%     resolution(1) = floor(10*min([header.xvoxelsize header.yvoxelsize header.zvoxelsize]))/10;

if ( resolution(1) < 0.25 )
    resolution(1) = 0.25;
end
resolution(2) = resolution(1);
resolution(3) = resolution(2);

newname = [name '_' num2str(resolution(1)) ext];
pPositions = find(newname == '.');
newname(pPositions(1)) = 'p';
newfullname = fullfile(home, outputDir, newname);

command = ['resample_Real' ' ' fullname ' ' newfullname ' ' '-size ' num2str(resolution(1)) ' ' num2str(resolution(2)) ' ' num2str(resolution(3)) ' ' '-bspline'];
command
[s, w] = dos(command, '-echo');

% end
disp('ResampleRun finished...');
return;