
function FormateFileNames_Prefix(home, subdirectory, suffix, prefix)
% find the first hdr file in the home/subdirectory and change its name to
% home_suffix[1..n].hdr
% home is the absolute path and subdirectory is only a subdirectory name
path_abo = fullfile(home, subdirectory, '*.hdr' );
indir = dir(path_abo) ;

num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end

% place = find(home(length(home):-1:1) == filesep );
% prefix = home(length(home)-place(1)+2:length(home));

for i = 1:num
    fullname = fullfile(home, subdirectory, indir(i).name);
    [pathstr,name,ext,versn] = fileparts(fullname);
    newname = [prefix '_' suffix '_' num2str(i) ext];
    newfullname = fullfile(home, subdirectory, newname);
    command = ['rename' ' ' fullname ' ' newname ];
    [s, w] = dos(command);
    
    % img
    name = [name '.img'];
    fullname2 = fullfile(home, subdirectory, name);
    newname = [prefix '_' suffix '_' num2str(i) '.img'];
    command = ['rename' ' ' fullname2 ' ' newname ];
    [s, w] = dos(command);
end

disp('FormateFileNames finished...');

return;