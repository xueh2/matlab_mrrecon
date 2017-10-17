function flag = isFileExist(filename)

flag = isempty(dir(filename));
flag = ~flag;
