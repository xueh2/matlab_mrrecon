function data = readGadgetronBatchResult_h5(filename)
% read gadgetron h5 outputs
% all series will be read

d_info = h5info(filename);

num_of_data = numel(d_info.Groups.Groups)

data = cell(num_of_data,4);

for n=1:num_of_data
    
    d_data = h5read(filename, [d_info.Groups.Groups(n).Name '/data']);
    d_attributes = h5read(filename, [d_info.Groups.Groups(n).Name '/attributes']);
    d_header = h5read(filename, [d_info.Groups.Groups(n).Name '/header']);
    
    data{n, 1} = d_info.Groups.Groups(n).Name;
    data{n, 2} = d_data;
    data{n, 3} = d_attributes;
    data{n, 4} = d_header;    
end