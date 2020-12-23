function csv_polygon = generate_VIA_polygon_str(curr_C, file_name, info)
% csv_polygon = generate_VIA_polygon_str(curr_C, file_name, info)
% info is the dir(file_name) for file size in bytes

csv_str = [file_name ',' num2str(info.bytes) ',' '"{}",'];

if(~isempty(curr_C))
    csv_str=[csv_str '1,0,'];
    csv_str=[csv_str '"{""name"":""polygon"",""all_points_x"":['];
    N = size(curr_C,1);
    for pt=1:N
        X = round(curr_C(pt,1)-1);
        if(pt<N)
            csv_str=[csv_str num2str(X) ','];
        else
            csv_str=[csv_str num2str(X) '],""all_points_y"":['];
        end
    end
    for pt=1:N
        Y = round(curr_C(pt,2)-1);
        if(pt<N)
            csv_str=[csv_str num2str(Y) ','];
        else
            csv_str=[csv_str num2str(Y) ']}","{}"'];
        end
    end
else
    csv_str=[csv_str '0,0,"{}","{}"'];
end

csv_polygon = csv_str;