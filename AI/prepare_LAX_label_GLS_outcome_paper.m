function labels = prepare_LAX_label_GLS_outcome_paper(lax_label)
% labels = prepare_LAX_label_GLS_outcome_paper(lax_label)
% if no landmarks for a image, landmarks is -1

    fields = fieldnames(lax_label)
    max_N=numel(fields);

    labels = [];
    for k=1:max_N
        v = getfield(lax_label, fields{k});
        [v_path, v_name, ext] = fileparts(v.filename);

        disp(['load ' num2str(k) ' - ' v_name]);

        ind = find(v_name=='_');
        pid = str2num(v_name(1:ind(1)-1));
        phs = str2num(v_name(ind(end)+1:end));
        series = v_name(ind(2)+1:ind(end)-1);
        
        pt1 = [-1 -1];
        pt2 = [-1 -1];
        pt3 = [-1 -1];
        start_e1 = 0;
        start_ro = 0;

        if(numel(v.regions)==3)
            pt1 = [v.regions(1).shape_attributes.cx v.regions(1).shape_attributes.cy];
            pt2 = [v.regions(2).shape_attributes.cx v.regions(2).shape_attributes.cy];
            pt3 = [v.regions(3).shape_attributes.cx v.regions(3).shape_attributes.cy];

            if(start_e1>0)
                pt1(1) = pt1(1) - start_e1 -1;
                pt2(1) = pt2(1) - start_e1 -1;
                pt3(1) = pt3(1) - start_e1 -1;
            end

            if(start_ro>0)
                pt1(2) = pt1(2) - start_ro -1;
                pt2(2) = pt2(2) - start_ro -1;
                pt3(2) = pt3(2) - start_ro -1;
            end
        end

        labels = [labels; {v_name}, pid, {series}, phs, pt1, pt2, pt3];
    end
end
