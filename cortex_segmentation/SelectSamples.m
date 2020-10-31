
function [selected_x, selected_label] = SelectSamples(x_I, x_label, post, minimalP, typeofselection, maxN)
% select samples for the kernel density estimation

x_post = post(x_label(:));

if ( strcmp(typeofselection, 'probability') == 1 )
    
    index = find(x_post>=minimalP);
    while ( length(index)>maxN )
        minimalP = minimalP + 0.0005;
        if ( minimalP >= 1.0 )
            break;
        end
        index = find(x_post>=minimalP);
    end
    
    while ( length(index)<maxN/2 )
        minimalP = minimalP - 0.0005;
        if ( minimalP <= 0.02 )
            break;
        end
        index = find(x_post>=minimalP);
    end
    
    if ( length(index)>maxN  )
        % if this happans, we need turn to the uniform sampling
        num = length(index);
        if ( num>maxN )
            step = ceil(num/maxN);
        else
            step = 1;
        end
        selected_x = x_I(index(1:step:end));
        selected_label = x_label(index(1:step:end));
    else    
        selected_x = x_I(index(:));
        selected_label = x_label(index(:));
    end
    
end

if ( strcmp(typeofselection, 'uniform') == 1 )
    
    num = size(x_I, 1);
    
    if ( num>maxN )
        step = ceil(num/maxN);
    else
        step = 1;
    end
    
    selected_x = x_I(1:step:end);
    selected_label = x_label(1:step:end);
    
end

return