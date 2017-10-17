
function inFlag = inCellArray(arrays, item)

inFlag = false;

num = numel(arrays);

if ( isstr(item) )

    for k=1:num

        if ( strcmp(arrays{k}, item) )
            inFlag = true;
            return;
        end

    end

else
    for k=1:num
    
        if ( arrays{k} == item )
            inFlag = true;
            return;
        end
    
    end
end