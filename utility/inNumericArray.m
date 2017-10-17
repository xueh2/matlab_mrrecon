
function inFlag = inNumericArray(arrays, item)

inFlag = false;

num = numel(arrays);

for k=1:num

    if ( arrays(k) == item )
        inFlag = true;
        return;
    end

end