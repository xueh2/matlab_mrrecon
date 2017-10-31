
function TransformationRun_Binary(target, source, output, dof)
% rigidly register the image from t2 to tof

if (isempty(dir(output))==0)
    disp([' exists : ' output ' ... ' ]);
    return;
end

command = ['transformation' ' ' source ' ' output ' ' '-dofin' ' ' dof ' ' '-target' ' ' target ' '  '-bspline' ];
command
[s, w] = dos(command, '-echo');

[data, header] = LoadAnalyze(output, 'Grey');
se = strel('ball', 3, 3, 0);
pp = imclose(data, se);
pp(find(pp>0)) = 1;
SaveAnalyze(uint32(pp), header, output, 'Grey');

return