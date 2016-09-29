
function tUsed = PerformGadgetronRecon_WholeStudy_List(dataDir)
% PerformGadgetronRecon_WholeStudy_List(dataDir)

[names, num] = findFILE(dataDir, '*Adj*.dat'); 
[names, num_dat] = findFILE(dataDir, '201*.dat'); 

if(num==0 & num_dat>0)
    disp('Auto saved study detected ...');
    [names, num] = findFILE(dataDir, '201*.dat');
    auto_saved = 1;
    hasAdj = 1;
else
    disp('Manually saved study detected ...');   
    auto_saved = 0;    
    hasAdj = 0;    
    [names, num] = findFILE(dataDir, 'meas_*.dat');               
end

for n=1:num
    disp([num2str(n) ' - ' names{n}])
end
