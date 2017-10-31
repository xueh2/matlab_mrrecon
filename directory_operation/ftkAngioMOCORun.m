
function ftkAngioMOCORun(pre, post, frameOfReferenceSame, regStrategy, rigid, stepDivRigid, histogramSize, subsamplings, rigidInitialFlag, nonRigid, invFlag, iters, sigma, neighbor, stepDiv, algo, numOfRep, adaptiveSigmaForNonRigid, minimalSigma, intensityClamp, volumePreserving, BatFileName, exePath)
% MrPDE core rigid/nonRigid

if ( isempty(exePath) )    
    command = ['ftkAngioMOCO' '  ' ];
else
    command = [exePath 'ftkAngioMOCO' '  ' ];       
end

N = numel(pre);
M = numel(post);

command = [command ' ' num2str(N) ' ' num2str(M) ' '];

for kk=1:N    
    command = [command ' ' pre{kk} ' '];        
end

for kk=1:M
    command = [command ' ' post{kk} ' '];        
end

if ( frameOfReferenceSame > 0 )
    command = [command ' ' '-frameOfReferenceSame true'];
else
    command = [command ' ' '-frameOfReferenceSame false'];
end

if ( isempty(regStrategy) == 0 )
    command = [command '  ' '-regStrategy' '  ' regStrategy ];
end

if ( rigid > 0 )
    command = [command ' ' '-rigid true'];
else
    command = [command ' ' '-rigid false'];
end

if ( isempty(stepDivRigid) == 0 )
    command = [command '  ' '-stepDivRigid' '  ' num2str(stepDivRigid) ];
end

if ( isempty(subsamplings) == 0 )
    command = [command '  ' '-subsamplings' ' ' num2str(numel(subsamplings)) ' ' num2str(subsamplings) ];
end

if ( rigidInitialFlag > 0 )
    command = [command ' ' '-rigidInitialFlag true'];
else
    command = [command ' ' '-rigidInitialFlag false'];
end

if ( nonRigid > 0 )
    command = [command ' ' '-nonRigid true'];
else
    command = [command ' ' '-nonRigid false'];
end

if ( invFlag > 0 )
    command = [command ' ' '-invFlag true'];
else
    command = [command ' ' '-invFlag false'];
end

if ( isempty(iters) == 0 )
    command = [command '  ' '-iters' '  ' num2str(numel(iters)) ' ' num2str(iters) ];
end

if ( isempty(sigma) == 0 )
    command = [command '  ' '-sigma' '  ' num2str(sigma) ];
end

if ( isempty(neighbor) == 0 )
    command = [command '  ' '-neighbor' '  ' num2str(neighbor) ];
end

if ( isempty(stepDiv) == 0 )
    command = [command '  ' '-stepDiv' '  ' num2str(stepDiv) ];
end

if ( isempty(algo) == 0 )
    command = [command '  ' '-algo' '  ' algo ];
end

if ( isempty(numOfRep) == 0 )
    command = [command '  ' '-numOfRep' '  ' num2str(numOfRep) ];
end

if ( adaptiveSigmaForNonRigid > 0 )
    command = [command ' ' '-adaptiveSigmaForNonRigid true'];
else
    command = [command ' ' '-adaptiveSigmaForNonRigid false'];
end

if ( isempty(minimalSigma) == 0 )
    command = [command '  ' '-minimalSigma' '  ' num2str(minimalSigma) ];
end

if ( intensityClamp > 0 )
    command = [command ' ' '-intensityClamp true'];
else
    command = [command ' ' '-intensityClamp false'];
end

if ( volumePreserving > 0 )
    command = [command ' ' '-volumePreserving true'];
else
    command = [command ' ' '-volumePreserving false'];
end

disp(command)

if ( isempty(BatFileName) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end
return
