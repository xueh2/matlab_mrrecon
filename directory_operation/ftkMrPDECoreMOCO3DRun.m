
function ftkMrPDECoreMOCO3DRun(target, source, transform, rigid, nonRigid, inv, dofIn, dofOut, dxFile, dyFile, dzFile, ...
                                    dxInvFile, dyInvFile, dzInvFile, iters, sigma, neighbor, stepDiv, algoRigid, algo, interp, rigidInterp, BatFileName, exePath)
% MrPDE core rigid/nonRigid

if ( isempty(exePath) )    
    command = ['ftkMrPDECoreMOCO3D' '  ' target '  ' source ];
else
    command = [exePath 'ftkMrPDECoreMOCO3D' '  ' target '  ' source ];
end

if ( isempty(transform) == 0 )
    command = [command '  ' '-transform' '  ' transform ];
end

if ( isempty(rigid) == 0 )
    command = [command '  ' '-rigid' '  ' rigid ];
end

if ( isempty(nonRigid) == 0 )
    command = [command '  ' '-nonRigid' '  ' nonRigid ];
end

if ( isempty(inv) == 0 )
    command = [command '  ' '-inv' '  ' inv ];
end

if ( isempty(dofIn) == 0 )
    command = [command '  ' '-dofIn' '  ' dofIn ];
end

if ( isempty(dofOut) == 0 )
    command = [command '  ' '-dofOut' '  ' dofOut ];
end

if ( (isempty(dxFile)==0) & (isempty(dyFile)==0) & (isempty(dzFile)==0) )
    command = [command '  ' '-deform' '  ' dxFile '  ' dyFile '  ' dzFile ];
end

if ( (isempty(dxInvFile)==0) & (isempty(dyInvFile)==0) & (isempty(dzInvFile)==0) )
    command = [command '  ' '-deformInv' '  ' dxInvFile '  ' dyInvFile '  ' dzInvFile ];
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

if ( isempty(algoRigid) == 0 )
    command = [command '  ' '-algoRigid' '  ' algoRigid ];
end

if ( isempty(algo) == 0 )
    command = [command '  ' '-algo' '  ' algo ];
end

if ( isempty(interp) == 0 )
    command = [command '  ' '-interp' '  ' interp ];
end

if ( isempty(rigidInterp) == 0 )
    command = [command '  ' '-rigidInterp' '  ' rigidInterp ];
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