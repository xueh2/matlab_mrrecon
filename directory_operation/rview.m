
function rview(target, source, dofin)

if ( nargin == 1 )
    command = ['rview' ' ' target];
end

if ( nargin == 2 )
    command = ['rview' ' ' target ' ' source];
end

if ( nargin == 3 )
    command = ['rview' ' ' target ' ' source ' ' dofin];
end

[s, w] = dos(command, '-echo');

return;