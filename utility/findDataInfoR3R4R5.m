function [R, isStress] = findDataInfoR3R4R5(currDir)

isStress = 1;
if ( (isempty(strfind(currDir, 'STRESS'))==1) & (isempty(strfind(currDir, 'stress'))==1) )
    isStress = 0;
end

R = -1;
if ( isempty(strfind(currDir, 'R3_'))~=1 )        
    R = 3;
end

if ( isempty(strfind(currDir, '_R3_'))~=1 )        
    R = 3;
end

if ( isempty(strfind(currDir, '_3_'))~=1 )        
    R = 3;
end

if ( isempty(strfind(currDir, '_R4_'))~=1 )        
    R = 4;
end

if ( isempty(strfind(currDir, '_4_'))~=1 )        
    R = 4;
end

if ( isempty(strfind(currDir, '_R5_'))~=1 )
    R = 5;
end

if ( isempty(strfind(currDir, '_r5_'))~=1 )
    R = 5;
end

if ( isempty(strfind(currDir, '_5_'))~=1 )        
    R = 5;
end
